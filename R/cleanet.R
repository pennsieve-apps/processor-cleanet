# Cleanet on one sample: back to linear values, simulate doublets, flag events whose
# neighbourhood is dominated by simulated doublets, write the singlets.

suppressPackageStartupMessages({
  library(Cleanet)
  library(ggplot2)
})

cleanet_params <- function() {
  cof <- env_string("COFACTOR", "panel")
  if (tolower(cof) != "panel") {
    v <- suppressWarnings(as.numeric(cof))
    if (is.na(v) || v <= 0) stop(sprintf("COFACTOR must be 'panel' or a positive number, got '%s'", cof))
    cof <- v
  } else cof <- "panel"
  list(
    KEEP_DOUBLETS = env_bool("KEEP_DOUBLETS", FALSE),
    THRESH = env_number("THRESH", 5, min = 0, max = 14, integer = TRUE),
    COFACTOR = cof,
    DEBRIS_SUM_MAX = env_number("DEBRIS_SUM_MAX", 0, min = 0),
    MIN_EVENTS = env_number("MIN_EVENTS", 1000, min = 10, integer = TRUE),
    PLOTS = env_bool("PLOTS", TRUE),
    THREADS = env_number("THREADS", 4, min = 1, max = 64, integer = TRUE),
    SEED = env_number("SEED", 42, integer = TRUE)
  )
}

# Per-channel cofactors of the marker channels from panel.csv (NA = not transformed).
marker_cofactors <- function(panel, channels) {
  if (is.null(panel) || !("cofactor" %in% names(panel))) return(stats::setNames(rep(NA_real_, length(channels)), channels))
  cf <- suppressWarnings(as.numeric(panel$cofactor[match(channels, panel$channel)]))
  stats::setNames(cf, channels)
}

# fcs-prep also records the cofactors in the file: PENNSIEVE_COFACTORS=<channel>=<c>,...
cofactors_from_keyword <- function(ff) {
  kw <- flowCore::keyword(ff)
  i <- which(toupper(names(kw)) == "PENNSIEVE_COFACTORS")  # keywords are case-insensitive
  if (!length(i)) return(numeric())
  v <- as.character(kw[[i[1]]])
  if (!nzchar(v)) return(numeric())
  parts <- strsplit(strsplit(v, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
  parts <- parts[lengths(parts) == 2]
  stats::setNames(suppressWarnings(as.numeric(vapply(parts, `[`, "", 2))), vapply(parts, `[`, "", 1))
}

# Cleanet adds two events' values to simulate a doublet, so it needs linear values, and
# arcsinh(x / cofactor)-transforms internally (one cofactor for all columns). fcs-prep
# already transformed the markers, each with its panel cofactor c_j. Undoing that gives
# raw_j = c_j * sinh(x_j). With COFACTOR=panel the data handed to Cleanet are
# raw_j / c_j = sinh(x_j) with cofactor 1, so its distances are computed on exactly the
# arcsinh(raw_j / c_j) scale the rest of the workflow uses; a numeric COFACTOR instead
# hands raw_j with that one cofactor. Markers without a panel cofactor are taken as
# linear already (divided by 6000 under COFACTOR=panel).
linear_matrix <- function(x, cofactors, cofactor_param) {
  out <- x
  for (j in colnames(x)) {
    cf <- cofactors[[j]]
    raw_over_c <- if (is.na(cf)) x[, j] / 6000 else sinh(x[, j])
    raw <- if (is.na(cf)) x[, j] else sinh(x[, j]) * cf
    out[, j] <- if (identical(cofactor_param, "panel")) raw_over_c else raw
  }
  out
}

diagnostic_plot <- function(path, sid, res, ff, keep_debris, p) {
  nb <- data.frame(
    n = c(res$n_nb_doub, res$n_nb_sim),
    kind = rep(c("events in the file", "simulated doublets"),
               c(length(res$n_nb_doub), length(res$n_nb_sim))))
  g1 <- ggplot(nb, aes(x = .data$n, fill = .data$kind)) +
    geom_bar(aes(y = after_stat(prop), group = .data$kind), position = "dodge", width = 0.8) +
    geom_vline(xintercept = p$THRESH + 0.5, linetype = "dashed") +
    scale_fill_manual(values = c("events in the file" = "#4C78A8", "simulated doublets" = "#E45756"), name = NULL) +
    scale_x_continuous(breaks = 0:15) +
    labs(x = "simulated doublets among the 15 nearest neighbours",
         y = "fraction of events",
         title = sprintf("%s: %d doublets (%.2f%%)", sid, sum(res$status == "Doublet"),
                         100 * mean(res$status == "Doublet")),
         subtitle = sprintf("%.1f%% of simulated doublets recognised; dashed line: threshold (> %d)",
                            100 * res$sensitivity, p$THRESH)) +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  plots <- list(g1)
  e <- flowCore::exprs(ff)
  if (all(c("FSC-A", "FSC-H") %in% colnames(e))) {
    set.seed(p$SEED)
    singlets <- which(res$status != "Doublet")
    doublets <- which(res$status == "Doublet")
    idx <- c(sample(singlets, min(length(singlets), 20000)), doublets)
    df <- data.frame(a = e[idx, "FSC-A"], h = e[idx, "FSC-H"],
                     status = factor(res$status[idx], levels = c("Singlet", "Doublet")))
    df <- df[order(df$status), ]
    g2 <- ggplot(df, aes(x = .data$a, y = .data$h, colour = .data$status)) +
      geom_point(size = 0.3, alpha = 0.5) +
      scale_colour_manual(values = c(Singlet = "grey55", Doublet = "#E45756"), name = NULL) +
      labs(x = "FSC-A", y = "FSC-H", title = "Cleanet calls on forward scatter (not used by Cleanet)") +
      theme_bw(base_size = 11) + theme(legend.position = "top") +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
    plots <- c(plots, list(g2))
  }
  grDevices::png(path, width = 6 * length(plots), height = 5, units = "in", res = 110, type = "cairo")
  on.exit(grDevices::dev.off())
  if (length(plots) == 1) print(plots[[1]]) else gridExtra::grid.arrange(grobs = plots, ncol = length(plots))
}

# Process one FCS file. Returns a one-row list for doublets.csv. Never throws.
cleanet_sample <- function(path, panel, p, fcs_dir, plot_dir) {
  sid <- sample_id_of(path)
  t0 <- proc.time()[["elapsed"]]
  notes <- character()
  row <- list(sample_id = sid, events_in = NA_integer_, doublets = NA_integer_,
              pct_doublets = NA_real_, events_out = NA_integer_, sensitivity = NA_real_,
              debris_excluded = NA_integer_, n_channels = NA_integer_, status = "failed", notes = "")
  res <- tryCatch(withCallingHandlers({
    ff <- read_fcs(path)
    n_in <- nrow(ff)
    row$events_in <- n_in
    markers <- panel_channels(ff, panel, "marker")
    if (length(markers) < 2) stop("fewer than 2 marker channels (panel.csv type == marker) in this file")
    row$n_channels <- length(markers)
    if (n_in < p$MIN_EVENTS) {
      notes <- c(notes, sprintf("passed through unchanged: only %d events (< MIN_EVENTS %d)", n_in, p$MIN_EVENTS))
      write_fcs(ff, file.path(fcs_dir, basename(path)),
                list(PENNSIEVE_DOUBLETS = "cleanet;removed=0;passed_through=true"))
      row$doublets <- 0L
      row$pct_doublets <- 0
      row$events_out <- n_in
      row$status <- "passed_through"
      row
    } else {
      cofactors <- marker_cofactors(panel, markers)
      kw_cf <- cofactors_from_keyword(ff)
      fill <- is.na(cofactors) & markers %in% names(kw_cf)
      cofactors[fill] <- kw_cf[markers[fill]]
      if (identical(p$COFACTOR, "panel") && any(is.na(cofactors)))
        notes <- c(notes, sprintf("no panel cofactor for %s: taken as linear, scaled by 1/6000",
                                  paste(markers[is.na(cofactors)], collapse = ", ")))
      x <- linear_matrix(flowCore::exprs(ff)[, markers, drop = FALSE], cofactors, p$COFACTOR)
      df <- as.data.frame(x, check.names = FALSE)
      is_debris <- NULL
      if (p$DEBRIS_SUM_MAX > 0) {
        e <- flowCore::exprs(ff)
        if (all(c("FSC-A", "SSC-A") %in% colnames(e))) {
          # Cleanet::filter_debris_flow's rule, without its interactive plot.
          is_debris <- e[, "FSC-A"] + e[, "SSC-A"] < p$DEBRIS_SUM_MAX & e[, "FSC-A"] > 0 & e[, "SSC-A"] > 0
          row$debris_excluded <- sum(is_debris)
        } else notes <- c(notes, "DEBRIS_SUM_MAX set but no FSC-A/SSC-A: no debris excluded")
      }
      set.seed(p$SEED)
      res <- suppressMessages(Cleanet::cleanet(df, cols = markers,
                                               cofactor = if (identical(p$COFACTOR, "panel")) 1 else p$COFACTOR,
                                               thresh = p$THRESH, is_debris = is_debris))
      rm(df, x)
      is_doublet <- res$status == "Doublet"
      n_doub <- sum(is_doublet)
      row$doublets <- n_doub
      row$pct_doublets <- round(100 * n_doub / n_in, 3)
      row$sensitivity <- round(res$sensitivity, 4)
      keep <- if (p$KEEP_DOUBLETS) rep(TRUE, n_in) else !is_doublet
      row$events_out <- sum(keep)
      kw <- sprintf("cleanet;doublets=%d;pct_doublets=%.2f;removed=%d;thresh=%d;cofactor=%s",
                    n_doub, 100 * n_doub / n_in, n_in - sum(keep), p$THRESH, as.character(p$COFACTOR))
      write_fcs(ff[keep, ], file.path(fcs_dir, basename(path)), list(PENNSIEVE_DOUBLETS = kw))
      if (p$PLOTS) tryCatch(
        diagnostic_plot(file.path(plot_dir, paste0(sid, ".png")), sid, res, ff, is_debris, p),
        error = function(e) notes <<- c(notes, paste("plot failed:", conditionMessage(e))))
      row$status <- "ok"
      row
    }
  }, warning = function(w) {
    notes <<- c(notes, gsub("\\s+", " ", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    notes <<- c(notes, conditionMessage(e))
    row
  })
  res$notes <- paste(unique(notes), collapse = "; ")
  res$seconds <- round(proc.time()[["elapsed"]] - t0, 1)
  res
}

# The Cleanet commit the image was built from (remotes records it).
cleanet_sha <- function() {
  sha <- utils::packageDescription("Cleanet")$RemoteSha
  if (is.null(sha)) "unknown" else sha
}
