# processor-cleanet entry point: Rscript /app/R/main.R
# Reads $INPUT_DIR/{fcs,metadata}, writes $OUTPUT_DIR/{fcs,qc}.

app_dir <- function() {
  if (nzchar(Sys.getenv("APP_ROOT"))) return(file.path(Sys.getenv("APP_ROOT"), "R"))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "R"
}
source(file.path(app_dir(), "common.R"))
source(file.path(app_dir(), "cleanet.R"))

rel_path <- function(paths, root) {
  root <- paste0(normalizePath(root), "/")
  p <- normalizePath(paths)
  ifelse(startsWith(p, root), substring(p, nchar(root) + 1), p)
}

run <- function(input_dir = Sys.getenv("INPUT_DIR", "/data/input"),
                output_dir = Sys.getenv("OUTPUT_DIR", "/data/output")) {
  options(bitmapType = "cairo", warn = 1)
  started <- Sys.time()
  p <- cleanet_params()
  log_line("processor-cleanet %s (Cleanet %s @ %s): THRESH=%d COFACTOR=%s KEEP_DOUBLETS=%s DEBRIS_SUM_MAX=%g THREADS=%d",
           app_version(), as.character(utils::packageVersion("Cleanet")), cleanet_sha(), p$THRESH,
           as.character(p$COFACTOR), p$KEEP_DOUBLETS, p$DEBRIS_SUM_MAX, p$THREADS)

  files <- find_fcs(input_dir)
  if (!length(files)) {
    log_line("ERROR: no .fcs files under %s", input_dir)
    return(1L)
  }
  panel <- read_panel(input_dir)
  if (is.null(panel)) log_line("note: no panel.csv; checking channels with a $PnS (not scatter/time)")
  else log_line("panel: %s (%d marker channels)", rel_path(attr(panel, "path"), input_dir),
                sum(panel$type == "marker"))

  # Two files with the same name would overwrite each other: keep the first.
  dup <- duplicated(basename(files))
  for (f in files[dup]) log_line("%s: SKIPPED, another input file has the same name", f)
  files <- files[!dup]

  fcs_dir <- file.path(output_dir, "fcs")
  qc_dir <- file.path(output_dir, "qc")
  plot_dir <- file.path(qc_dir, "plots")
  for (d in c(fcs_dir, qc_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  log_line("%d file(s), %d at a time", length(files), min(p$THREADS, length(files)))
  one <- function(f) {
    r <- cleanet_sample(f, panel, p, fcs_dir, plot_dir)
    log_line("%s: %s in=%s doublets=%s (%s%%) out=%s sensitivity=%s %.1fs%s",
             r$sample_id, r$status, r$events_in, r$doublets, r$pct_doublets, r$events_out,
             r$sensitivity, r$seconds,
             if (nzchar(r$notes)) paste0(" [", r$notes, "]") else "")
    r
  }
  rows <- if (p$THREADS > 1 && length(files) > 1)
    parallel::mclapply(files, one, mc.cores = p$THREADS, mc.preschedule = FALSE)
  else lapply(files, one)
  # A worker that died (e.g. killed for memory) comes back as an error or NULL.
  rows <- Map(function(r, f) {
    if (!is.list(r) || is.null(r$status)) {
      log_line("%s: failed, worker died: %s", sample_id_of(f), paste(as.character(r), collapse = " "))
      list(sample_id = sample_id_of(f), status = "failed",
           notes = paste("worker died:", paste(as.character(r), collapse = " ")))
    } else r
  }, rows, files)

  cols <- c("sample_id", "events_in", "doublets", "pct_doublets", "events_out", "sensitivity",
            "debris_excluded", "n_channels", "status", "notes")
  tab <- do.call(rbind, lapply(rows, function(r) {
    as.data.frame(lapply(stats::setNames(cols, cols),
                         function(c) if (is.null(r[[c]]) || length(r[[c]]) == 0) NA else r[[c]]),
                  stringsAsFactors = FALSE)
  }))
  ok <- tab$status != "failed"
  # Failed samples have no output FCS; drop anything half-written.
  for (s in tab$sample_id[!ok]) unlink(file.path(fcs_dir, paste0(s, c(".fcs", ".fcs.part"))))
  utils::write.csv(tab, file.path(qc_dir, "doublets.csv"), row.names = FALSE, na = "")

  params <- p
  params$COFACTOR <- as.character(p$COFACTOR)
  prov <- list(
    app = "processor-cleanet", version = app_version(),
    tools = c(pkg_versions(c("Cleanet", "flowCore", "RcppHNSW", "ggplot2")),
              list(Cleanet_commit = cleanet_sha())),
    parameters = params,
    inputs = list(fcs = rel_path(files, input_dir),
                  panel = if (is.null(panel)) NULL else rel_path(attr(panel, "path"), input_dir)),
    started = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    samples = lapply(rows, function(r) r[intersect(c("sample_id", "status", "events_in", "events_out",
                                                     "doublets", "pct_doublets", "seconds", "notes"),
                                                   names(r))])
  )
  for (d in c(fcs_dir, qc_dir)) write_json_file(prov, file.path(d, "provenance.json"))

  log_line("done: %d ok, %d passed through, %d failed; %d of %d events called doublets (%.2f%%); %.0fs",
           sum(tab$status == "ok"), sum(tab$status == "passed_through"), sum(!ok),
           sum(tab$doublets, na.rm = TRUE), sum(tab$events_in, na.rm = TRUE),
           100 * sum(tab$doublets, na.rm = TRUE) / max(1, sum(tab$events_in[ok], na.rm = TRUE)), prov$seconds)
  if (!any(ok)) {
    log_line("ERROR: every sample failed")
    return(1L)
  }
  0L
}

if (sys.nframe() == 0L) quit(status = run(), save = "no")
