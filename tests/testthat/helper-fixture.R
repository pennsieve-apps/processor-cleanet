# A tiny synthetic fixture in the fcs-prep contract, built in R so the tests need no data.
app_root <- Sys.getenv("APP_ROOT", normalizePath(file.path("..", "..")))
source(file.path(app_root, "R", "common.R"))
source(file.path(app_root, "R", "cleanet.R"))
source(file.path(app_root, "R", "main.R"))

markers <- c(`BV421-A` = "CD3", `PE-A` = "CD4", `APC-A` = "CD8", `FITC-A` = "CD19",
             `PE-Cy7-A` = "CD56", `BV605-A` = "CD16")
COF <- 6000

# n singlets from four populations plus n_doublets made by adding two singlets' linear
# values (what a doublet is), then arcsinh(x / 6000) on the markers as fcs-prep does.
# The doublets are flagged in the returned attribute "is_doublet".
make_frame <- function(n, n_doublets = 0, seed = 1) {
  set.seed(seed)
  centres <- rbind(c(4e4, 3e4, 0, 0, 0, 0),      # CD4 T
                   c(4e4, 0, 5e4, 0, 0, 0),      # CD8 T
                   c(0, 0, 0, 6e4, 0, 0),        # B
                   c(0, 0, 0, 0, 3e4, 5e4))      # NK
  pop <- sample(1:4, n, replace = TRUE)
  lin <- centres[pop, ] * exp(matrix(rnorm(n * 6, sd = 0.15), n)) + matrix(rnorm(n * 6, sd = 400), n)
  fsc <- rnorm(n, 8e5, 8e4)
  fsch <- fsc * rnorm(n, 0.8, 0.02)
  if (n_doublets > 0) {
    a <- sample(n, n_doublets, replace = TRUE)
    b <- sample(n, n_doublets, replace = TRUE)
    lin <- rbind(lin, lin[a, ] + lin[b, ])
    fsc <- c(fsc, fsc[a] + fsc[b])
    fsch <- c(fsch, (fsch[a] + fsch[b]) * 0.6)
  }
  N <- nrow(lin)
  ord <- sample(N)  # doublets arrive mixed in with the singlets
  x <- cbind(`FSC-A` = fsc, `FSC-H` = fsch, `SSC-A` = rnorm(N, 3e5, 5e4),
             asinh(lin / COF), Time = sort(runif(N, 0, 1e5)))[ord, ]
  x[, "Time"] <- sort(x[, "Time"])
  colnames(x)[4:9] <- names(markers)
  ff <- flowCore::flowFrame(x)
  pd <- flowCore::pData(flowCore::parameters(ff))
  pd$desc <- c(NA, NA, NA, unname(markers), NA)
  flowCore::pData(flowCore::parameters(ff)) <- pd
  kw <- flowCore::keyword(ff)
  kw[["$DATE"]] <- "01-Jan-2025"
  kw[["PENNSIEVE_PREP"]] <- "transform=arcsinh;cofactor=6000"
  kw[["PENNSIEVE_QC"]] <- "peacoqc;removed=0"
  ff@description <- kw
  attr(ff, "is_doublet") <- (seq_len(N) > n)[ord]
  ff
}

write_panel <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- data.frame(channel = c("FSC-A", "FSC-H", "SSC-A", names(markers), "Time"),
                  marker = c("FSC-A", "FSC-H", "SSC-A", unname(markers), "Time"),
                  type = c("scatter", "scatter", "scatter", rep("marker", 6), "time"),
                  cofactor = c("", "", "", rep(COF, 6), ""),
                  clustering = c("false", "false", "false", rep("true", 6), "false"),
                  state = "false")
  utils::write.csv(p, file.path(dir, "panel.csv"), row.names = FALSE, quote = FALSE)
}

make_input <- function(layout = "linked") {
  root <- tempfile("in")
  fcs <- if (layout == "linked") file.path(root, "fcs") else root
  meta <- if (layout == "linked") file.path(root, "metadata") else root
  dir.create(fcs, recursive = TRUE)
  write_fcs(make_frame(8000, n_doublets = 400, seed = 1), file.path(fcs, "s01.fcs"))
  write_fcs(make_frame(6000, n_doublets = 0, seed = 2), file.path(fcs, "s02.fcs"))
  write_fcs(make_frame(300, seed = 3), file.path(fcs, "tiny.fcs"))
  write_panel(meta)
  root
}
