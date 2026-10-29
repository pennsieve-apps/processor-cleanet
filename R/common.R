# Shared helpers: parameters, input discovery, FCS I/O, provenance.
# (The same file is used by processor-cleanet; keep them in step.)

suppressPackageStartupMessages({
  library(flowCore)
  library(jsonlite)
})

log_line <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...), "\n", sep = "")
  flush(stdout())
}

# ---- parameters ---------------------------------------------------------------------

env_string <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (identical(v, "")) default else trimws(v)
}

env_bool <- function(name, default) {
  v <- tolower(env_string(name, if (default) "true" else "false"))
  if (v %in% c("true", "1", "yes", "y", "t")) return(TRUE)
  if (v %in% c("false", "0", "no", "n", "f")) return(FALSE)
  stop(sprintf("%s must be true or false, got '%s'", name, v))
}

env_number <- function(name, default, min = -Inf, max = Inf, integer = FALSE) {
  raw <- env_string(name, as.character(default))
  v <- suppressWarnings(as.numeric(raw))
  if (is.na(v) || v < min || v > max || (integer && v != round(v)))
    stop(sprintf("%s must be a %s in [%s, %s], got '%s'", name,
                 if (integer) "whole number" else "number", min, max, raw))
  if (integer) as.integer(v) else v
}

# ---- inputs -------------------------------------------------------------------------

# Files under dir (following symlinks), skipping hidden files such as macOS "._x.fcs".
list_files <- function(dir, pattern) {
  if (!dir.exists(dir)) return(character())
  f <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE,
                  ignore.case = TRUE)
  f[!startsWith(basename(f), ".")]
}

# A port's files: INPUT_DIR/<port>/ when it exists (linked inputs, placement subfolder),
# else anywhere under INPUT_DIR (no links: everything side by side).
port_files <- function(input_dir, port, pattern) {
  sub <- file.path(input_dir, port)
  if (dir.exists(sub)) {
    f <- list_files(sub, pattern)
    if (length(f)) return(sort(f))
  }
  sort(list_files(input_dir, pattern))
}

find_fcs <- function(input_dir) port_files(input_dir, "fcs", "\\.fcs$")

# panel.csv from the metadata port (or anywhere); NULL when there is none.
read_panel <- function(input_dir) {
  f <- port_files(input_dir, "metadata", "^panel\\.csv$")
  if (!length(f)) return(NULL)
  if (length(f) > 1) log_line("note: several panel.csv found, using %s", f[1])
  p <- utils::read.csv(f[1], check.names = FALSE, stringsAsFactors = FALSE,
                       colClasses = "character")
  need <- c("channel", "type")
  if (!all(need %in% names(p)))
    stop(sprintf("%s lacks column(s) %s", f[1], paste(setdiff(need, names(p)), collapse = ", ")))
  p$type <- tolower(trimws(p$type))
  attr(p, "path") <- f[1]
  p
}

# The channels of one type in a flowFrame. With a panel: its rows of that type that the
# file has. Without one (fallback): markers are channels with a $PnS that are not
# scatter/time; scatter is FSC*/SSC*.
panel_channels <- function(ff, panel, type) {
  pd <- pData(parameters(ff))
  names <- as.character(pd$name)
  if (!is.null(panel)) return(intersect(panel$channel[panel$type == type], names))
  desc <- as.character(pd$desc)
  desc[is.na(desc)] <- ""
  scatter <- grepl("^(FSC|SSC)", names, ignore.case = TRUE)
  time <- grepl("^time$", names, ignore.case = TRUE)
  switch(type,
         marker = names[nzchar(desc) & desc != names & !scatter & !time],
         scatter = names[scatter],
         time = names[time],
         character())
}

sample_id_of <- function(path) sub("\\.fcs$", "", basename(path), ignore.case = TRUE)

# ---- FCS I/O ------------------------------------------------------------------------

read_fcs <- function(path) {
  # No transformation and no truncation at $PnR: values stay exactly as written upstream.
  suppressWarnings(flowCore::read.FCS(path, transformation = FALSE,
                                      truncate_max_range = FALSE, alter.names = FALSE))
}

# Keywords flowCore adds on read (not part of the file) and on write (recomputed).
.reader_keywords <- c("FILENAME", "GUID", "ORIGINALGUID", "transformation", "FCSversion",
                      "$BEGINDATA", "$ENDDATA", "$BEGINSTEXT", "$ENDSTEXT",
                      "$BEGINANALYSIS", "$ENDANALYSIS", "$NEXTDATA")

# Write ff with its own keywords plus `extra` (named list) as FCS 3.1 float32.
# flowCore writes an "FCS3.0" header; its TEXT segment (float, $PnE 0,0, $BYTEORD,
# '|' delimiter with escaping) is also valid FCS 3.1, so only the version is restamped.
write_fcs <- function(ff, path, extra = list()) {
  kw <- flowCore::keyword(ff)
  kw <- kw[!(names(kw) %in% .reader_keywords) & !grepl("^flowCore_", names(kw))]
  for (k in names(extra)) kw[[k]] <- as.character(extra[[k]])
  kw[["$FIL"]] <- basename(path)
  ff@description <- kw
  tmp <- paste0(path, ".part")
  suppressWarnings(flowCore::write.FCS(ff, tmp))
  con <- file(tmp, "r+b")
  stopifnot(identical(readChar(con, 6, useBytes = TRUE), "FCS3.0"))
  seek(con, 0, rw = "write")
  writeChar("FCS3.1", con, eos = NULL, useBytes = TRUE)
  close(con)
  file.rename(tmp, path)
  invisible(path)
}

# ---- provenance ---------------------------------------------------------------------

app_version <- function() {
  for (f in c("/app/VERSION", "VERSION")) if (file.exists(f)) return(trimws(readLines(f, 1)))
  "unknown"
}

pkg_versions <- function(pkgs) {
  v <- lapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                         error = function(e) NA_character_))
  names(v) <- pkgs
  c(list(R = paste(R.version$major, R.version$minor, sep = ".")), v)
}

write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null",
                       na = "null")
}
