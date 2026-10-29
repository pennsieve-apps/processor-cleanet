# Unit tests: Rscript tests/testthat.R from the repo root (or /app/tests/testthat.R in the
# image). Needs R with Cleanet, flowCore, jsonlite and testthat; `make test` runs them in the image.
library(testthat)
if (!nzchar(Sys.getenv("APP_ROOT"))) {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(APP_ROOT = normalizePath(file.path(dirname(f[1]), "..")))
}
test_dir(file.path(Sys.getenv("APP_ROOT"), "tests", "testthat"), reporter = "summary",
         stop_on_failure = TRUE)
