run_app <- function(input, env = character()) {
  out <- tempfile("out")
  old <- Sys.getenv(names(env), unset = NA, names = TRUE)
  if (length(env)) do.call(Sys.setenv, as.list(env))
  on.exit(for (k in names(env)) {
    if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k))
  })
  status <- run(input, out)
  list(status = status, out = out,
       qc = utils::read.csv(file.path(out, "qc", "doublets.csv"), stringsAsFactors = FALSE))
}

for (layout in c("linked", "flat")) test_that(paste("end to end,", layout, "inputs"), {
  r <- run_app(make_input(layout), c(THREADS = "1"))
  expect_identical(r$status, 0L)
  expect_setequal(list.files(file.path(r$out, "fcs")), c("s01.fcs", "s02.fcs", "tiny.fcs", "provenance.json"))
  qc <- r$qc
  expect_identical(names(qc)[1:4], c("sample_id", "events_in", "doublets", "pct_doublets"))
  s01 <- qc[qc$sample_id == "s01", ]
  expect_identical(s01$status, "ok")
  expect_identical(s01$events_in, 8400L)
  # 400 heterotypic + homotypic doublets were planted (4.8%); homotypic ones of the same
  # population look like bright singlets, so not all are found.
  expect_gt(s01$doublets, 150)
  expect_lt(s01$doublets, 800)
  expect_identical(s01$events_out, s01$events_in - s01$doublets)
  ff <- read_fcs(file.path(r$out, "fcs", "s01.fcs"))
  expect_identical(nrow(ff), as.integer(s01$events_out))
  expect_identical(colnames(ff), colnames(make_frame(10)))
  expect_match(flowCore::keyword(ff)[["PENNSIEVE_DOUBLETS"]], "^cleanet;doublets=")
  expect_identical(flowCore::keyword(ff)[["PENNSIEVE_QC"]], "peacoqc;removed=0")
  expect_false(is.unsorted(flowCore::exprs(ff)[, "Time"]))
  # Few calls where there are none.
  expect_lt(qc$pct_doublets[qc$sample_id == "s02"], 2)
  tiny <- qc[qc$sample_id == "tiny", ]
  expect_identical(tiny$status, "passed_through")
  expect_identical(tiny$events_out, 300L)
  expect_true(file.exists(file.path(r$out, "qc", "plots", "s01.png")))
  prov <- jsonlite::read_json(file.path(r$out, "qc", "provenance.json"))
  expect_identical(prov$app, "processor-cleanet")
  expect_length(prov$samples, 3)
})

test_that("the planted doublets are the ones called", {
  ff <- make_frame(8000, n_doublets = 400, seed = 1)
  truth <- attr(ff, "is_doublet")
  x <- linear_matrix(flowCore::exprs(ff)[, names(markers)], setNames(rep(COF, 6), names(markers)), "panel")
  set.seed(42)
  res <- suppressMessages(Cleanet::cleanet(as.data.frame(x), names(markers), cofactor = 1))
  called <- res$status == "Doublet"
  expect_gt(sum(called & truth) / sum(called), 0.8)  # precision
})

test_that("same result with files in parallel, and SEED makes it reproducible", {
  input <- make_input("linked")
  a <- run_app(input, c(THREADS = "1", PLOTS = "false"))$qc
  b <- run_app(input, c(THREADS = "3", PLOTS = "false"))$qc
  expect_identical(a[, 1:5], b[, 1:5])
})

test_that("KEEP_DOUBLETS=true keeps every event", {
  r <- run_app(make_input("linked"), c(THREADS = "1", PLOTS = "false", KEEP_DOUBLETS = "true"))
  s01 <- r$qc[r$qc$sample_id == "s01", ]
  expect_gt(s01$doublets, 0)
  expect_identical(s01$events_out, s01$events_in)
})

test_that("no input fails the run", {
  empty <- tempfile("empty")
  dir.create(empty)
  expect_identical(run(empty, tempfile("out")), 1L)
})
