test_that("write_fcs keeps channels and keywords, stamps FCS3.1 and adds keywords", {
  ff <- make_frame(500)
  f <- tempfile(fileext = ".fcs")
  write_fcs(ff, f, list(PENNSIEVE_DOUBLETS = "cleanet;removed=0"))
  expect_identical(readChar(f, 6, useBytes = TRUE), "FCS3.1")
  back <- read_fcs(f)
  expect_identical(colnames(back), colnames(ff))
  expect_equal(as.vector(flowCore::exprs(back)), as.vector(flowCore::exprs(ff)), tolerance = 1e-6)
  kw <- flowCore::keyword(back)
  expect_identical(kw[["PENNSIEVE_DOUBLETS"]], "cleanet;removed=0")
  expect_identical(kw[["PENNSIEVE_QC"]], "peacoqc;removed=0")
  expect_identical(kw[["$P5S"]], "CD4")
})

test_that("inputs are found in subfolders and flat, following symlinks", {
  linked <- make_input("linked")
  expect_length(find_fcs(linked), 3)
  flat <- make_input("flat")
  expect_length(find_fcs(flat), 3)
  sl <- tempfile("sl")
  dir.create(sl)
  file.symlink(file.path(linked, "fcs"), file.path(sl, "fcs"))
  file.symlink(file.path(linked, "metadata"), file.path(sl, "metadata"))
  expect_length(find_fcs(sl), 3)
  expect_setequal(panel_channels(make_frame(10), read_panel(sl), "marker"), names(markers))
})

test_that("markers go back to linear values before Cleanet adds them", {
  raw <- matrix(c(-500, 0, 1e3, 5e4, 2e5, 10), ncol = 2, dimnames = list(NULL, c("a", "b")))
  x <- asinh(raw / COF)
  cf <- c(a = COF, b = COF)
  expect_equal(linear_matrix(x, cf, 1000), raw)
  expect_equal(linear_matrix(x, cf, "panel"), raw / COF)
  # A channel without a cofactor is already linear.
  expect_equal(linear_matrix(raw, c(a = NA, b = NA), 1000), raw)
})

test_that("bad parameters are rejected", {
  Sys.setenv(COFACTOR = "-3"); on.exit(Sys.unsetenv("COFACTOR"))
  expect_error(cleanet_params(), "COFACTOR")
  Sys.setenv(COFACTOR = "500")
  expect_identical(cleanet_params()$COFACTOR, 500)
})

test_that("cofactors can come from the PENNSIEVE_COFACTORS keyword", {
  ff <- make_frame(10)
  kw <- flowCore::keyword(ff)
  kw[["PENNSIEVE_COFACTORS"]] <- "BV421-A=6000,PE-A=150"
  ff@description <- kw
  expect_identical(cofactors_from_keyword(ff), c(`BV421-A` = 6000, `PE-A` = 150))
  expect_length(cofactors_from_keyword(make_frame(10)), 0)
})
