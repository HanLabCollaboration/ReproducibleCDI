test_that("computeOriginalCDI stops for years after 2022 before any API call", {
  expect_error(computeOriginalCDI(2023), "2017-2022")
  expect_error(computeOriginalCDI(2024), "B09019")
  expect_error(computeOriginalCDI(2023), "Use computeCDI\\(\\) for 2023")
})

test_that("computeOriginalCDI stops for years before 2017", {
  expect_error(computeOriginalCDI(2016), "B28002")
})

test_that("computeCDI requires year >= 2017", {
  expect_error(computeCDI("tract", 2016, "RI"), "2017")
})

test_that("Puerto Rico is part of the CMS national run only", {
  expect_true("PR" %in% ReproducibleCDI:::resolveCDIStates("all", include_pr = TRUE))
  expect_false("PR" %in% ReproducibleCDI:::resolveCDIStates("all"))
  expect_error(ReproducibleCDI:::resolveCDIStates("PR"), "Invalid")
})
