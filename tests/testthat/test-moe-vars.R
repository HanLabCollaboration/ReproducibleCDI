test_that("proportion MOE follows the CMS formula, both branches", {
  # x_sub > 0: sqrt(10^2 - 0.5^2 * 4^2) / 20
  expect_equal(ReproducibleCDI:::cdiMoeRatio(10, 20, 10, 4), sqrt(100 - 0.25 * 16) / 20)
  # x_sub <= 0: falls back to the additive formula
  expect_equal(ReproducibleCDI:::cdiMoeRatio(10, 20, 1, 10), sqrt(1 + 0.25 * 100) / 20)
  # zero or missing denominator -> NA
  expect_true(is.na(ReproducibleCDI:::cdiMoeRatio(0, 0, 5, 5)))
  expect_true(is.na(ReproducibleCDI:::cdiMoeRatio(1, NA, 5, 5)))
})

test_that("combined MOE is the root sum of squares (no zero-cell rule)", {
  df <- data.frame(aM = 3, bM = 4)
  expect_equal(ReproducibleCDI:::cdiCombinedMoe(df, c("a", "b")), 5)
  expect_equal(ReproducibleCDI:::cdiCombinedMoe(df, "a"), 3)
})

test_that("zero denominators: 0 for CMS, NA for computeCDI", {
  df <- toy_acs(1)
  df$B17010_001E <- 0
  df$B17010_002E <- 0
  zero <- ReproducibleCDI:::calculateCDIVars(df, "zero")
  na   <- ReproducibleCDI:::calculateCDIVars(df, "na")
  expect_equal(zero$fam_bel_pov_perc, 0)
  expect_true(is.na(na$fam_bel_pov_perc))
})

test_that("rates, sums, and income disparity on the ratio scale", {
  df <- toy_acs(2)
  out <- ReproducibleCDI:::calculateCDIVars(df, "zero", with_moe = TRUE)
  # 15 cells of 100 over 100
  expect_equal(out$educ_12less_perc, c(15, 15))
  # 100 * (3 * 100) / (4 * 100)
  expect_equal(out$inc_dis_upd, c(75, 75))
  expect_equal(out$tot_pop, c(100, 100))

  # MOE of the ratio is x100 and uses multi-cell numerator and denominator
  p <- 300 / 400
  expected <- 100 * sqrt(3 * 100 - p^2 * 4 * 100) / 400
  expect_equal(out$inc_dis_upd_moe, c(expected, expected))

  # missing when either the low or the high count is 0, in both modes
  df$B19001_002E <- 0; df$B19001_003E <- 0; df$B19001_004E <- 0
  expect_true(all(is.na(ReproducibleCDI:::calculateCDIVars(df, "zero")$inc_dis_upd)))
  expect_true(all(is.na(ReproducibleCDI:::calculateCDIVars(df, "na")$inc_dis_upd)))
})

test_that("medians and their MOEs are passed through", {
  df <- toy_acs(1)
  df$B19013_001E <- 250001
  df$B19013_001M <- NA
  out <- ReproducibleCDI:::calculateCDIVars(df, "zero", with_moe = TRUE)
  expect_equal(out$med_hhld_inc, 250001)
  expect_true(is.na(out$med_hhld_inc_moe))
})

test_that("finalizeIncDis takes the log and drops non-positive ratios", {
  out <- ReproducibleCDI:::finalizeIncDis(data.frame(inc_dis_upd = c(100, 0, NA)))
  expect_equal(out$inc_dis_upd, c(log(100), NA, NA))
})
