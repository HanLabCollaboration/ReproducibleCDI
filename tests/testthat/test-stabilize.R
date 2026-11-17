toy_tract <- function() {
  bg <- data.frame(
    GEOID = c("010010201001", "010010201002", "010010201003", "010010201004", "010010202001"),
    x     = c(0.10, 0.30, NA, 0.20, 0.40),
    x_moe = c(0.05, 0.10, 0.10, NA, 0.05)
  )
  tract <- data.frame(GEOID = c("01001020100", "01001020200"), x = c(0.20, 0.50))
  list(bg = bg, tract = tract)
}

test_that("stabilizeCDI reproduces the CMS formula", {
  d <- toy_tract()
  out <- stabilizeCDI(d$bg, d$tract, "x", geo_parent_key_len = 11)

  # Tract 1: 4 rows (missing value still counted); sqdiff over observed rows only
  sqdiff <- (0.10 - 0.20)^2 + (0.30 - 0.20)^2 + (0.20 - 0.20)^2
  t2 <- sqdiff / (4 - 1)
  s2 <- (0.05 / 1.645)^2
  w <- (1 / s2) / (1 / s2 + 1 / t2)

  expect_equal(out$parent_n, c(4, 4, 4, 4, 1))
  expect_equal(out$x_sqdiff[1], sqdiff)
  expect_equal(out$x_t2[1], t2)
  expect_equal(out$x_wgt[1], w)
  expect_equal(out$x_stabilized[1], w * 0.10 + (1 - w) * 0.20)
  expect_equal(out$stabilizedflg_x[1:2], c(1L, 1L))
})

test_that("stabilizeCDI leaves values unchanged in every CMS exception", {
  d <- toy_tract()
  out <- stabilizeCDI(d$bg, d$tract, "x", geo_parent_key_len = 11)
  # row 3: value missing; row 4: MOE missing; row 5: only unit in its parent
  expect_true(is.na(out$x_stabilized[3]))
  expect_equal(out$x_stabilized[4], 0.20)
  expect_equal(out$x_stabilized[5], 0.40)
  expect_equal(out$stabilizedflg_x[3:5], c(0L, 0L, 0L))
  expect_true(all(is.na(out$x_wgt[3:5])))
  expect_equal(out$nvar_stabilized, c(1L, 1L, 0L, 0L, 0L))

  # parent value missing -> unchanged
  d$tract$x[1] <- NA
  out2 <- stabilizeCDI(d$bg, d$tract, "x", geo_parent_key_len = 11)
  expect_equal(out2$x_stabilized[1:2], c(0.10, 0.30))

  # zero MOE (s2 = 0) -> unchanged
  d <- toy_tract()
  d$bg$x_moe[1] <- 0
  out3 <- stabilizeCDI(d$bg, d$tract, "x", geo_parent_key_len = 11)
  expect_equal(out3$x_stabilized[1], 0.10)
  expect_equal(out3$stabilizedflg_x[1], 0L)
})

test_that("stabilizeCDI requires MOE columns", {
  d <- toy_tract()
  expect_error(stabilizeCDI(d$bg[, 1:2], d$tract, "x", 11), "MOE")
})
