test_that("hierarchical imputation falls back tract -> county -> state", {
  bg <- data.frame(
    GEOID = c("010010201001", "010010201002", "010010202001", "010030301001"),
    x = c(1, NA, NA, NA),
    inc_dis_upd = c(0, 5, NA, NA)
  )
  parents <- list(
    tract  = data.frame(GEOID = c("01001020100", "01001020200"), x = c(2, NA), inc_dis_upd = c(7, 0)),
    county = data.frame(GEOID = c("01001", "01003"), x = c(3, NA), inc_dis_upd = c(8, NA)),
    state  = data.frame(GEOID = "01", x = 4, inc_dis_upd = 9)
  )
  out <- ReproducibleCDI:::imputeHierarchical(bg, parents, c("x", "inc_dis_upd"))
  expect_equal(out$x, c(1, 2, 3, 4))
  expect_equal(out$x_imp_src, c("self", "tract", "county", "state"))
  # a 0 ratio is not usable: falls back to the next level
  expect_equal(out$inc_dis_upd, c(7, 5, 8, 9))
  expect_equal(out$inc_dis_upd_imp_src, c("tract", "self", "county", "state"))
  expect_equal(out$nvar_imputed, c(1L, 1L, 2L, 2L))
  expect_equal(out$nvar_still_missing, c(0L, 0L, 0L, 0L))
})

test_that("CMS suppression flags", {
  df <- data.frame(
    GEOID = as.character(1:5),
    B01001_001E = c(500, 50, 500, 50, NA),
    B25032_001E = c(200, 200, 200, 200, 200),
    B09019_001E = c(500, 50, 500, 50, 500),
    B09019_026E = c(0, 0, 200, 40, 0)
  )
  out <- ReproducibleCDI:::flagCDISuppression(df, 2019)
  expect_equal(out$suppression_flag, c(NA, "PH", "GQ", "GQ-PH", "PH"))

  # threshold is 0.3333, not 1/3
  df2 <- data.frame(GEOID = "1", B01001_001E = 10000, B25032_001E = 500,
                    B09019_001E = 10000, B09019_038E = 3333)
  expect_equal(ReproducibleCDI:::flagCDISuppression(df2, 2018)$suppression_flag, "GQ")
})

test_that("rankGroupsSAS uses floor(rank * k / (n + 1)) + 1", {
  x <- c(5, 1, 3, NA, 3)
  r <- rank(x, ties.method = "average", na.last = "keep")
  expect_equal(ReproducibleCDI:::rankGroupsSAS(x, 100), as.integer(floor(r * 100 / 5) + 1))
  expect_true(is.na(ReproducibleCDI:::rankGroupsSAS(x, 100)[4]))
  expect_equal(range(ReproducibleCDI:::rankGroupsSAS(1:1000, 100)), c(1L, 100L))
})

test_that("PCA coefficients are loadings / eigenvalue, sign-aligned", {
  set.seed(1)
  n <- 200
  f <- rnorm(n)
  vars <- ReproducibleCDI:::cdiVarNames(TRUE)
  df <- as.data.frame(sapply(vars, function(v) f * runif(1, 0.5, 1) + rnorm(n, sd = 0.5)))
  df$fam_bel_pov_perc <- -df$fam_bel_pov_perc # force a negative raw loading
  pca <- ReproducibleCDI:::fitCDIPCA(df, vars, quiet = TRUE)
  expect_gt(pca$loadings[["fam_bel_pov_perc"]], 0)
  expect_equal(unname(abs(pca$coefficients)), unname(abs(pca$loadings / pca$fit$values[1])))
  expect_equal(unname(abs(pca$coefficients)), unname(abs(pca$fit$weights[, 1])), tolerance = 1e-6)
  expect_equal(mean(pca$score_std), 100)
  expect_equal(stats::sd(pca$score_std), 20)
})
