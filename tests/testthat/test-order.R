# 6 x 6 grid of block groups in 6 tracts (one tract per grid row).
toy_scope <- function(missing = TRUE) {
  set.seed(42)
  vars <- ReproducibleCDI:::cdiVarNames(TRUE)
  n <- 36
  tract <- sprintf("010010%05d", rep(1:6, each = 6) * 100)
  df <- data.frame(GEOID = paste0(tract, rep(1:6, times = 6)), tot_pop = 1000)
  f <- rnorm(n)
  for (v in vars) {
    df[[v]] <- abs(f * runif(1, 0.5, 1) + rnorm(n, sd = 0.5)) + 1
    df[[paste0(v, "_moe")]] <- runif(n, 0.05, 0.5)
  }
  if (missing) {
    df$med_hhld_inc[c(3, 17, 30)] <- NA
    df$fam_bel_pov_perc[c(8, 22)] <- NA
  }
  df$state <- "AL"
  sf_df <- sf::st_sf(df, geometry = toy_grid(6, 6))

  parent <- data.frame(GEOID = unique(tract))
  for (v in vars) parent[[v]] <- tapply(df[[v]], tract, mean, na.rm = TRUE)[parent$GEOID]
  list(data = sf_df, parent = parent, vars = vars)
}

run_order <- function(toy, order) {
  ReproducibleCDI:::processCDIScope(
    toy$data, impute_sub = TRUE, stabilize_sub = TRUE, order_sub = order,
    CDI_parent = toy$parent, geo_parent_key_len = 11, quiet_sub = TRUE,
    weight_var_sub = "tot_pop"
  )
}

test_that("with no missing data the two orders are identical", {
  toy <- toy_scope(missing = FALSE)
  a <- run_order(toy, "stabilize_first")$final_data
  b <- run_order(toy, "impute_first")$final_data
  expect_equal(a$CDI_std, b$CDI_std)
})

test_that("with missing data the orders differ, and imputed cells are never stabilized", {
  toy <- toy_scope(missing = TRUE)
  sf_first <- run_order(toy, "stabilize_first")
  if_first <- run_order(toy, "impute_first")

  expect_false(isTRUE(all.equal(sf_first$final_data$CDI_std, if_first$final_data$CDI_std)))
  expect_true(all(!is.na(sf_first$final_data$CDI_std)))
  expect_true(all(!is.na(if_first$final_data$CDI_std)))

  # impute_first: originally-missing cells were imputed but not stabilized
  st <- if_first$data_stabilized
  miss <- c(3, 17, 30)
  expect_true(all(!is.na(st$med_hhld_inc[miss])))
  expect_equal(st$stabilizedflg_med_hhld_inc[miss], c(0L, 0L, 0L))

  # stabilize_first: missing cells stay missing during stabilization
  st1 <- sf_first$data_stabilized
  expect_true(all(is.na(st1$med_hhld_inc_stabilized[miss])))

  ann <- sf_first$final_data$annotation
  expect_true(all(ann[miss] == "Stabilized + imputation completed"))
})
