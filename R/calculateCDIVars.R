#' Calculate the CDI component variables
#'
#' @description
#' Internal helper that derives the CDI variables from raw ACS estimates
#' returned by `fetchCDIAcs()`, following the CMS reference definitions in
#' `cdiVarSpec()`:
#' \itemize{
#'   \item `educ_12less_perc`  - pop 25+ with < 12 years schooling (no diploma)
#'   \item `educ_hs_pl_perc`   - pop 25+ with a high-school diploma or more
#'     (computed by CMS but **not** used in the PCA)
#'   \item `educ_ba_pl_perc`   - pop 25+ with a bachelor's degree or higher
#'   \item `emp_white_col_perc`- civilian employed 16+ in mgmt/bus/sci/arts or
#'     sales/office occupations
#'   \item `fam_bel_pov_perc`  - families below the poverty line
#'   \item `hhld_1pl_room_perc`- occupied units with > 1 person per room
#'   \item `hhld_no_int_perc`  - households with no internet access
#'   \item `hhld_no_veh_perc`  - households with no vehicle available
#'   \item `hous_no_plumb_perc`- units lacking complete plumbing
#'   \item `inc_dis_upd`       - income disparity, 100 x households < $20k /
#'     households > $100k, **on the ratio scale** (missing when either count is
#'     0). It is log-transformed by `finalizeIncDis()` after stabilization and
#'     imputation, as in the CMS code.
#'   \item `med_gross_rent`    - median gross rent
#'   \item `med_hhld_inc`      - median household income
#'   \item `med_home_val`      - median home value
#'   \item `med_month_mort`    - median monthly owner cost with a mortgage
#'   \item `no_ins_perc`       - population with no health insurance
#'   \item `one_par_hhld_perc` - single-parent family households
#'   \item `own_occ_hous_perc` - owner-occupied housing units
#'   \item `pop_bel_150_pov_perc`  - population below 150% of poverty
#'   \item `unemp_perc`        - unemployed share of the civilian labor force
#' }
#'
#' Rates are proportions (0-1). Medians are used exactly as tidycensus returns
#' them (Census jam codes are already `NA`; top/bottom-coded medians are numeric
#' caps).
#'
#' @param data The raw (optionally name-cleaned) wide data from `fetchCDIAcs()`.
#' @param zero_denominator How to treat a rate whose denominator is 0:
#'   `"zero"` sets it to 0 (CMS rule, used by `computeOriginalCDI()`); `"na"`
#'   sets it to missing (rule of the sibling ReproducibleX packages, used by
#'   `computeCDI()`).
#' @param with_moe Logical. If `TRUE`, also return `<var>_moe` columns computed
#'   with the CMS MOE formulas. Requires the `<cell>M` columns.
#'
#' @return A data frame with `GEOID`, any name / geography columns present,
#'   `tot_pop`, the 19 variables, optional `_moe` columns, and geometry if
#'   present.
#'
#' @importFrom dplyr select any_of all_of bind_cols
calculateCDIVars <- function(data, zero_denominator = c("zero", "na"), with_moe = FALSE) {
  zero_denominator <- match.arg(zero_denominator)
  spec <- cdiVarSpec()
  plain <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else data
  plain <- as.data.frame(plain)
  E <- function(x) paste0(x, "E")

  out <- list()
  for (i in seq_len(nrow(spec))) {
    v   <- spec$var[i]
    num <- spec$num[[i]]
    den <- spec$den[[i]]

    if (spec$type[i] == "median") {
      out[[v]] <- as.numeric(plain[[E(num)]])
      if (with_moe) out[[paste0(v, "_moe")]] <- as.numeric(plain[[paste0(num, "M")]])
      next
    }

    n_est <- sasRowSum(plain, E(num))
    d_est <- sasRowSum(plain, E(den))

    if (spec$type[i] == "rate") {
      zero_value <- if (zero_denominator == "zero") 0 else NA_real_
      out[[v]] <- ifelse(is.na(d_est), NA_real_,
                         ifelse(d_est == 0, zero_value, n_est / d_est))
    } else {
      # Inc_dis_upd: defined only when both the low- and high-income counts are
      # non-zero (the CMS code leaves it missing otherwise).
      ok <- !is.na(n_est) & !is.na(d_est) & n_est != 0 & d_est != 0
      out[[v]] <- ifelse(ok, 100 * n_est / d_est, NA_real_)
    }

    if (with_moe) {
      moe <- cdiMoeRatio(n_est, d_est, cdiCombinedMoe(plain, num), cdiCombinedMoe(plain, den))
      # CMS multiplies the income-disparity MOE by 100 to match the ratio scale.
      if (spec$type[i] == "ratio") moe <- moe * 100
      out[[paste0(v, "_moe")]] <- moe
    }
  }

  vars <- dplyr::bind_cols(out)
  # Order: the 19 values first, then their MOEs.
  vars <- vars[, c(spec$var, if (with_moe) paste0(spec$var, "_moe")), drop = FALSE]

  keep <- data |>
    dplyr::select(
      GEOID,
      dplyr::any_of(c("NAME", "block_group", "tract", "county", "state"))
    )
  keep$tot_pop <- as.numeric(plain[["B01001_001E"]])

  data_calculated <- dplyr::bind_cols(keep, vars)
  if (inherits(data, "sf")) {
    data_calculated <- sf::st_as_sf(data_calculated)
  }
  data_calculated
}

#' Log-transform the income-disparity variable
#'
#' @description
#' Final step for `inc_dis_upd`, applied after stabilization and imputation as
#' in the CMS code: natural log of the ratio. Non-positive ratios become `NA`.
#'
#' @param data A data frame with an `inc_dis_upd` column on the ratio scale.
#' @return `data` with `inc_dis_upd` on the log scale.
#' @noRd
finalizeIncDis <- function(data) {
  x <- data$inc_dis_upd
  data$inc_dis_upd <- ifelse(!is.na(x) & x > 0, log(x), NA_real_)
  data
}
