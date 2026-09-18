#' Calculate the 18 CDI component variables
#'
#' @description
#' Internal helper that derives the 18 Community Deprivation Index component
#' variables from the raw ACS estimates returned by `fetchCDIAcs()`.
#'
#' Direction of each variable (positive = more deprived, negative = less
#' deprived) is *not* enforced here - it is learned from the data by the
#' principal component analysis in `processCDIScope()`. The percentage/level
#' variables kept here match the naming used in the CDI supplement:
#' \itemize{
#'   \item `educ_12less_perc`  - pop 25+ with < 12 years schooling (no diploma)
#'   \item `educ_ba_pl_perc`   - pop 25+ with a bachelor's degree or higher
#'   \item `emp_white_col_perc`- civilian employed 16+ in mgmt/bus/sci/arts
#'   \item `fam_bel_pov_perc`  - families below the poverty line
#'   \item `hhld_1pl_room_perc`- occupied units with > 1 person per room
#'   \item `hhld_no_int_perc`  - households with no internet access
#'   \item `hhld_no_veh_perc`  - households with no vehicle available
#'   \item `hous_no_plumb_perc`- units lacking complete plumbing
#'   \item `inc_dis_upd_imp`   - income disparity (log ratio, updated)
#'   \item `med_gross_rent_imp`- median gross rent
#'   \item `med_hhld_inc_imp`  - median household income
#'   \item `med_home_val_imp`  - median home value
#'   \item `med_month_mort_imp`- median monthly owner cost with a mortgage
#'   \item `no_ins_perc`       - population with no health insurance
#'   \item `one_par_hhld_perc_imp` - single-parent family households
#'   \item `own_occ_hous_perc` - owner-occupied housing units
#'   \item `pop_bel_150_pov_perc`  - population below 150% of poverty
#'   \item `unemp_perc`        - unemployed share of the civilian labor force
#' }
#'
#' @param data The cleaned data frame from `cleanAcsNamesCDI()`.
#' @param acs_vars The `acs_vars` list attached by `fetchCDIAcs()`.
#'
#' @return A data frame with GEOID, geography columns, `tot_pop`, and the 18
#'   component variables (plus geometry if present).
#'
#' @importFrom dplyr mutate across all_of any_of where select case_when
calculateCDIVars <- function(data, acs_vars) {

  # Helper: sum a set of estimate columns row-wise (na.rm mimics SAS `sum()`,
  # which skips missing cells).
  rs <- function(df, cols) rowSums(dplyr::across(dplyr::all_of(cols)), na.rm = TRUE)

  # Helper: form a rate following the CMS reference convention - a zero
  # denominator yields 0 (not missing), a missing denominator yields missing.
  pct <- function(num, denom) dplyr::case_when(
    is.na(denom) ~ NA_real_,
    denom == 0 ~ 0,
    TRUE ~ num / denom
  )

  # Helper: median/level columns use Census "jam" codes for non-computable
  # cells (large negative sentinels). Treat any negative value as missing.
  clean_median <- function(x) ifelse(!is.na(x) & x < 0, NA_real_, x)

  data_calculated <- data |>
    dplyr::mutate(
      tot_pop = .data[[acs_vars$total_pop]],

      # --- percentage variables (CMS reference definitions) ---
      educ_12less_perc = pct(rs(data, acs_vars$educ_12less_num), .data[[acs_vars$educ_denom]]),
      educ_ba_pl_perc  = pct(rs(data, acs_vars$educ_ba_pl_num),  .data[[acs_vars$educ_denom]]),
      emp_white_col_perc = pct(rs(data, acs_vars$emp_white_col_num), .data[[acs_vars$emp_denom]]),
      fam_bel_pov_perc = pct(.data[[acs_vars$fam_pov_num]], .data[[acs_vars$fam_pov_denom]]),
      hhld_1pl_room_perc = pct(rs(data, acs_vars$crowd_num), .data[[acs_vars$crowd_denom]]),
      hhld_no_int_perc = pct(.data[[acs_vars$int_no_access]], .data[[acs_vars$int_denom]]),
      hhld_no_veh_perc = pct(rs(data, acs_vars$noveh_num), .data[[acs_vars$noveh_denom]]),
      hous_no_plumb_perc = pct(rs(data, acs_vars$noplumb_num), .data[[acs_vars$noplumb_denom]]),
      no_ins_perc = pct(rs(data, acs_vars$no_ins_num), .data[[acs_vars$no_ins_denom]]),
      one_par_hhld_perc_imp = pct(rs(data, acs_vars$onepar_num), .data[[acs_vars$onepar_denom]]),
      own_occ_hous_perc = pct(.data[[acs_vars$ownocc_num]], .data[[acs_vars$ownocc_denom]]),
      pop_bel_150_pov_perc = pct(rs(data, acs_vars$pov150_num), .data[[acs_vars$pov150_denom]]),
      unemp_perc = pct(.data[[acs_vars$unemp_num]], .data[[acs_vars$unemp_denom]]),

      # --- income disparity (updated), per the CMS reference code ---
      # natural log of 100 * (households < $20k) / (households > $100k).
      # Missing when either the low- or high-income count is zero.
      .incdis_low  = rs(data, acs_vars$incdis_low),
      .incdis_high = rs(data, acs_vars$incdis_high),
      inc_dis_upd_imp = ifelse(
        .data[[".incdis_high"]] > 0 & .data[[".incdis_low"]] > 0,
        log(100 * (.data[[".incdis_low"]] / .data[[".incdis_high"]])),
        NA_real_
      ),

      # --- median / level variables (jam codes -> NA) ---
      med_gross_rent_imp = clean_median(.data[[acs_vars$med_gross_rent]]),
      med_hhld_inc_imp   = clean_median(.data[[acs_vars$med_hhld_inc]]),
      med_home_val_imp   = clean_median(.data[[acs_vars$med_home_val]]),
      med_month_mort_imp = clean_median(.data[[acs_vars$med_month_mort]])
    ) |>
    # NaN (0/0) -> NA so imputation and PCA treat them as missing.
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ ifelse(is.nan(.), NA, .))) |>
    dplyr::select(
      GEOID, NAME, tot_pop,
      dplyr::any_of(c("block_group", "tract", "county", "state")),
      educ_12less_perc, educ_ba_pl_perc, emp_white_col_perc, fam_bel_pov_perc,
      hhld_1pl_room_perc, hhld_no_int_perc, hhld_no_veh_perc, hous_no_plumb_perc,
      inc_dis_upd_imp, med_gross_rent_imp, med_hhld_inc_imp, med_home_val_imp,
      med_month_mort_imp, no_ins_perc, one_par_hhld_perc_imp, own_occ_hous_perc,
      pop_bel_150_pov_perc, unemp_perc,
      dplyr::any_of("geometry")
    )

  return(data_calculated)
}
