#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data :=
## usethis namespace: end

utils::globalVariables(c(
  # shared
  "GEOID", "NAME", "CDI", "annotation", "block_group", "county", "tract",
  "state", "tot_pop", "geometry", "nvar_imputed", "nvar_still_missing",
  "nvar_stabilized", "suppression_flag",
  # the 18 CDI PCA variables (+ educ_hs_pl_perc, computed but not in the PCA)
  "educ_12less_perc", "educ_hs_pl_perc", "educ_ba_pl_perc", "emp_white_col_perc",
  "fam_bel_pov_perc", "hhld_1pl_room_perc", "hhld_no_int_perc",
  "hhld_no_veh_perc", "hous_no_plumb_perc", "inc_dis_upd",
  "med_gross_rent", "med_hhld_inc", "med_home_val",
  "med_month_mort", "no_ins_perc", "one_par_hhld_perc",
  "own_occ_hous_perc", "pop_bel_150_pov_perc", "unemp_perc",
  # intermediate / output
  "CDI_score", "CDI_std", "our_CDI", "cdi_rank", "CDI_NATRANK"
))
