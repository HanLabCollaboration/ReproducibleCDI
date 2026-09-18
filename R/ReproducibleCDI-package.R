#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data :=
## usethis namespace: end

utils::globalVariables(c(
  # shared
  "GEOID", "NAME", "CDI", "annotation", "block_group", "county", "tract",
  "state", "tot_pop", "geometry", "nvar_imputed", "nvar_still_missing",
  # the 18 CDI component variables
  "educ_12less_perc", "educ_ba_pl_perc", "emp_white_col_perc",
  "fam_bel_pov_perc", "hhld_1pl_room_perc", "hhld_no_int_perc",
  "hhld_no_veh_perc", "hous_no_plumb_perc", "inc_dis_upd_imp",
  "med_gross_rent_imp", "med_hhld_inc_imp", "med_home_val_imp",
  "med_month_mort_imp", "no_ins_perc", "one_par_hhld_perc_imp",
  "own_occ_hous_perc", "pop_bel_150_pov_perc", "unemp_perc",
  # intermediate / output
  "CDI_score", "CDI_std", "our_CDI"
))
