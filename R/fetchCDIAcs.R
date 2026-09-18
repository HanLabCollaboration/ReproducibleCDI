#' Fetch raw ACS data for Community Deprivation Index (CDI) variables
#'
#' @description
#' Internal helper that queries the Census API for all variables required to
#' compute the 18 CDI component variables. The CDI follows the Singh / Area
#' Deprivation Index (ADI) methodology and is designed to be computed natively
#' at the **census block group** level. For that reason every variable is drawn
#' from a **detailed (`B`/`C`) table**: the Census Bureau does not tabulate the
#' Subject (`S`) or Data Profile (`DP`) table families below the census-tract
#' level, so those families cannot be used for a block-group index. (Health
#' insurance uses B27010, which is tabulated at the block-group level, rather
#' than B27001, which stops at the tract.)
#'
#' @param geo Geography level ("county", "tract", "block group"/"cbg").
#' @param year The ACS5 end year.
#' @param states A vector of state abbreviations or 'all'.
#' @param get_geometry Logical, whether to download geometry.
#'
#' @return A tidycensus data frame (wide) with an `acs_vars` attribute that
#'   records the variable groupings used by `calculateCDIVars()`.
#'
#' @importFrom dplyr select ends_with
#' @importFrom stringr str_pad
fetchCDIAcs <- function(geo, year, states, get_geometry) {

  # --- Total population ------------------------------------------------------
  total_pop <- "B01003_001"

  # (1) educ_12less_perc: pop 25+ with < 12 years of schooling / no diploma
  #     B15003 = Educational Attainment for pop 25+.
  educ_denom <- "B15003_001"
  educ_12less_num <- paste0("B15003_", stringr::str_pad(2:16, 3, pad = "0"))

  # (2) educ_ba_pl_perc: pop 25+ with a bachelor's degree or higher
  educ_ba_pl_num <- paste0("B15003_", stringr::str_pad(22:25, 3, pad = "0"))

  # (3) emp_white_col_perc: civilian employed 16+ in "white collar" occupations,
  #     defined (per the CMS reference code) as management/business/science/arts
  #     PLUS sales and office occupations, for both sexes. C24010 = Sex by
  #     Occupation for the civilian employed population 16+.
  #       _003 male mgmt/bus/sci/arts, _027 male sales & office,
  #       _039 female mgmt/bus/sci/arts, _063 female sales & office.
  emp_denom <- "C24010_001"
  emp_white_col_num <- c("C24010_003", "C24010_039", "C24010_027", "C24010_063")

  # (4) fam_bel_pov_perc: families with income below poverty. B17010.
  fam_pov_num <- "B17010_002"
  fam_pov_denom <- "B17010_001"

  # (5) hhld_1pl_room_perc: occupied units with > 1 occupant per room
  #     (crowding). B25014 = Tenure by Occupants per Room.
  crowd_num <- paste0("B25014_", stringr::str_pad(c(5, 6, 7, 11, 12, 13), 3, pad = "0"))
  crowd_denom <- "B25014_001"

  # (6) hhld_no_int_perc: households with no internet access at all (per the CMS
  #     reference code). B28002 = Presence and Types of Internet Subscriptions in
  #     Household. _013 = "No Internet access".
  #     NOTE: B28002 is only published from ACS5 2016 onward.
  int_no_access <- "B28002_013"
  int_denom <- "B28002_001"

  # (7) hhld_no_veh_perc: occupied units with no vehicle available. B25044.
  noveh_num <- c("B25044_003", "B25044_010") # owner + renter, no vehicle
  noveh_denom <- "B25044_001"

  # (8) hous_no_plumb_perc: occupied units lacking complete plumbing (per the
  #     CMS reference code). B25016 = Tenure by Plumbing Facilities by Occupants
  #     per Room. _007 = owner-occupied lacking complete plumbing; _016 =
  #     renter-occupied lacking complete plumbing; _001 = total occupied.
  noplumb_num <- c("B25016_007", "B25016_016")
  noplumb_denom <- "B25016_001"

  # (9) inc_dis_upd_imp: income disparity (updated), per the CMS reference code:
  #     natural log of 100 * (households with income < $20,000) /
  #     (households with income > $100,000). B19001 = Household Income.
  #     _002 = < $10k, _003 = $10-15k, _004 = $15-20k;
  #     _014..017 = $100-125k, $125-150k, $150-200k, $200k+.
  incdis_low <- paste0("B19001_", stringr::str_pad(2:4, 3, pad = "0"))
  incdis_high <- paste0("B19001_", stringr::str_pad(14:17, 3, pad = "0"))

  # (10) med_gross_rent_imp: median gross rent. B25064_001.
  med_gross_rent <- "B25064_001"

  # (11) med_hhld_inc_imp: median household income. B19013_001.
  med_hhld_inc <- "B19013_001"

  # (12) med_home_val_imp: median value of owner-occupied units. B25077_001.
  med_home_val <- "B25077_001"

  # (13) med_month_mort_imp: median selected monthly owner costs (with a
  #      mortgage). B25088_002.
  med_month_mort <- "B25088_002"

  # (14) no_ins_perc: no health insurance coverage. B27010 = Types of Health
  #      Insurance Coverage by Age. This table (unlike B27001, which stops at
  #      the tract level) IS tabulated at the block-group level, so the CDI can
  #      use a native block-group uninsured rate. Numerator = the "No health
  #      insurance coverage" cell for each of the four age brackets.
  no_ins_num <- c("B27010_017", "B27010_033", "B27010_050", "B27010_066")
  no_ins_denom <- "B27010_001"

  # (15) one_par_hhld_perc_imp: single-parent family households with own
  #      children < 18. B11003. _010 = male householder, no spouse present;
  #      _016 = female householder, no spouse present.
  onepar_num <- c("B11003_010", "B11003_016")
  onepar_denom <- "B11003_001"

  # (16) own_occ_hous_perc: owner-occupied housing units. B25003.
  ownocc_num <- "B25003_002"
  ownocc_denom <- "B25003_001"

  # (17) pop_bel_150_pov_perc: population below 150% of the poverty line.
  #      C17002 = Ratio of Income to Poverty Level. _002..005 = under 1.50.
  pov150_num <- paste0("C17002_", stringr::str_pad(2:5, 3, pad = "0"))
  pov150_denom <- "C17002_001"

  # (18) unemp_perc: unemployed share of the labor force 16+ (per the CMS
  #      reference code). B23025. _005 = unemployed; _002 = total in labor force.
  unemp_num <- "B23025_005"
  unemp_denom <- "B23025_002"

  # --- Assemble the full variable request ------------------------------------
  vars <- unique(c(
    total_pop,
    educ_denom, educ_12less_num, educ_ba_pl_num,
    emp_denom, emp_white_col_num,
    fam_pov_num, fam_pov_denom,
    crowd_num, crowd_denom,
    int_no_access, int_denom,
    noveh_num, noveh_denom,
    noplumb_num, noplumb_denom,
    incdis_low, incdis_high,
    med_gross_rent, med_hhld_inc, med_home_val, med_month_mort,
    no_ins_num, no_ins_denom,
    onepar_num, onepar_denom,
    ownocc_num, ownocc_denom,
    pov150_num, pov150_denom,
    unemp_num, unemp_denom
  ))

  # Keep the groupings (with the "E" estimate suffix) for calculateCDIVars().
  E <- function(x) paste0(x, "E")
  acs_vars <- list(
    vars = E(vars),
    total_pop = E(total_pop),
    educ_denom = E(educ_denom),
    educ_12less_num = E(educ_12less_num),
    educ_ba_pl_num = E(educ_ba_pl_num),
    emp_denom = E(emp_denom),
    emp_white_col_num = E(emp_white_col_num),
    fam_pov_num = E(fam_pov_num),
    fam_pov_denom = E(fam_pov_denom),
    crowd_num = E(crowd_num),
    crowd_denom = E(crowd_denom),
    int_no_access = E(int_no_access),
    int_denom = E(int_denom),
    noveh_num = E(noveh_num),
    noveh_denom = E(noveh_denom),
    noplumb_num = E(noplumb_num),
    noplumb_denom = E(noplumb_denom),
    incdis_low = E(incdis_low),
    incdis_high = E(incdis_high),
    med_gross_rent = E(med_gross_rent),
    med_hhld_inc = E(med_hhld_inc),
    med_home_val = E(med_home_val),
    med_month_mort = E(med_month_mort),
    no_ins_num = E(no_ins_num),
    no_ins_denom = E(no_ins_denom),
    onepar_num = E(onepar_num),
    onepar_denom = E(onepar_denom),
    ownocc_num = E(ownocc_num),
    ownocc_denom = E(ownocc_denom),
    pov150_num = E(pov150_num),
    pov150_denom = E(pov150_denom),
    unemp_num = E(unemp_num),
    unemp_denom = E(unemp_denom)
  )

  # --- Fetch -----------------------------------------------------------------
  # Every table used here (including B27010 for health insurance) is a detailed
  # B/C table that the Census Bureau tabulates down to the block-group level, so
  # a single request serves all supported geographies.
  CDI_data_raw <- tidycensus::get_acs(
    geography = geo,
    state = states,
    year = year,
    output = "wide",
    variables = vars,
    survey = "acs5",
    geometry = get_geometry,
    cache_table = TRUE
  )

  CDI_data_raw <- CDI_data_raw |>
    dplyr::select(GEOID, NAME, dplyr::ends_with("E"))
  attr(CDI_data_raw, "acs_vars") <- acs_vars
  return(CDI_data_raw)
}
