#' CMS Community Deprivation Index variable specification
#'
#' @description
#' Internal single source of truth for the CDI component variables, transcribed
#' cell-for-cell from the official CMS reference implementation
#' (`cdi_construction_s1.sas`, CMSgov/community-deprivation-index). Fetching,
#' variable calculation, margins of error, stabilization, and the PCA variable
#' list are all driven from this table.
#'
#' Columns:
#' \itemize{
#'   \item `var`      - variable name (lowercase CMS name).
#'   \item `type`     - `"rate"` (numerator / denominator), `"median"` (a single
#'     published cell), or `"ratio"` (`Inc_dis_upd`: 100 x low / high income
#'     households, kept on the ratio scale until after stabilization and
#'     imputation, then log-transformed).
#'   \item `num`, `den` - ACS cells (without the `E`/`M` suffix). `den` is empty
#'     for medians.
#'   \item `moe_type` - the CMS MOE macro used (`1` published MOE, `2` single
#'     numerator and denominator cell, `3` sums of cells).
#'   \item `in_pca`   - whether the variable enters the PCA. CMS computes
#'     `Educ_hs_pl_perc` but does not use it.
#' }
#'
#' @return A tibble with one row per CDI variable.
#' @noRd
cdiVarSpec <- function() {
  b <- function(table, idx) paste0(table, "_", sprintf("%03d", idx))

  rows <- list(
    list("educ_12less_perc",     "rate",   b("B15003", 2:16),               "B15003_001", 3, TRUE),
    list("educ_hs_pl_perc",      "rate",   b("B15003", 17:25),              "B15003_001", 3, FALSE),
    list("educ_ba_pl_perc",      "rate",   b("B15003", 22:25),              "B15003_001", 3, TRUE),
    list("emp_white_col_perc",   "rate",   b("C24010", c(3, 39, 27, 63)),   "C24010_001", 3, TRUE),
    list("fam_bel_pov_perc",     "rate",   "B17010_002",                    "B17010_001", 2, TRUE),
    list("hhld_1pl_room_perc",   "rate",   b("B25014", c(5:7, 11:13)),      "B25014_001", 3, TRUE),
    list("hhld_no_int_perc",     "rate",   "B28002_013",                    "B28002_001", 2, TRUE),
    list("hhld_no_veh_perc",     "rate",   b("B25044", c(3, 10)),           "B25044_001", 3, TRUE),
    list("hous_no_plumb_perc",   "rate",   b("B25016", c(7, 16)),           "B25016_001", 3, TRUE),
    list("inc_dis_upd",          "ratio",  b("B19001", 2:4),                b("B19001", 14:17), 3, TRUE),
    list("med_gross_rent",       "median", "B25064_001",                    character(0), 1, TRUE),
    list("med_hhld_inc",         "median", "B19013_001",                    character(0), 1, TRUE),
    list("med_home_val",         "median", "B25077_001",                    character(0), 1, TRUE),
    list("med_month_mort",       "median", "B25088_002",                    character(0), 1, TRUE),
    list("no_ins_perc",          "rate",   b("B27010", c(17, 33, 50, 66)),  "B27010_001", 3, TRUE),
    list("one_par_hhld_perc",    "rate",   b("B11003", c(10, 16)),          "B11003_001", 3, TRUE),
    list("own_occ_hous_perc",    "rate",   "B25003_002",                    "B25003_001", 2, TRUE),
    list("pop_bel_150_pov_perc", "rate",   b("C17002", 2:5),                "C17002_001", 3, TRUE),
    list("unemp_perc",           "rate",   "B23025_005",                    "B23025_002", 2, TRUE)
  )

  dplyr::tibble(
    var      = vapply(rows, `[[`, character(1), 1),
    type     = vapply(rows, `[[`, character(1), 2),
    num      = lapply(rows, `[[`, 3),
    den      = lapply(rows, `[[`, 4),
    moe_type = vapply(rows, function(r) as.integer(r[[5]]), integer(1)),
    in_pca   = vapply(rows, `[[`, logical(1), 6)
  )
}

#' Names of the CDI variables
#'
#' @param pca_only Logical. If `TRUE` (default) return only the 18 variables
#'   used in the PCA; otherwise all 19 computed variables.
#' @return A character vector of variable names, in CMS order.
#' @noRd
cdiVarNames <- function(pca_only = TRUE) {
  spec <- cdiVarSpec()
  if (pca_only) spec$var[spec$in_pca] else spec$var
}

#' ACS cells needed for the CDI variables
#' @return A character vector of ACS cell names (no `E`/`M` suffix).
#' @noRd
cdiAcsCells <- function() {
  spec <- cdiVarSpec()
  unique(c(unlist(spec$num), unlist(spec$den)))
}

#' ACS cells used for the CMS block-group suppression flags
#'
#' @description
#' `PH` uses total population (`B01001_001`) and housing units (`B25032_001`);
#' `GQ` uses the group-quarters share from `B09019`. The group-quarters cell
#' moved when the table was restructured: `B09019_038` up to ACS5 2018,
#' `B09019_026` from ACS5 2019 on.
#'
#' @param year ACS5 end year.
#' @return A named character vector (`pop`, `hu`, `gq_total`, `gq`).
#' @noRd
cdiSuppressionCells <- function(year) {
  c(
    pop      = "B01001_001",
    hu       = "B25032_001",
    gq_total = "B09019_001",
    gq       = if (year <= 2018) "B09019_038" else "B09019_026"
  )
}
