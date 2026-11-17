#' Fetch raw ACS data for Community Deprivation Index (CDI) variables
#'
#' @description
#' Internal helper that queries the Census API for every ACS cell needed to
#' compute the CDI variables defined in `cdiVarSpec()` (the CMS reference
#' definitions), plus total population (`B01001_001`). Every cell comes from a
#' detailed (`B`/`C`) ACS table; the Census Bureau does not tabulate the Subject
#' (`S`) or Data Profile (`DP`) families below the tract, so they are never used.
#'
#' Values are returned exactly as tidycensus provides them.
#'
#' @param geo Geography level: `"state"`, `"county"`, `"tract"`, or
#'   `"block group"` (`"cbg"`).
#' @param year The ACS5 end year.
#' @param states A vector of state abbreviations.
#' @param get_geometry Logical, whether to download geometry.
#' @param moe Logical. If `TRUE`, also keep the margin-of-error (`M`) columns.
#' @param suppression Logical. If `TRUE`, also fetch the cells used for the CMS
#'   block-group suppression flags (see `cdiSuppressionCells()`).
#'
#' @return A wide tidycensus data frame with `GEOID`, `NAME`, the `<cell>E`
#'   columns (and `<cell>M` columns when `moe = TRUE`), plus geometry if
#'   requested.
#'
#' @importFrom dplyr select ends_with any_of
fetchCDIAcs <- function(geo, year, states, get_geometry = FALSE, moe = FALSE,
                        suppression = FALSE) {

  geo <- if (geo == "cbg") "block group" else geo

  vars <- c("B01001_001", cdiAcsCells())
  if (suppression) vars <- c(vars, unname(cdiSuppressionCells(year)))
  vars <- unique(vars)

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

  keep_cols <- c(paste0(vars, "E"), if (moe) paste0(vars, "M"))
  CDI_data_raw |>
    dplyr::select(GEOID, NAME, dplyr::any_of(keep_cols))
}

#' Stop if a CDI table is not published at the requested geography
#'
#' @description
#' The CDI uses exactly the CMS ACS cells; no other table is ever substituted.
#' If every value of a needed cell is missing across the whole request, the
#' Census Bureau does not publish that table at this geography / year, and the
#' function stops (at block group, the CDI for that year can only be computed
#' at the tract level).
#'
#' @param data Raw wide data from `fetchCDIAcs()`.
#' @param geo Geography level requested.
#' @param year ACS5 end year.
#' @param cells ACS cells to check (default: all CDI cells).
#' @return `invisible(TRUE)`, or an error.
#' @noRd
checkCDIAvailability <- function(data, geo, year, cells = cdiAcsCells()) {
  if (nrow(data) == 0) return(invisible(TRUE))
  plain <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else data
  cols <- intersect(paste0(cells, "E"), colnames(plain))
  missing_cols <- sub("E$", "", cols[vapply(cols, function(cl) all(is.na(plain[[cl]])), logical(1))])
  if (length(missing_cols) > 0) {
    tables <- unique(sub("_\\d+$", "", missing_cols))
    hint <- if (geo %in% c("block group", "cbg")) " Use geo = 'tract' for this year." else ""
    stop(glue::glue(
      "Table(s) {glue::glue_collapse(tables, sep = ', ')} are not published at the ",
      "'{geo}' level for ACS5 {year} (every value is missing). The CDI uses the CMS ",
      "variables exactly and does not substitute other tables.{hint}"
    ), call. = FALSE)
  }
  invisible(TRUE)
}
