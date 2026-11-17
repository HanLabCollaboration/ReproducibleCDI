#' Flag block groups for CMS suppression
#'
#' @description
#' Port of the suppression step in the CMS reference code
#' (`cdi_construction_s1.sas`):
#' \itemize{
#'   \item `PH` (low population / housing): total population `B01001_001` < 100
#'     or housing units `B25032_001` < 30.
#'   \item `GQ` (group quarters): population in group quarters / total
#'     population (`B09019`) >= 0.3333.
#'   \item `GQ-PH`: both.
#' }
#' Missing counts follow SAS comparison rules (a missing value is smaller than
#' any number, so a missing population or housing count flags `PH`).
#'
#' @param data Raw wide block-group data from `fetchCDIAcs(suppression = TRUE)`.
#' @param year ACS5 end year (selects the group-quarters cell, see
#'   `cdiSuppressionCells()`).
#'
#' @return A data frame with `GEOID` and `suppression_flag` (`NA`, `"PH"`,
#'   `"GQ"`, or `"GQ-PH"`).
flagCDISuppression <- function(data, year) {
  cells <- paste0(cdiSuppressionCells(year), "E")
  names(cells) <- names(cdiSuppressionCells(year))
  plain <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else data

  pop      <- plain[[cells[["pop"]]]]
  hu       <- plain[[cells[["hu"]]]]
  gq_total <- plain[[cells[["gq_total"]]]]
  gq       <- plain[[cells[["gq"]]]]

  gq_share <- ifelse(!is.na(gq_total) & gq_total > 0 & !is.na(gq), gq / gq_total, NA_real_)

  ph_flag <- is.na(pop) | pop < 100 | is.na(hu) | hu < 30
  gq_flag <- (is.na(gq) | gq != 0 | is.na(gq_total) | gq_total != 0) &
    !is.na(gq_share) & gq_share >= 0.3333

  data.frame(
    GEOID = plain$GEOID,
    suppression_flag = dplyr::case_when(
      ph_flag & gq_flag ~ "GQ-PH",
      gq_flag ~ "GQ",
      ph_flag ~ "PH",
      TRUE ~ NA_character_
    ),
    stringsAsFactors = FALSE
  )
}
