#' Clean ACS geographic names
#'
#' @description
#' Internal helper that splits the `NAME` column returned by tidycensus into
#' clean columns for block group, tract, county, and state.
#'
#' @param dframe The raw data from `fetchCDIAcs()`.
#' @param geo The geography level ("county", "tract", "block group"/"cbg").
#'
#' @return A data frame with cleaned name columns.
#'
#' @importFrom tidyr separate
#' @importFrom dplyr mutate across all_of any_of
#' @importFrom stringr str_squish str_extract str_remove
cleanAcsNamesCDI <- function(dframe, geo) {

  geo_norm <- ifelse(geo == "cbg", "block group", geo)
  sep_cols <- switch(
    geo_norm,
    "block group" = c("block_group", "tract", "county", "state"),
    "tract" = c("tract", "county", "state"),
    "county" = c("county", "state"),
    stop("Unknown geo type: ", geo_norm)
  )

  data_cleaned <- dframe |>
    tidyr::separate(
      col    = NAME,
      into   = sep_cols,
      sep    = "; |, ",
      remove = FALSE,
      fill   = "right",
      extra  = "merge"
    ) |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(sep_cols),
        .fns  = stringr::str_squish
      ),
      dplyr::across(
        .cols = dplyr::any_of(c("block_group", "tract")),
        .fns  = ~ stringr::str_extract(.x, "\\d+\\.?\\d*")
      ),
      county = stringr::str_remove(county, "\\s*County$")
    )

  return(data_cleaned)
}
