#' Validate and expand requested states
#'
#' @param states A vector of 2-letter state abbreviations, or `"all"`.
#' @param include_pr Logical. If `TRUE`, Puerto Rico (`"PR"`) is valid and is
#'   part of `"all"` (the CMS national ranking includes it). Otherwise `"all"`
#'   is the 50 states and DC, as in the sibling ReproducibleX packages.
#' @return A character vector of state abbreviations.
#' @importFrom glue glue glue_collapse
#' @noRd
resolveCDIStates <- function(states, include_pr = FALSE) {
  all_geos <- unique(tidycensus::fips_codes$state)
  non_states <- c("PR", "AS", "GU", "MP", "VI", "UM")
  if (include_pr) non_states <- setdiff(non_states, "PR")
  valid_states <- setdiff(all_geos, non_states)

  if (states[1] == "all") return(valid_states)

  invalid_states <- setdiff(states, valid_states)
  if (length(invalid_states) > 0) {
    invalid_states_str <- glue::glue_collapse(invalid_states, sep = ", ")
    allowed <- if (include_pr) "the 50 US states, DC, and PR" else "the 50 US states and DC"
    stop(glue::glue(
      "Invalid state abbreviations provided: {invalid_states_str}.
       Please use standard 2-letter postal codes for {allowed}."
    ))
  }
  states
}
