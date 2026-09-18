#' Compute the Community Deprivation Index (CDI)
#'
#' @description
#' Fetches and processes ACS data to calculate the Community Deprivation Index
#' (CDI), a composite measure of neighborhood-level socioeconomic deprivation.
#'
#' The CDI follows the Singh / Area Deprivation Index (ADI) methodology and adds
#' a percentage-uninsured variable, for a total of 18 standardized ACS
#' variables. The function automates: 1) data fetching from `tidycensus`,
#' 2) calculation of the 18 component variables, 3) optional spatial imputation,
#' 4) a principal component analysis whose first-component loadings provide the
#' factor weights, 5) construction of the factor-weighted composite, 6)
#' rescaling to mean 100 / SD 20, and 7) ranking into national percentiles 1-100.
#'
#' The CDI is designed to be computed natively at the **census block group**
#' level; every component variable is drawn from a detailed (`B`/`C`) ACS table
#' so that block-group tabulation is available.
#'
#' @param geo The geography level. One of `"county"`, `"tract"`,
#'   `"block group"` (or `"cbg"`). The CDI is conventionally a block-group index.
#' @param year The 5-year ACS survey end year (e.g., 2022). Must be >= 2016
#'   because the household-internet variable (B28002) is not published earlier.
#' @param states A vector of state abbreviations (e.g., `c("CA", "NY")`) or
#'   `"all"` for the entire US.
#' @param impute Logical. If `TRUE`, performs spatial (queen-contiguity) mean
#'   imputation for missing data, as in the sibling Yost/NDI packages.
#' @param quiet Logical. If `TRUE`, suppresses messages and warnings.
#' @param return_format `"detailed"` (default) or `"minimal"`. `"minimal"`
#'   returns only `GEOID`, `CDI_std`, `our_CDI`, and `annotation`.
#' @param scope `"national"` (default) or `"state"`. `"national"` ranks all
#'   geographies against each other; `"state"` fits the PCA and ranks within
#'   each state separately.
#' @param weight_var Variable for weighted spatial imputation: `"tot_pop"`
#'   (default) or `"none"`. Only relevant when `impute = TRUE`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' If `return_format = "minimal"`, a data frame of `GEOID`, `CDI_std`,
#' `our_CDI`, and `annotation`. Otherwise a list with:
#' \itemize{
#'   \item `df_cdi`: primary data frame with the component variables, the raw
#'     composite `CDI`, the rescaled `CDI_std` (mean 100 / SD 20), the national
#'     percentile `our_CDI` (1-100), and `annotation`.
#'   \item `data_harmonized`: the standardized component data fed into the PCA.
#'   \item `data_without_imputation`: the raw, non-imputed component data.
#'   \item `obj_pca`: the `psych::principal` object (or a named list of objects
#'     when `scope = "state"`).
#' }
#'
#' @references
#' Robst, J., Forlines, G., Kautter, J., et al. (2024). The development of the
#' Community Deprivation Index and its application to Accountable Care
#' Organizations. \emph{Health Affairs Scholar}. \doi{10.1093/haschl/qxae161}
#'
#' Singh, G.K. (2003). Area deprivation and widening inequalities in US
#' mortality, 1969-1998. \emph{American Journal of Public Health}, 93(7),
#' 1137-1143. \doi{10.2105/ajph.93.7.1137}
#'
#' Kind, A.J.H., & Buckingham, W.R. (2018). Making neighborhood-disadvantage
#' metrics accessible - the Neighborhood Atlas. \emph{New England Journal of
#' Medicine}, 378(26), 2456-2458. \doi{10.1056/NEJMp1802313}
#'
#' @importFrom glue glue glue_collapse
#' @importFrom dplyr group_split select any_of arrange
#' @importFrom purrr map map_dfr map_chr keep
#' @export
computeCDI <- function(
    geo,
    year,
    states,
    impute = FALSE,
    quiet = FALSE,
    return_format = "detailed",
    scope = "national",
    weight_var = "tot_pop",
    ...) {

  # --- 1. Argument checks ----------------------------------------------------
  geo <- rlang::arg_match(geo, values = c("county", "tract", "block group", "cbg"))
  return_format <- rlang::arg_match(return_format, values = c("minimal", "detailed"))
  scope <- rlang::arg_match(scope, values = c("national", "state"))
  weight_var <- rlang::arg_match(weight_var, values = c("tot_pop", "none"))
  stopifnot(is.numeric(year), year >= 2016)

  get_geometry <- impute

  # Validate / expand the requested states.
  all_geos <- unique(tidycensus::fips_codes$state)
  non_states <- c("PR", "AS", "GU", "MP", "VI", "UM")
  valid_states <- setdiff(all_geos, non_states)

  if (states[1] == "all") {
    states <- valid_states
  } else {
    invalid_states <- setdiff(states, valid_states)
    if (length(invalid_states) > 0) {
      invalid_states_str <- glue::glue_collapse(invalid_states, sep = ", ")
      stop(glue::glue(
        "Invalid state abbreviations provided: {invalid_states_str}.
         Please use standard 2-letter postal codes for the 50 US states."
      ))
    }
  }

  # --- 2. Data fetching and processing ---------------------------------------
  if (!quiet) message("Fetching raw ACS data...")

  CDI_data_raw <- fetchCDIAcs(geo, year, states, get_geometry)
  acs_vars <- attr(CDI_data_raw, "acs_vars")

  CDI_data_clean <- cleanAcsNamesCDI(dframe = CDI_data_raw, geo = geo)

  if (!quiet) message("Calculating CDI component variables...")
  CDI_data <- CDI_data_clean |>
    calculateCDIVars(acs_vars = acs_vars)

  # --- 3. Main branching logic -----------------------------------------------
  if (scope == "state" && length(states) > 1) {
    if (!quiet) message("Scope is 'state'. Processing each state separately.")

    list_of_state_data <- dplyr::group_split(CDI_data, state)

    results_list <- purrr::map(list_of_state_data, ~processCDIScope(
      CDI_data_sub = .x,
      impute_sub = impute,
      quiet_sub = quiet,
      return_format_sub = return_format,
      weight_var_sub = weight_var
    ))

    results_list <- purrr::keep(results_list, ~ !is.null(.x))
    CDI_output_data <- purrr::map_dfr(results_list, "final_data")
    data_harmonized <- purrr::map_dfr(results_list, "data_harmonized")
    data_without_imputation <- purrr::map_dfr(results_list, "data_with_no_imputation")
    list_of_pca_objects <- purrr::map(results_list, "pca_object")

    names(list_of_pca_objects) <- purrr::map_chr(results_list, ~ unique(.x$final_data$state)[1])
    out_obj_pca <- list_of_pca_objects

  } else {
    if (!quiet) message("Scope is 'national' or only one state selected. Processing all geographies together.")

    single_result <- processCDIScope(
      CDI_data_sub = CDI_data,
      impute_sub = impute,
      quiet_sub = quiet,
      return_format_sub = return_format,
      weight_var_sub = weight_var
    )

    CDI_output_data <- single_result$final_data
    data_harmonized <- single_result$data_harmonized
    data_without_imputation <- single_result$data_with_no_imputation
    out_obj_pca <- single_result$pca_object
  }

  # --- 4. Final output preparation -------------------------------------------
  if (return_format == "minimal") {
    return(
      CDI_output_data |>
        dplyr::select(GEOID, CDI_std, our_CDI, annotation) |>
        dplyr::arrange(GEOID)
    )
  }

  out <- list()
  out$df_cdi <- CDI_output_data |> dplyr::arrange(GEOID)
  out$data_harmonized <- data_harmonized
  out$data_without_imputation <- data_without_imputation
  out$obj_pca <- out_obj_pca

  return(out)
}
