#' Compute the Community Deprivation Index (CDI)
#'
#' @description
#' Fetches and processes ACS data to calculate the Community Deprivation Index
#' (CDI), a composite measure of neighborhood-level socioeconomic deprivation,
#' following the conventions of the sibling ReproducibleX packages
#' (ReproducibleYostIndex, ReproducibleNDI, ReproducibleSDI).
#'
#' The variables are exactly the 18 used by the official CMS reference
#' implementation (CMSgov/community-deprivation-index). The function automates:
#' 1) data fetching from `tidycensus`, 2) calculation of the variables (a rate
#' with a zero denominator is missing, as in the sibling packages), 3) optional
#' CMS stabilization and/or spatial imputation, in the order you choose, 4) a
#' one-component principal component analysis whose scoring coefficients give
#' the factor-weighted composite, 5) rescaling to mean 100 / SD 20, and 6)
#' ranking into percentiles 1-100.
#'
#' To reproduce every CMS step exactly (hierarchical imputation, suppression,
#' SAS ranking), use [computeOriginalCDI()] instead.
#'
#' @param geo The geography level. One of `"county"`, `"tract"`,
#'   `"block group"` (or `"cbg"`). The CDI is conventionally a block-group index.
#' @param year The 5-year ACS survey end year (e.g., 2022). Must be >= 2017
#'   because the household-internet table B28002 is not published earlier.
#' @param states A vector of state abbreviations (e.g., `c("CA", "NY")`) or
#'   `"all"` for the entire US.
#' @param stabilize Logical. If `TRUE`, applies the CMS MOE-based empirical
#'   Bayes stabilization ([stabilizeCDI()]) toward the parent level: block group
#'   toward tract (as CMS), tract toward county, or county toward state (as in
#'   ReproducibleYostIndex). Only observed values with a margin of error are
#'   stabilized.
#' @param impute Logical. If `TRUE`, performs spatial (queen-contiguity) mean
#'   imputation of missing data ([imputeMissing()]), as in the sibling packages.
#' @param order Order of the two steps when both `stabilize` and `impute` are
#'   `TRUE`: `"stabilize_first"` (default) stabilizes the observed values and
#'   then imputes the remaining gaps from the stabilized neighbors;
#'   `"impute_first"` imputes from the raw neighbors and then stabilizes
#'   (imputed cells are not stabilized themselves, but they enter the parent
#'   heterogeneity used to stabilize the observed cells).
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
#'   \item `df_cdi`: primary data frame with the component variables (after
#'     stabilization / imputation), the raw composite `CDI`, the rescaled
#'     `CDI_std` (mean 100 / SD 20), the percentile `our_CDI` (1-100),
#'     `nvar_stabilized`, `nvar_imputed`, and `annotation`. With
#'     `stabilize = TRUE` it also has `<var>_raw`, `<var>_raw_moe`, and
#'     `stabilizedflg_<var>` (the `inc_dis_upd` MOE is on the ratio scale).
#'   \item `data_harmonized`: the complete-case data fed into the PCA.
#'   \item `data_without_imputation`: the raw component data.
#'   \item `data_stabilized`: the [stabilizeCDI()] output with diagnostics
#'     (`NULL` unless `stabilize = TRUE`).
#'   \item `obj_pca`: the `psych::principal` object (or a named list of objects
#'     when `scope = "state"`).
#'   \item `settings`: the method choices used.
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
    stabilize = FALSE,
    impute = FALSE,
    order = c("stabilize_first", "impute_first"),
    quiet = FALSE,
    return_format = "detailed",
    scope = "national",
    weight_var = "tot_pop",
    ...) {

  # --- 1. Argument checks ----------------------------------------------------
  geo <- rlang::arg_match(geo, values = c("county", "tract", "block group", "cbg"))
  if (geo == "cbg") geo <- "block group"
  order <- rlang::arg_match(order, values = c("stabilize_first", "impute_first"))
  return_format <- rlang::arg_match(return_format, values = c("minimal", "detailed"))
  scope <- rlang::arg_match(scope, values = c("national", "state"))
  weight_var <- rlang::arg_match(weight_var, values = c("tot_pop", "none"))
  stopifnot(is.logical(stabilize), length(stabilize) == 1, !is.na(stabilize))
  stopifnot(is.logical(impute), length(impute) == 1, !is.na(impute))
  if (!is.numeric(year) || length(year) != 1 || is.na(year) || year < 2017) {
    stop("`year` must be an ACS5 end year >= 2017: table B28002 (household internet), ",
         "one of the 18 CDI variables, is not published before ACS5 2017.", call. = FALSE)
  }

  states <- resolveCDIStates(states)

  # --- 2. Data fetching and processing ---------------------------------------
  if (!quiet) message("Fetching raw ACS data...")

  CDI_data_raw <- fetchCDIAcs(geo, year, states, get_geometry = impute, moe = stabilize)
  checkCDIAvailability(CDI_data_raw, geo, year)

  if (!quiet) message("Calculating CDI component variables...")
  CDI_data <- cleanAcsNamesCDI(dframe = CDI_data_raw, geo = geo) |>
    calculateCDIVars(zero_denominator = "na", with_moe = stabilize)

  CDI_parent <- NULL
  geo_parent_key_len <- NULL
  if (stabilize) {
    geo_hierarchy <- list(
      "block group" = list(parent = "tract",  key_len = 11),
      "tract"       = list(parent = "county", key_len = 5),
      "county"      = list(parent = "state",  key_len = 2)
    )
    geo_parent <- geo_hierarchy[[geo]]$parent
    geo_parent_key_len <- geo_hierarchy[[geo]]$key_len
    if (!quiet) message(glue::glue("Fetching {geo_parent}-level data for stabilization..."))
    CDI_parent <- fetchCDIAcs(geo_parent, year, states) |>
      calculateCDIVars(zero_denominator = "na")
  }

  run_scope <- function(data_sub) {
    processCDIScope(
      CDI_data_sub = data_sub,
      impute_sub = impute,
      stabilize_sub = stabilize,
      order_sub = order,
      CDI_parent = CDI_parent,
      geo_parent_key_len = geo_parent_key_len,
      quiet_sub = quiet,
      return_format_sub = return_format,
      weight_var_sub = weight_var
    )
  }

  # --- 3. Main branching logic -----------------------------------------------
  if (scope == "state" && length(states) > 1) {
    if (!quiet) message("Scope is 'state'. Processing each state separately.")

    list_of_state_data <- dplyr::group_split(CDI_data, state)
    results_list <- purrr::map(list_of_state_data, run_scope)

    results_list <- purrr::keep(results_list, ~ !is.null(.x))
    CDI_output_data <- purrr::map_dfr(results_list, "final_data")
    data_harmonized <- purrr::map_dfr(results_list, "data_harmonized")
    data_without_imputation <- purrr::map_dfr(results_list, "data_with_no_imputation")
    data_stabilized <- if (stabilize) purrr::map_dfr(results_list, "data_stabilized") else NULL
    list_of_pca_objects <- purrr::map(results_list, "pca_object")

    names(list_of_pca_objects) <- purrr::map_chr(results_list, ~ unique(.x$final_data$state)[1])
    out_obj_pca <- list_of_pca_objects

  } else {
    if (!quiet) message("Scope is 'national' or only one state selected. Processing all geographies together.")

    single_result <- run_scope(CDI_data)

    CDI_output_data <- single_result$final_data
    data_harmonized <- single_result$data_harmonized
    data_without_imputation <- single_result$data_with_no_imputation
    data_stabilized <- single_result$data_stabilized
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
  out$data_stabilized <- data_stabilized
  out$obj_pca <- out_obj_pca
  out$settings <- list(
    geo = geo, year = year, stabilize = stabilize, impute = impute,
    order = if (stabilize && impute) order else NA_character_,
    scope = scope, weight_var = weight_var
  )

  return(out)
}
