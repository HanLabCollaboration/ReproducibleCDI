#' Process a single scope for the Community Deprivation Index (CDI)
#'
#' @description
#' Primary internal workhorse for `computeCDI()`. It processes a single "scope"
#' (all national data, or the data for a single state):
#' \enumerate{
#'   \item filter geographies with no population;
#'   \item optionally stabilize (CMS shrinkage toward the parent level, via
#'     [stabilizeCDI()]) and/or impute (spatial queen-contiguity, via
#'     [imputeMissing()]), in the order given by `order_sub`;
#'   \item log-transform the income-disparity ratio;
#'   \item on complete cases, standardize the 18 variables, fit a
#'     one-component PCA, and score with the scoring coefficients,
#'     sign-aligned so higher = more deprived ([fitCDIPCA()]);
#'   \item rescale to mean 100, SD 20 (`CDI_std`); and
#'   \item rank into percentiles 1-100 (`our_CDI`, `dplyr::ntile()`).
#' }
#'
#' Order of stabilization and imputation (only matters when both are on):
#' \itemize{
#'   \item `"stabilize_first"`: observed values are stabilized (missing values
#'     stay missing and add nothing to the parent heterogeneity), then the
#'     remaining gaps are imputed from the stabilized neighbors.
#'   \item `"impute_first"`: gaps are imputed from the raw neighbors, then the
#'     data are stabilized. Imputed cells have no MOE, so they are not
#'     stabilized themselves, but they enter the parent heterogeneity \eqn{t^2}
#'     used to stabilize their observed siblings.
#' }
#'
#' @param CDI_data_sub A data frame subset for a single scope (optionally `sf`).
#' @param impute_sub Logical, whether to perform spatial imputation.
#' @param stabilize_sub Logical, whether to stabilize toward the parent level.
#' @param order_sub `"stabilize_first"` or `"impute_first"`.
#' @param CDI_parent Parent-level variables (required when `stabilize_sub`).
#' @param geo_parent_key_len `GEOID` prefix length of the parent level.
#' @param quiet_sub Logical, whether to suppress messages and warnings.
#' @param return_format_sub The desired return format (handled by the parent).
#' @param weight_var_sub A string passed to `imputeMissing()` ("tot_pop" or
#'   "none"). Only used when `impute_sub = TRUE`.
#'
#' @return A list with:
#' * `final_data`: full data frame with `CDI`, `CDI_std`, `our_CDI`, annotation.
#' * `pca_object`: the `psych::principal` object (or `NULL` if PCA failed).
#' * `pca_coefficients`: the sign-aligned scoring coefficients.
#' * `data_with_no_imputation`: the raw, non-imputed input for the scope.
#' * `data_harmonized`: the complete-case data fed into the PCA, with
#'   standardized `std_` columns.
#' * `data_stabilized`: the [stabilizeCDI()] output with diagnostics (or
#'   `NULL`).
#'
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr filter select any_of mutate across all_of left_join
#' @importFrom dplyr everything bind_rows ntile case_when
#' @importFrom glue glue
processCDIScope <- function(
    CDI_data_sub,
    impute_sub = FALSE,
    stabilize_sub = FALSE,
    order_sub = "stabilize_first",
    CDI_parent = NULL,
    geo_parent_key_len = NULL,
    quiet_sub = FALSE,
    return_format_sub = "detailed",
    weight_var_sub = "tot_pop") {

  # The 18 CDI PCA variables, in CMS order.
  varlist <- cdiVarNames(pca_only = TRUE)
  moe_cols <- paste0(varlist, "_moe")

  # Hold the input before manipulating (income disparity on the log scale).
  raw_CDI_data_sub <- CDI_data_sub |> sf::st_drop_geometry() |> finalizeIncDis()

  # Filter out no-population areas for this subset.
  CDI_data_sub <- CDI_data_sub |> dplyr::filter(tot_pop > 0)

  rows_with_no_pop <- nrow(raw_CDI_data_sub) - nrow(CDI_data_sub)
  if (!quiet_sub) message(glue::glue("... there are {rows_with_no_pop} geoids with no population."))
  if (nrow(CDI_data_sub) == 0) return(NULL)

  if (impute_sub && !inherits(CDI_data_sub, "sf")) {
    stop("Imputation requires geometry. Please re-run with impute = TRUE in computeCDI().")
  }
  geometry_tbl <- if (impute_sub) CDI_data_sub |> dplyr::select(GEOID) else NULL

  work <- CDI_data_sub |>
    sf::st_drop_geometry() |>
    dplyr::select(GEOID, tot_pop, dplyr::all_of(varlist), dplyr::any_of(moe_cols))
  work$nvar_stabilized <- 0L
  work$nvar_imputed <- 0L
  data_stabilized <- NULL

  # --- Stabilization (CMS shrinkage toward the parent level) -----------------
  stabilize_step <- function(w) {
    st <- stabilizeCDI(w, CDI_parent, varlist, geo_parent_key_len)
    for (v in varlist) w[[v]] <- st[[paste0(v, "_stabilized")]]
    w$nvar_stabilized <- st$nvar_stabilized
    data_stabilized <<- st
    w
  }

  # --- Spatial (queen-contiguity) imputation ---------------------------------
  impute_step <- function(w) {
    if (!quiet_sub) message("... imputing missing data for scope.")
    to_impute <- geometry_tbl |>
      dplyr::left_join(w |> dplyr::select(GEOID, tot_pop, dplyr::all_of(varlist)), by = "GEOID")
    imp <- imputeMissing(dfGeo = to_impute, weight_var_sub = weight_var_sub, quiet = quiet_sub)
    idx <- match(w$GEOID, imp$GEOID)
    for (v in varlist) w[[v]] <- imp[[v]][idx]
    w$nvar_imputed <- as.integer(imp$nvar_imputed[idx])
    w
  }

  if (order_sub == "impute_first") {
    if (impute_sub) {
      was_missing <- is.na(as.matrix(work[, varlist]))
      work <- impute_step(work)
      # Imputed cells carry no MOE, so they are never stabilized themselves.
      if (stabilize_sub) {
        for (j in seq_along(varlist)) work[[moe_cols[j]]][was_missing[, j]] <- NA_real_
      }
    }
    if (stabilize_sub) work <- stabilize_step(work)
  } else {
    if (stabilize_sub) work <- stabilize_step(work)
    if (impute_sub) work <- impute_step(work)
  }

  work <- finalizeIncDis(work)
  work$nvar_still_missing <- as.integer(rowSums(is.na(as.matrix(work[, varlist]))))
  work <- work |> dplyr::select(-dplyr::any_of(moe_cols))

  # --- Standardize + PCA on complete cases -----------------------------------
  current_scope <- "Current Scope"
  if ("state" %in% colnames(raw_CDI_data_sub)) {
    current_scope <- unique(raw_CDI_data_sub$state)[1]
  }

  CDI_inp_data_fa <- work |>
    dplyr::filter(nvar_still_missing == 0) |>
    dplyr::select(GEOID, tot_pop, dplyr::all_of(varlist))

  pca <- fitCDIPCA(CDI_inp_data_fa, varlist, quiet = quiet_sub, scope_label = current_scope)

  CDI_inp_data_fa <- CDI_inp_data_fa |>
    dplyr::mutate(dplyr::across(dplyr::all_of(varlist), ~ scale(.)[, 1], .names = "std_{.col}"))
  CDI_inp_data_fa$CDI <- pca$score
  CDI_inp_data_fa$CDI_std <- pca$score_std

  # --- Assemble output -------------------------------------------------------
  raw_cols <- raw_CDI_data_sub |> dplyr::select(-dplyr::any_of(moe_cols))
  cols_from_raw <- setdiff(colnames(raw_cols), colnames(work))

  CDI_output_data <- work |>
    dplyr::left_join(
      CDI_inp_data_fa |> dplyr::select(GEOID, CDI, CDI_std),
      by = "GEOID"
    ) |>
    dplyr::left_join(
      raw_cols |> dplyr::select(GEOID, dplyr::all_of(cols_from_raw)),
      by = "GEOID"
    ) |>
    dplyr::select(GEOID, dplyr::all_of(cols_from_raw), dplyr::everything())

  if (stabilize_sub) {
    # Yost naming: pre-stabilization values and MOEs as <var>_raw / _raw_moe.
    raw_vals <- raw_CDI_data_sub[match(CDI_output_data$GEOID, raw_CDI_data_sub$GEOID), ]
    for (v in varlist) {
      CDI_output_data[[paste0(v, "_raw")]] <- raw_vals[[v]]
      CDI_output_data[[paste0(v, "_raw_moe")]] <- raw_vals[[paste0(v, "_moe")]]
    }
    flags <- data_stabilized |>
      dplyr::select(GEOID, dplyr::all_of(paste0("stabilizedflg_", varlist)))
    CDI_output_data <- CDI_output_data |> dplyr::left_join(flags, by = "GEOID")
  }

  # Bring back the no-population geographies.
  CDI_output_data <- CDI_output_data |>
    dplyr::bind_rows(raw_cols |> dplyr::filter(tot_pop == 0))

  # National percentile ranking 1-100 (unweighted; each geography = one unit).
  CDI_output_data <- CDI_output_data |>
    dplyr::mutate(our_CDI = dplyr::ntile(CDI_std, 100))

  # Annotate GEOIDs by data quality / stabilization / imputation status.
  CDI_output_data <- CDI_output_data |>
    dplyr::mutate(
      annotation = dplyr::case_when(
        tot_pop == 0 ~ "No population",
        !is.na(CDI) & nvar_imputed == 0 & nvar_stabilized == 0 ~ "Complete data",
        !is.na(CDI) & nvar_imputed == 0 & nvar_stabilized > 0 ~ "Stabilized",
        !is.na(CDI) & nvar_imputed > 0 & nvar_stabilized == 0 ~ "Imputation completed",
        !is.na(CDI) & nvar_imputed > 0 & nvar_stabilized > 0 ~ "Stabilized + imputation completed",
        is.na(CDI) & nvar_still_missing > 0 ~ "Imputation incomplete",
        is.na(CDI) & nvar_still_missing == 0 ~ "PCA failed",
        TRUE ~ "Other"
      )
    )

  return(list(
    final_data = CDI_output_data,
    pca_object = pca$fit,
    pca_coefficients = pca$coefficients,
    data_with_no_imputation = raw_CDI_data_sub,
    data_harmonized = CDI_inp_data_fa,
    data_stabilized = data_stabilized
  ))
}
