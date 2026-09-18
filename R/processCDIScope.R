#' Process a single scope for the Community Deprivation Index (CDI)
#'
#' @description
#' Primary internal workhorse for `computeCDI()`. It processes a single "scope"
#' (all national data, or the data for a single state) and follows the Singh /
#' ADI methodology described in the CDI supplement:
#' \enumerate{
#'   \item filter geographies with no population;
#'   \item optionally impute missing data with spatial (queen-contiguity)
#'     neighbors, exactly as the sibling Yost/NDI packages do;
#'   \item standardize each of the 18 component variables to mean 0, SD 1;
#'   \item run a principal component analysis and take the first component's
#'     loadings as the factor weights \eqn{W_j};
#'   \item form the factor-weighted composite \eqn{S_i = \sum_j X'_{ij} W_j},
#'     sign-aligned so that higher scores mean more deprivation;
#'   \item rescale the composite to mean 100, SD 20 (`CDI_std`); and
#'   \item rank the rescaled score into national percentiles 1-100 (`our_CDI`).
#' }
#'
#' @param CDI_data_sub A data frame subset for a single scope (optionally `sf`).
#' @param impute_sub Logical, whether to perform spatial imputation.
#' @param quiet_sub Logical, whether to suppress messages and warnings.
#' @param return_format_sub The desired return format (handled by the parent).
#' @param weight_var_sub A string passed to `imputeMissing()` ("tot_pop" or
#'   "none"). Only used when `impute_sub = TRUE`.
#'
#' @return A list with:
#' * `final_data`: full data frame with `CDI`, `CDI_std`, `our_CDI`, annotation.
#' * `pca_object`: the `psych::principal` object (or `NULL` if PCA failed).
#' * `data_with_no_imputation`: the raw, non-imputed input for the scope.
#' * `data_harmonized`: the standardized component data fed into the PCA.
#'
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr filter select any_of mutate across all_of left_join
#' @importFrom dplyr everything bind_rows ntile case_when
#' @importFrom glue glue
#' @importFrom psych principal
#' @importFrom stats sd cor
processCDIScope <- function(
    CDI_data_sub,
    impute_sub = FALSE,
    quiet_sub = FALSE,
    return_format_sub = "detailed",
    weight_var_sub = "tot_pop") {

  # The 18 CDI component variables, in a fixed order.
  varlist <- c(
    "educ_12less_perc", "educ_ba_pl_perc", "emp_white_col_perc",
    "fam_bel_pov_perc", "hhld_1pl_room_perc", "hhld_no_int_perc",
    "hhld_no_veh_perc", "hous_no_plumb_perc", "inc_dis_upd_imp",
    "med_gross_rent_imp", "med_hhld_inc_imp", "med_home_val_imp",
    "med_month_mort_imp", "no_ins_perc", "one_par_hhld_perc_imp",
    "own_occ_hous_perc", "pop_bel_150_pov_perc", "unemp_perc"
  )

  # Hold the input before manipulating.
  raw_CDI_data_sub <- CDI_data_sub |> sf::st_drop_geometry()

  # Filter out no-population areas for this subset.
  CDI_data_sub <- CDI_data_sub |> dplyr::filter(tot_pop > 0)

  rows_with_no_pop <- nrow(raw_CDI_data_sub) - nrow(CDI_data_sub)
  if (!quiet_sub) message(glue::glue("... there are {rows_with_no_pop} geoids with no population."))
  if (nrow(CDI_data_sub) == 0) return(NULL)

  CDI_input_data <- CDI_data_sub |>
    dplyr::select(GEOID, dplyr::all_of(varlist), tot_pop, dplyr::any_of("geometry"))

  # --- Imputation (spatial neighbors) or complete cases ----------------------
  if (impute_sub) {
    if (!"geometry" %in% colnames(CDI_input_data)) {
      stop("Imputation requires geometry. Please re-run with impute = TRUE in computeCDI().")
    }
    if (!quiet_sub) message("... imputing missing data for scope.")

    CDI_input_data <- imputeMissing(
      dfGeo = CDI_input_data,
      weight_var_sub = weight_var_sub,
      quiet = quiet_sub
    )
  } else {
    vars_impute <- setdiff(
      colnames(CDI_input_data), c("GEOID", "geometry", "impute", "has_geometry")
    )
    CDI_input_data <- CDI_input_data |>
      dplyr::mutate(
        nvar_still_missing = rowSums(
          is.na(dplyr::across(dplyr::all_of(vars_impute))), na.rm = FALSE
        ),
        nvar_imputed = 0
      )
  }

  # --- Standardize the 18 components (mean 0, SD 1) across the scope ----------
  CDI_inp_data_fa <- CDI_input_data |>
    dplyr::filter(nvar_still_missing == 0) |>
    dplyr::select(GEOID, tot_pop, dplyr::all_of(varlist))

  CDI_inp_data_fa <- CDI_inp_data_fa |>
    dplyr::mutate(dplyr::across(dplyr::all_of(varlist), ~ scale(.)[, 1], .names = "std_{.col}"))

  std_vars <- paste0("std_", varlist)

  # --- Principal Component Analysis (factor weights) with error handling ------
  current_scope <- "Current Scope"
  if ("state" %in% colnames(raw_CDI_data_sub)) {
    current_scope <- unique(raw_CDI_data_sub$state)[1]
  }

  mat <- as.matrix(CDI_inp_data_fa[, std_vars, drop = FALSE])
  n_obs <- nrow(mat)
  n_vars <- ncol(mat)

  CDI_pca <- NULL
  CDI_score <- rep(NA_real_, n_obs)

  if (n_obs < n_vars + 2) {
    if (!quiet_sub) {
      warning(glue::glue(
        "Skipping PCA for '{current_scope}': ",
        "Insufficient complete observations ({n_obs}) for {n_vars} variables."
      ))
    }
  } else {
    set.seed(5757)
    pca_result <- tryCatch({
      fit <- psych::principal(mat, nfactors = 1, rotate = "none")

      # Factor weights W_j = first-component loadings.
      W <- as.numeric(fit$loadings[, 1])

      # Sign-align so % families below poverty loads positive
      # (higher CDI = more deprivation).
      pov_idx <- match("std_fam_bel_pov_perc", std_vars)
      if (!is.na(pov_idx) && is.finite(W[pov_idx]) && W[pov_idx] != 0) {
        factor_sign <- if (W[pov_idx] < 0) -1 else 1
      } else {
        factor_sign <- if (sum(W < 0) > sum(W > 0)) -1 else 1
      }
      W <- W * factor_sign

      # Factor-weighted composite S_i = sum_j X'_ij * W_j.
      scores <- as.numeric(mat %*% W)

      list(scores = scores, object = fit, weights = W)
    }, error = function(e) {
      if (!quiet_sub) {
        warning(glue::glue(
          "PCA failed for '{current_scope}' which had {n_obs} complete ",
          "observations. Error: {e$message}"
        ))
      }
      return(NULL)
    })

    if (!is.null(pca_result)) {
      CDI_score <- pca_result$scores
      CDI_pca <- pca_result$object
    }
  }

  # --- Rescale to mean 100 / SD 20, then rank into national percentiles ------
  CDI_inp_data_fa$CDI <- CDI_score
  if (nrow(CDI_inp_data_fa) == 0 ||
      all(is.na(CDI_score)) ||
      stats::sd(CDI_score, na.rm = TRUE) == 0) {
    CDI_inp_data_fa$CDI_std <- rep(NA_real_, nrow(CDI_inp_data_fa))
  } else {
    s_mean <- mean(CDI_score, na.rm = TRUE)
    s_sd <- stats::sd(CDI_score, na.rm = TRUE)
    CDI_inp_data_fa$CDI_std <- 100 + 20 * (CDI_score - s_mean) / s_sd
  }

  # --- Assemble output -------------------------------------------------------
  cols_from_raw <- setdiff(colnames(raw_CDI_data_sub), colnames(CDI_input_data))

  CDI_output_data <- CDI_input_data |>
    dplyr::left_join(
      CDI_inp_data_fa |> dplyr::select(GEOID, CDI, CDI_std),
      by = "GEOID"
    ) |>
    dplyr::left_join(
      raw_CDI_data_sub |> dplyr::select(GEOID, dplyr::all_of(cols_from_raw)),
      by = "GEOID"
    ) |>
    dplyr::select(GEOID, dplyr::all_of(cols_from_raw), dplyr::everything())

  # Bring back the no-population geographies.
  CDI_output_data <- CDI_output_data |>
    dplyr::bind_rows(raw_CDI_data_sub |> dplyr::filter(tot_pop == 0))

  # National percentile ranking 1-100 (unweighted; each geography = one unit).
  CDI_output_data <- CDI_output_data |>
    dplyr::mutate(our_CDI = dplyr::ntile(CDI_std, 100))

  # Annotate GEOIDs by data quality / imputation status.
  CDI_output_data <- CDI_output_data |>
    dplyr::mutate(
      annotation = dplyr::case_when(
        tot_pop == 0 ~ "No population",
        !is.na(CDI) & nvar_imputed == 0 & nvar_still_missing == 0 ~ "Complete data",
        !is.na(CDI) & nvar_imputed >= 0 & nvar_still_missing == 0 ~ "Imputation completed",
        is.na(CDI) & nvar_imputed >= 0 & nvar_still_missing > 0 ~ "Imputation incomplete",
        is.na(CDI) & nvar_still_missing == 0 ~ "PCA failed",
        TRUE ~ "Other"
      )
    )

  if (inherits(CDI_output_data, "sf")) {
    CDI_output_data <- sf::st_drop_geometry(CDI_output_data)
  }

  return(list(
    final_data = CDI_output_data,
    pca_object = CDI_pca,
    data_with_no_imputation = raw_CDI_data_sub,
    data_harmonized = CDI_inp_data_fa
  ))
}
