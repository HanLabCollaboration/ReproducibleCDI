#' Compute the Community Deprivation Index exactly as CMS does
#'
#' @description
#' Reproduces every step of the official CMS Community Deprivation Index
#' reference implementation (`cdi_construction_s1.sas` and
#' `cdi_construction_s2.sas`, CMSgov/community-deprivation-index), at the
#' census block-group level, with no options that change the method:
#' \enumerate{
#'   \item fetch ACS5 estimates at block group (with margins of error), tract,
#'     county, and state;
#'   \item compute the 19 CMS variables at each level (a rate whose denominator
#'     is 0 is set to 0), with block-group MOEs from the CMS formulas;
#'   \item flag block groups for suppression: `PH` (population < 100 or
#'     housing units < 30), `GQ` (>= 1/3 of residents in group quarters), or
#'     `GQ-PH`;
#'   \item stabilize every block group toward its tract with the CMS
#'     MOE-based shrinkage ([stabilizeCDI()]), using all block groups;
#'   \item fill remaining gaps hierarchically from the tract, county, then
#'     state value, and take the log of the income-disparity ratio;
#'   \item standardize and fit a one-component PCA on the non-suppressed block
#'     groups, score them with the scoring coefficients, and rescale to mean
#'     100 / SD 20;
#'   \item rank into 1-100 groups exactly as SAS `proc rank groups=100` does.
#'     Suppressed block groups get their flag (`GQ`, `PH`, `GQ-PH`) instead of
#'     a rank.
#' }
#'
#' Two CMS behaviors are reproduced as-is: a non-suppressed block group whose
#' score is missing receives rank 1 (SAS `sum(missing, 1)`), and the
#' proportion MOE formula is also used for the income-disparity ratio.
#'
#' For the ReproducibleX-style index with spatial imputation and a choice of
#' stabilization order, see [computeCDI()].
#'
#' @param year The ACS 5-year end year. Must be between 2017 and 2022: the
#'   household-internet table B28002 does not exist before ACS5 2017, and the
#'   group-quarters table B09019 needed for the `GQ` suppression is not
#'   published at the block-group level from ACS5 2023 on. CMS published the
#'   index for 2019.
#' @param states A vector of state abbreviations or `"all"` (default; the 50
#'   states, DC, and Puerto Rico, as in the CMS national ranking). Ranks are
#'   relative to the states requested.
#' @param quiet Logical. If `TRUE`, suppresses messages and warnings.
#' @param return_format `"detailed"` (default) or `"minimal"`.
#'
#' @return
#' If `return_format = "minimal"`, a data frame with `GEOID`, `CDI_NATRANK`
#' (character, as in the CMS file: `"1"`-`"100"`, `"GQ"`, `"PH"`, `"GQ-PH"`),
#' and `cdi_rank` (integer, `NA` for suppressed block groups). Otherwise a list:
#' \itemize{
#'   \item `df_cdi`: one row per block group with the 19 final variables
#'     (`inc_dis_upd` on the log scale), `suppression_flag`, `CDI` (raw score),
#'     `CDI_std` (mean 100 / SD 20), `cdi_rank`, `CDI_NATRANK`,
#'     `nvar_stabilized`, `nvar_imputed`, and `annotation`.
#'   \item `data_raw`: block-group values and MOEs before stabilization
#'     (`inc_dis_upd` on the ratio scale).
#'   \item `data_stabilized`: output of [stabilizeCDI()] with all diagnostics.
#'   \item `data_imputed`: values after hierarchical imputation, with
#'     `<var>_imp_src` (`"self"`, `"tract"`, `"county"`, `"state"`).
#'   \item `parents`: tract, county, and state values used.
#'   \item `pca`: scoring `coefficients`, `loadings`, the `psych::principal`
#'     `fit`, and `n_obs`.
#' }
#'
#' @references
#' Robst, J., Forlines, G., Kautter, J., et al. (2024). The development of the
#' Community Deprivation Index and its application to Accountable Care
#' Organizations. \emph{Health Affairs Scholar}. \doi{10.1093/haschl/qxae161}
#'
#' CMS Community Deprivation Index reference implementation,
#' \url{https://github.com/CMSgov/community-deprivation-index}.
#'
#' @examples
#' \dontrun{
#' cdi_2019 <- computeOriginalCDI(year = 2019, states = "all")
#' }
#'
#' @importFrom dplyr left_join select any_of all_of arrange case_when
#' @export
computeOriginalCDI <- function(year, states = "all", quiet = FALSE,
                               return_format = "detailed") {

  # --- 1. Year check (before any Census API call) ----------------------------
  if (!is.numeric(year) || length(year) != 1 || is.na(year)) {
    stop("`year` must be a single ACS5 end year, e.g. 2019.", call. = FALSE)
  }
  if (year > 2022) {
    stop(glue::glue(
      "computeOriginalCDI() supports ACS5 2017-2022 only. For {year}, the Census ",
      "Bureau does not publish table B09019 (group quarters) at the block-group ",
      "level, so the CMS GQ suppression step cannot be reproduced. Use ",
      "computeCDI() for {year}."
    ), call. = FALSE)
  }
  if (year < 2017) {
    stop(glue::glue(
      "computeOriginalCDI() supports ACS5 2017-2022 only. Table B28002 (household ",
      "internet), one of the 18 CDI variables, is not published before ACS5 2017."
    ), call. = FALSE)
  }
  return_format <- rlang::arg_match(return_format, values = c("minimal", "detailed"))
  states <- resolveCDIStates(states, include_pr = TRUE)
  if (!quiet && length(states) < 52) {
    message("Ranks are relative to the ", length(states),
            " state(s) requested; CMS ranks all block groups nationally.")
  }

  vars_all <- cdiVarNames(pca_only = FALSE)
  vars_pca <- cdiVarNames(pca_only = TRUE)

  # --- 2. Fetch block group (with MOE + suppression cells) and parents -------
  if (!quiet) message("Fetching ACS data: block group, tract, county, state...")
  bg_raw <- fetchCDIAcs("block group", year, states, moe = TRUE, suppression = TRUE)
  checkCDIAvailability(bg_raw, "block group", year,
                       cells = c(cdiAcsCells(), unname(cdiSuppressionCells(year))))
  parents <- list(
    tract  = calculateCDIVars(fetchCDIAcs("tract",  year, states), zero_denominator = "zero"),
    county = calculateCDIVars(fetchCDIAcs("county", year, states), zero_denominator = "zero"),
    state  = calculateCDIVars(fetchCDIAcs("state",  year, states), zero_denominator = "zero")
  )

  # --- 3. Variables, MOEs, suppression ---------------------------------------
  if (!quiet) message("Calculating CDI variables and margins of error...")
  bg <- cleanAcsNamesCDI(bg_raw, geo = "block group") |>
    calculateCDIVars(zero_denominator = "zero", with_moe = TRUE) |>
    dplyr::left_join(flagCDISuppression(bg_raw, year), by = "GEOID")

  # --- 4. Stabilization (all block groups) -----------------------------------
  if (!quiet) message("Stabilizing block groups toward their tract...")
  bg_stab <- stabilizeCDI(bg, parents$tract, vars_all, geo_parent_key_len = 11)

  # --- 5. Hierarchical imputation, then log of the income-disparity ratio ----
  if (!quiet) message("Imputing remaining gaps from tract, county, state...")
  bg_values <- bg_stab[, c("GEOID", paste0(vars_all, "_stabilized"))]
  colnames(bg_values) <- c("GEOID", vars_all)
  bg_imp <- imputeHierarchical(bg_values, parents, vars_all) |> finalizeIncDis()

  # --- 6. PCA on non-suppressed block groups ---------------------------------
  universe <- is.na(bg$suppression_flag)
  if (!quiet) message(glue::glue(
    "Fitting the PCA on {sum(universe)} non-suppressed block groups ",
    "({sum(!universe)} suppressed)..."
  ))
  pca <- fitCDIPCA(bg_imp[universe, , drop = FALSE], vars_pca, quiet = quiet,
                   scope_label = "Selected states")

  CDI <- rep(NA_real_, nrow(bg))
  CDI_std <- rep(NA_real_, nrow(bg))
  CDI[universe] <- pca$score
  CDI_std[universe] <- pca$score_std

  # --- 7. SAS proc rank groups=100, + 1 --------------------------------------
  cdi_rank <- rep(NA_integer_, nrow(bg))
  cdi_rank[universe] <- rankGroupsSAS(pca$score, groups = 100)
  # CMS: CDI_NATRANK_NUM = sum(CDI_NATRANK_NUM, 1), so a missing score becomes 1.
  score_missing <- universe & is.na(CDI) & !is.null(pca$fit)
  cdi_rank[score_missing] <- 1L
  CDI_NATRANK <- ifelse(universe, as.character(cdi_rank), bg$suppression_flag)

  # --- 8. Output -------------------------------------------------------------
  flags <- as.matrix(bg_stab[, paste0("stabilizedflg_", vars_pca)])
  src   <- as.matrix(bg_imp[, paste0(vars_pca, "_imp_src")])
  nvar_stabilized <- as.integer(rowSums(flags))
  nvar_imputed    <- as.integer(rowSums(!is.na(src) & src != "self"))

  df_cdi <- bg |>
    dplyr::select(GEOID, dplyr::any_of(c("NAME", "block_group", "tract", "county", "state")),
                  tot_pop, suppression_flag) |>
    dplyr::bind_cols(bg_imp[, vars_all]) |>
    dplyr::mutate(
      CDI = CDI, CDI_std = CDI_std, cdi_rank = cdi_rank, CDI_NATRANK = CDI_NATRANK,
      nvar_stabilized = nvar_stabilized, nvar_imputed = nvar_imputed,
      annotation = dplyr::case_when(
        !universe ~ paste("Suppressed:", suppression_flag),
        is.null(pca$fit) ~ "PCA failed",
        score_missing ~ "Score missing (CMS assigns rank 1)",
        nvar_stabilized > 0 & nvar_imputed > 0 ~ "Stabilized + imputed",
        nvar_imputed > 0 ~ "Imputed",
        nvar_stabilized > 0 ~ "Stabilized",
        TRUE ~ "Complete data"
      )
    ) |>
    dplyr::arrange(GEOID)

  if (return_format == "minimal") {
    return(df_cdi |> dplyr::select(GEOID, CDI_NATRANK, cdi_rank))
  }

  list(
    df_cdi = df_cdi,
    data_raw = bg,
    data_stabilized = bg_stab,
    data_imputed = bg_imp,
    parents = parents,
    pca = pca[c("coefficients", "loadings", "fit", "n_obs")]
  )
}
