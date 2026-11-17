#' Fit the CDI principal component and score
#'
#' @description
#' Shared PCA step of `computeOriginalCDI()` and `computeCDI()`, following
#' `cdi_construction_s2.sas`:
#' \enumerate{
#'   \item standardize each variable to mean 0 / SD 1 over the rows supplied,
#'     using every available value of that variable (`proc standard`);
#'   \item fit a one-component principal component analysis on the correlation
#'     matrix of the complete rows (`proc factor method=principal`);
#'   \item score each row with the standardized scoring coefficients
#'     (loadings / eigenvalue), sign-aligned so that `fam_bel_pov_perc` loads
#'     positive (higher = more deprived);
#'   \item rescale the score to mean 100 / SD 20.
#' }
#'
#' @param data A plain data frame of the rows forming the PCA universe.
#' @param varlist The PCA variables.
#' @param quiet Logical, suppress warnings.
#' @param scope_label Label used in warnings.
#'
#' @return A list with `score` (raw composite, aligned with the rows of
#'   `data`), `score_std` (mean 100 / SD 20), `coefficients`, `loadings`,
#'   `fit` (the `psych::principal` object or `NULL`), and `n_obs` (complete
#'   rows used in the PCA).
#'
#' @importFrom psych principal
#' @importFrom stats sd complete.cases
fitCDIPCA <- function(data, varlist, quiet = FALSE, scope_label = "Current Scope") {
  x <- as.matrix(as.data.frame(data)[, varlist, drop = FALSE])
  z <- apply(x, 2, function(col) (col - mean(col, na.rm = TRUE)) / stats::sd(col, na.rm = TRUE))
  z <- matrix(z, nrow = nrow(x), dimnames = list(NULL, varlist))

  complete <- stats::complete.cases(z)
  n_obs <- sum(complete)
  n_vars <- length(varlist)

  empty <- list(
    score = rep(NA_real_, nrow(x)), score_std = rep(NA_real_, nrow(x)),
    coefficients = NULL, loadings = NULL, fit = NULL, n_obs = n_obs
  )

  if (n_obs < n_vars + 2) {
    if (!quiet) {
      warning(glue::glue(
        "Skipping PCA for '{scope_label}': ",
        "Insufficient complete observations ({n_obs}) for {n_vars} variables."
      ))
    }
    return(empty)
  }

  set.seed(5757)
  fit <- tryCatch(
    psych::principal(z[complete, , drop = FALSE], nfactors = 1, rotate = "none"),
    error = function(e) {
      if (!quiet) {
        warning(glue::glue(
          "PCA failed for '{scope_label}' which had {n_obs} complete ",
          "observations. Error: {e$message}"
        ))
      }
      NULL
    }
  )
  if (is.null(fit)) return(empty)

  loadings <- stats::setNames(as.numeric(fit$loadings[, 1]), varlist)
  coefficients <- loadings / fit$values[1]

  # Sign-align so % families below poverty loads positive.
  pov <- loadings[["fam_bel_pov_perc"]]
  factor_sign <- if (!is.na(pov) && pov != 0) sign(pov) else if (sum(loadings < 0) > sum(loadings > 0)) -1 else 1
  loadings <- loadings * factor_sign
  coefficients <- coefficients * factor_sign

  score <- as.numeric(z %*% coefficients)
  s_sd <- stats::sd(score, na.rm = TRUE)
  score_std <- if (is.na(s_sd) || s_sd == 0) {
    rep(NA_real_, length(score))
  } else {
    100 + 20 * (score - mean(score, na.rm = TRUE)) / s_sd
  }

  list(
    score = score, score_std = score_std, coefficients = coefficients,
    loadings = loadings, fit = fit, n_obs = n_obs
  )
}
