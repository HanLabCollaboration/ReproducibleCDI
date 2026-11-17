#' SAS-style row sum of ACS cells
#'
#' @description
#' Mirrors SAS `sum(of ...)`: missing cells are skipped, but the result is
#' missing when every cell is missing.
#'
#' @param data A plain data frame.
#' @param cols Column names to sum.
#' @return A numeric vector.
#' @noRd
sasRowSum <- function(data, cols) {
  if (length(cols) == 0) return(rep(NA_real_, nrow(data)))
  m <- as.matrix(data[, cols, drop = FALSE])
  out <- rowSums(m, na.rm = TRUE)
  out[rowSums(!is.na(m)) == 0] <- NA_real_
  out
}

#' Margin of error of a derived proportion or ratio (CMS `moe2` / `moe3`)
#'
#' @description
#' Port of the `%moe2` / `%moe3` macros in `cdi_construction_s1.sas`. With
#' \eqn{p = num / den},
#' \deqn{x = MOE_{num}^2 - p^2 MOE_{den}^2}
#' and the MOE is \eqn{\sqrt{x} / den} when \eqn{x > 0}, otherwise
#' \eqn{\sqrt{MOE_{num}^2 + p^2 MOE_{den}^2} / den}. It is missing when the
#' denominator is zero or missing. The "subtract" (proportion) formula is used
#' for every variable, including the income-disparity ratio, exactly as CMS does.
#'
#' `num_moe` / `den_moe` are the already-combined MOEs of the numerator and
#' denominator (see [cdiCombinedMoe()]).
#'
#' @param num,den Numerator and denominator estimates.
#' @param num_moe,den_moe Their margins of error.
#' @return A numeric vector of MOEs.
#' @noRd
cdiMoeRatio <- function(num, den, num_moe, den_moe) {
  ok <- !is.na(den) & den != 0
  p <- num / den
  x_sub <- num_moe^2 - p^2 * den_moe^2
  x_add <- num_moe^2 + p^2 * den_moe^2
  # SAS treats a missing x_sub as "not > 0", so it falls through to x_add.
  root <- ifelse(!is.na(x_sub) & x_sub > 0, sqrt(pmax(x_sub, 0)), sqrt(x_add))
  ifelse(ok, root / den, NA_real_)
}

#' Combined MOE of a sum of ACS cells
#'
#' @description
#' Root sum of squares of the cell MOEs (SAS `sqrt(sum(m1**2, m2**2, ...))`),
#' with no special handling of zero-estimate cells, as in the CMS code. A single
#' cell returns its own MOE.
#'
#' @param data A plain data frame holding the `<cell>M` columns.
#' @param cells ACS cells (no suffix).
#' @return A numeric vector.
#' @noRd
cdiCombinedMoe <- function(data, cells) {
  cols <- paste0(cells, "M")
  if (length(cols) == 1) return(as.numeric(data[[cols]]))
  sq <- as.data.frame(lapply(data[, cols, drop = FALSE], function(x) x^2))
  sqrt(sasRowSum(sq, names(sq)))
}
