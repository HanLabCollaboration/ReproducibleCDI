#' Rank into groups the way SAS `proc rank groups=` does
#'
#' @description
#' Reproduces `proc rank groups = k` (default `ties = mean`) followed by `+ 1`,
#' as in the CMS reference code (`cdi_construction_s2.sas`):
#' \deqn{group = \lfloor r \cdot k / (n + 1) \rfloor + 1}
#' where \eqn{r} is the rank (average rank for ties) and \eqn{n} the number of
#' non-missing values. This differs slightly from `dplyr::ntile()`.
#'
#' @param x A numeric vector.
#' @param groups Number of groups (100 for national percentiles).
#'
#' @return An integer vector in `1..groups`, `NA` where `x` is missing.
rankGroupsSAS <- function(x, groups = 100) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  n <- sum(!is.na(x))
  as.integer(floor(r * groups / (n + 1)) + 1)
}
