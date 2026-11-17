#' Stabilize CDI Variables Using Parent-Level Information
#'
#' Applies the empirical Bayes-style shrinkage of the CMS Community Deprivation
#' Index reference code (`%shrink` macro in `cdi_construction_s1.sas`): each
#' lower-level estimate (e.g., block group) is pulled toward its parent-level
#' estimate (e.g., tract) with a weight that compares the estimate's sampling
#' variance (from its margin of error) to the heterogeneity of lower-level
#' units within the parent. The interface and output columns mirror
#' `ReproducibleYostIndex::stabilizeYostIndex()`.
#'
#' @param cdi_data A data frame of lower-level units with a `GEOID` column, the
#'   variables in `varlist`, and their margins of error named `<var>_moe`.
#' @param cdi_data_parent A data frame of parent-level units with a `GEOID`
#'   column and the same variables. Parent values come from the parent level's
#'   own ACS tables, not from aggregating the lower-level units.
#' @param varlist A character vector of variable names to stabilize.
#' @param geo_parent_key_len Integer. Number of leading `GEOID` characters that
#'   identify the parent (11 = tract, 5 = county, 2 = state).
#'
#' @details
#' For variable \eqn{v} and lower-level unit \eqn{i} in parent \eqn{T}:
#' \enumerate{
#'   \item \eqn{n_T} is the number of rows of `cdi_data` in \eqn{T}. Every row
#'     counts, including rows where \eqn{v} is missing (as in the CMS and Yost
#'     code). Which rows are passed in is decided by the caller.
#'   \item Within-parent heterogeneity:
#'   \deqn{t^2 = \frac{\sum_{i \in T} (v_i - v_T)^2}{n_T - 1}}
#'   summing over rows where both \eqn{v_i} and \eqn{v_T} are observed.
#'   \item Sampling variance from the 90\% MOE: \eqn{s^2 = (MOE_i / 1.645)^2}.
#'   \item Weight: \eqn{w = \frac{1/s^2}{1/s^2 + 1/t^2}}.
#'   \item Stabilized estimate: \eqn{v_i^{*} = w v_i + (1 - w) v_T}.
#' }
#' The original value is kept (not stabilized) when \eqn{v_i}, \eqn{v_T}, or
#' the MOE is missing, when \eqn{n_T = 1}, or when \eqn{s^2 \le 0} or
#' \eqn{t^2 \le 0}.
#'
#' @return `cdi_data` (same rows, same order; geometry kept if present) plus:
#' \itemize{
#'   \item `parent_geoid`, `parent_n`: the parent `GEOID` and \eqn{n_T}.
#'   \item `<var>_parent`, `<var>_sqdiff`, `<var>_t2`, `<var>_s2`, `<var>_wgt`
#'     (weight; `NA` when not stabilized).
#'   \item `<var>_stabilized`: the stabilized estimate.
#'   \item `stabilizedflg_<var>`: 1 when shrinkage was applied, 0 otherwise.
#'   \item `nvar_stabilized`: number of variables stabilized in the row.
#' }
#'
#' @references
#' CMS Community Deprivation Index reference implementation,
#' \url{https://github.com/CMSgov/community-deprivation-index}.
#'
#' @examples
#' \dontrun{
#' stabilized <- stabilizeCDI(
#'   cdi_data = bg_data,
#'   cdi_data_parent = tract_data,
#'   varlist = c("fam_bel_pov_perc", "med_hhld_inc"),
#'   geo_parent_key_len = 11
#' )
#' }
#'
#' @export
stabilizeCDI <- function(cdi_data, cdi_data_parent, varlist, geo_parent_key_len) {

  missing_moe <- setdiff(paste0(varlist, "_moe"), colnames(cdi_data))
  if (length(missing_moe) > 0) {
    stop("Missing MOE column(s): ", paste(missing_moe, collapse = ", "))
  }

  geometry <- if (inherits(cdi_data, "sf")) sf::st_geometry(cdi_data) else NULL
  out <- if (is.null(geometry)) cdi_data else sf::st_drop_geometry(cdi_data)
  parent <- if (inherits(cdi_data_parent, "sf")) sf::st_drop_geometry(cdi_data_parent) else cdi_data_parent

  out$parent_geoid <- substr(out$GEOID, 1, geo_parent_key_len)
  parent_idx <- match(out$parent_geoid, parent$GEOID)
  grp <- factor(out$parent_geoid)
  out$parent_n <- as.integer(table(grp)[as.character(grp)])

  for (var in varlist) {
    v  <- out[[var]]
    vp <- parent[[var]][parent_idx]
    moe <- out[[paste0(var, "_moe")]]

    # Sum of squared differences within the parent, over rows where both the
    # unit and the parent value are observed (missing if there are none).
    d <- (v - vp)^2
    has_d <- !is.na(d)
    s <- rowsum(ifelse(has_d, d, 0), grp, reorder = FALSE)[as.character(grp), 1]
    k <- rowsum(as.integer(has_d), grp, reorder = FALSE)[as.character(grp), 1]
    sqdiff <- ifelse(k > 0, s, NA_real_)

    t2 <- ifelse(out$parent_n > 1, sqdiff / (out$parent_n - 1), NA_real_)
    s2 <- (moe / 1.645)^2

    applied <- !is.na(v) & !is.na(vp) & !is.na(moe) & out$parent_n > 1 &
      !is.na(s2) & !is.na(t2) & s2 > 0 & t2 > 0
    w <- ifelse(applied, (1 / s2) / ((1 / s2) + (1 / t2)), NA_real_)

    out[[paste0(var, "_parent")]]     <- vp
    out[[paste0(var, "_sqdiff")]]     <- unname(sqdiff)
    out[[paste0(var, "_t2")]]         <- unname(t2)
    out[[paste0(var, "_s2")]]         <- s2
    out[[paste0(var, "_wgt")]]        <- unname(w)
    out[[paste0(var, "_stabilized")]] <- unname(ifelse(applied, w * v + (1 - w) * vp, v))
    out[[paste0("stabilizedflg_", var)]] <- as.integer(applied)
  }

  out$nvar_stabilized <- as.integer(rowSums(
    as.matrix(out[, paste0("stabilizedflg_", varlist), drop = FALSE])
  ))

  if (!is.null(geometry)) out <- sf::st_set_geometry(out, geometry)
  out
}
