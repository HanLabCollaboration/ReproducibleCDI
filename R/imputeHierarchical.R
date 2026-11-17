#' Hierarchical imputation (CMS)
#'
#' @description
#' Port of the `%impute` macro in the CMS reference code
#' (`cdi_construction_s1.sas`): each variable takes the first non-missing value
#' of the unit itself, then its tract, county, and state. Parent values come
#' from each parent level's own ACS tables.
#'
#' For ratio variables (`inc_dis_upd`), a value of 0 is treated like a missing
#' value (CMS takes the log of the first value that is non-missing *and*
#' non-zero). The log itself is applied afterwards by `finalizeIncDis()`.
#'
#' @param cdi_data A data frame with `GEOID` and the variables in `varlist`.
#' @param parents A named list of parent data frames, in fallback order. Names
#'   must be among `"tract"`, `"county"`, `"state"`; each needs `GEOID` and the
#'   variables.
#' @param varlist Variables to impute.
#' @param ratio_vars Variables for which 0 also triggers the fallback.
#'
#' @return `cdi_data` with the variables filled in, plus `<var>_imp_src`
#'   (`"self"`, `"tract"`, `"county"`, `"state"`, or `NA`), `nvar_imputed`, and
#'   `nvar_still_missing`.
imputeHierarchical <- function(cdi_data, parents, varlist, ratio_vars = "inc_dis_upd") {
  key_len <- c(tract = 11L, county = 5L, state = 2L)
  if (!all(names(parents) %in% names(key_len))) {
    stop("`parents` must be named among: tract, county, state.")
  }

  usable <- function(x, var) {
    if (var %in% ratio_vars) !is.na(x) & x != 0 else !is.na(x)
  }

  for (var in varlist) {
    val <- cdi_data[[var]]
    ok  <- usable(val, var)
    src <- ifelse(ok, "self", NA_character_)
    val[!ok] <- NA_real_

    for (lvl in names(parents)) {
      p  <- parents[[lvl]]
      pv <- p[[var]][match(substr(cdi_data$GEOID, 1, key_len[[lvl]]), p$GEOID)]
      fill <- is.na(src) & usable(pv, var)
      val[fill] <- pv[fill]
      src[fill] <- lvl
    }

    cdi_data[[var]] <- val
    cdi_data[[paste0(var, "_imp_src")]] <- src
  }

  src_mat <- as.matrix(as.data.frame(cdi_data)[, paste0(varlist, "_imp_src"), drop = FALSE])
  cdi_data$nvar_imputed       <- as.integer(rowSums(!is.na(src_mat) & src_mat != "self"))
  cdi_data$nvar_still_missing <- as.integer(rowSums(is.na(src_mat)))
  cdi_data
}
