#' Impute Missing Values using Spatial Neighbors (Weighted or Unweighted)
#'
#' @description
#' A highly optimized spatial imputation function that fills missing values using
#' the mean or weighted mean of contiguous neighbors (Queen contiguity).
#'
#' @param dfGeo An `sf` data frame. Must contain a `GEOID` column and valid geometries.
#' @param weight_var_sub A string. The name of the column to use for weighting
#'   (e.g., "tot_pop"). If set to `"none"` (default), `NULL`, or `""`, the function
#'   performs an unweighted arithmetic mean.
#' @param quiet Logical. If `FALSE` (default), prints progress messages to the console.
#'
#' @return A `tibble` (without geometry) containing the original data with missing
#'   values imputed. Includes two metadata columns: `nvar_imputed` and `nvar_still_missing`.
#'
#' @importFrom sf st_is_empty st_drop_geometry st_as_sf
#' @importFrom spdep poly2nb
#' @importFrom dplyr bind_rows select mutate all_of left_join
#' @importFrom glue glue
#' @importFrom stats weighted.mean
#' @export
imputeMissing <- function(dfGeo, weight_var_sub = "none", quiet = FALSE) {

  # --- 1. Validation & Setup ---
  if (!all(c("GEOID", "geometry") %in% colnames(dfGeo))) {
    stop("Data frame must contain 'GEOID' and 'geometry' columns.")
  }

  if (!inherits(dfGeo, "sf")) dfGeo <- sf::st_as_sf(dfGeo)

  # Setup Weights Logic based on the weight_var_sub parameter
  use_weights <- !is.null(weight_var_sub) && weight_var_sub != "none" && weight_var_sub != ""

  if (use_weights) {
    if (!weight_var_sub %in% colnames(dfGeo)) {
      stop(glue::glue("Weight variable '{weight_var_sub}' not found in data."))
    }
  }

  # Identify variables to impute
  meta_cols <- c("GEOID", "geometry", "nvar_imputed", "nvar_still_missing", "county_state")
  if (use_weights) meta_cols <- c(meta_cols, weight_var_sub)

  vars_impute <- setdiff(colnames(dfGeo), meta_cols)

  # --- 2. Split Data: Valid Geometry vs. Invalid ---
  has_geom_idx <- !sf::st_is_empty(dfGeo)

  df_process_sf <- dfGeo[has_geom_idx, ]
  df_skipped_sf <- dfGeo[!has_geom_idx, ]

  n_total <- nrow(dfGeo)
  n_process <- nrow(df_process_sf)

  # --- 3. Build Neighbors ---
  if (n_process == 0) {
    if (!quiet) message("No valid geometries found. Skipping imputation.")
    df_out <- dfGeo |>
      sf::st_drop_geometry() |>
      dplyr::mutate(
        nvar_imputed = 0,
        nvar_still_missing = rowSums(is.na(dplyr::select(dfGeo, dplyr::all_of(vars_impute))))
      )
    return(df_out)
  }

  if (!quiet) {
    message(glue::glue("Processing {n_process} rows with valid geometry out of {n_total} total."))
    if (use_weights) {
      message(glue::glue("Using '{weight_var_sub}' for weighted spatial averaging."))
    } else {
      message("Using unweighted spatial averaging.")
    }
  }

  # Try to build neighbors with error handling
  nb <- tryCatch({
    spdep::poly2nb(df_process_sf, queen = TRUE)
  }, error = function(e) {
    if (!quiet) {
      warning(glue::glue("Failed to build spatial neighbors: {e$message}. Skipping imputation."))
    }
    return(NULL)
  })

  # If neighbor building failed, return data without imputation
  if (is.null(nb)) {
    df_out <- dfGeo |>
      sf::st_drop_geometry() |>
      dplyr::mutate(
        nvar_imputed = 0,
        nvar_still_missing = rowSums(is.na(dplyr::select(dfGeo, dplyr::all_of(vars_impute))))
      )
    return(df_out)
  }

  # Drop geometry for speed
  dat_plain <- sf::st_drop_geometry(df_process_sf)

  # Extract weights vector once if needed
  if (use_weights) {
    weights_all <- dat_plain[[weight_var_sub]]
  }

  # Pre-calculate missingness matrix
  na_mat_pre <- is.na(dat_plain[, vars_impute, drop = FALSE])

  # --- 4. The Optimized Imputation Loop (Column-wise) ---
  for (var in vars_impute) {

    na_indices <- which(na_mat_pre[, var])
    if (length(na_indices) == 0) next

    col_values <- dat_plain[[var]]

    imputed_vals <- vapply(na_indices, function(i) {
      neighbor_ids <- nb[[i]]

      # Handle case where spdep returns 0 (no neighbors)
      if (length(neighbor_ids) == 0 || (length(neighbor_ids) == 1 && neighbor_ids[1] == 0)) {
        return(NA_real_)
      }

      val_neighbors <- col_values[neighbor_ids]

      # Toggle between weighted and unweighted
      if (use_weights) {
        w_neighbors <- weights_all[neighbor_ids]
        res <- stats::weighted.mean(val_neighbors, w = w_neighbors, na.rm = TRUE)
      } else {
        res <- mean(val_neighbors, na.rm = TRUE)
      }

      # Handle cases where all neighbors were NA (result is NaN)
      if (is.nan(res)) return(NA_real_)

      return(res)

    }, FUN.VALUE = numeric(1))

    dat_plain[[var]][na_indices] <- imputed_vals
  }

  # --- 5. Statistics & Reassemble ---
  na_mat_post <- is.na(dat_plain[, vars_impute, drop = FALSE])

  dat_plain$nvar_still_missing <- rowSums(na_mat_post)
  dat_plain$nvar_imputed <- rowSums(na_mat_pre) - dat_plain$nvar_still_missing

  # Handle the skipped rows (those with no geometry)
  dat_skipped <- sf::st_drop_geometry(df_skipped_sf)
  if (nrow(dat_skipped) > 0) {
    dat_skipped$nvar_imputed <- 0
    dat_skipped$nvar_still_missing <- rowSums(is.na(dat_skipped[, vars_impute, drop = FALSE]))
    df_final <- dplyr::bind_rows(dat_plain, dat_skipped)
  } else {
    df_final <- dat_plain
  }

  return(df_final)
}
