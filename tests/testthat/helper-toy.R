# Toy wide ACS data (estimates + MOEs) for every CDI cell. All counts are set
# from `base`; callers overwrite individual cells to create edge cases.
toy_acs <- function(n = 3, base = 100, moe = 10, geoid = NULL) {
  cells <- unique(c("B01001_001", ReproducibleCDI:::cdiAcsCells()))
  df <- data.frame(
    GEOID = if (is.null(geoid)) sprintf("440010301%03d", seq_len(n)) else geoid,
    NAME = "toy",
    stringsAsFactors = FALSE
  )
  for (cl in cells) {
    df[[paste0(cl, "E")]] <- base
    df[[paste0(cl, "M")]] <- moe
  }
  df
}

# A square grid of polygons as an sf object (for spatial imputation tests).
toy_grid <- function(nx, ny) {
  polys <- list()
  for (j in seq_len(ny)) for (i in seq_len(nx)) {
    polys[[length(polys) + 1]] <- sf::st_polygon(list(rbind(
      c(i - 1, j - 1), c(i, j - 1), c(i, j), c(i - 1, j), c(i - 1, j - 1)
    )))
  }
  sf::st_sfc(polys)
}
