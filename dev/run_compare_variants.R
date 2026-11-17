# ---------------------------------------------------------------------------
# Precompute the national 2019 block-group runs used by compare_to_cms.qmd.
#
#   cdi_original_2019.rds       computeOriginalCDI(2019, "all")  (CMS method)
#   cdi_variant_<name>_2019.rds computeCDI(...) variants (ReproducibleX method)
#
# Each run is cached; delete an .rds to recompute it. The three imputation
# variants download block-group geometry for the whole US and build the
# queen-contiguity neighbours nationally, so they are slow.
# ---------------------------------------------------------------------------

suppressMessages(devtools::load_all(here::here(), quiet = TRUE))
options(tigris_use_cache = TRUE)

dev_dir <- if (basename(here::here()) == "dev") here::here() else here::here("dev")

# --- CMS method -------------------------------------------------------------
original_rds <- file.path(dev_dir, "cdi_original_2019.rds")
if (!file.exists(original_rds)) {
  t0 <- Sys.time()
  res <- computeOriginalCDI(year = 2019, states = "all", quiet = TRUE)
  saveRDS(list(
    df_cdi = res$df_cdi,
    pca = res$pca[c("coefficients", "loadings", "n_obs")],
    imp_src = res$data_imputed[, c("GEOID", grep("_imp_src$", names(res$data_imputed), value = TRUE))],
    elapsed = Sys.time() - t0
  ), original_rds)
  rm(res); gc()
}

# --- ReproducibleX variants -------------------------------------------------
variants <- list(
  complete_case   = list(stabilize = FALSE, impute = FALSE, order = "stabilize_first"),
  stabilize_only  = list(stabilize = TRUE,  impute = FALSE, order = "stabilize_first"),
  impute_only     = list(stabilize = FALSE, impute = TRUE,  order = "stabilize_first"),
  stabilize_first = list(stabilize = TRUE,  impute = TRUE,  order = "stabilize_first"),
  impute_first    = list(stabilize = TRUE,  impute = TRUE,  order = "impute_first")
)

for (nm in names(variants)) {
  rds <- file.path(dev_dir, paste0("cdi_variant_", nm, "_2019.rds"))
  if (file.exists(rds)) next
  v <- variants[[nm]]
  message(format(Sys.time(), "%H:%M"), "  running ", nm, " ...")
  t0 <- Sys.time()
  res <- computeCDI(
    geo = "block group", year = 2019, states = "all",
    stabilize = v$stabilize, impute = v$impute, order = v$order,
    scope = "national", return_format = "minimal", quiet = TRUE
  )
  attr(res, "elapsed") <- Sys.time() - t0
  attr(res, "settings") <- v
  saveRDS(res, rds)
  message(format(Sys.time(), "%H:%M"), "  done ", nm, " (", format(attr(res, "elapsed")), ")")
  rm(res); gc()
}
