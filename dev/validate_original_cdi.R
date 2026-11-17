# ---------------------------------------------------------------------------
# Validate computeOriginalCDI() against the official CMS published index.
#
# Reference: https://github.com/CMSgov/community-deprivation-index
#   CDI_NATRANK.csv = block-group national rank (1-100) from the CMS SAS code,
#   with GQ / PH / GQ-PH labels for suppressed block groups (50 states, DC, PR).
#
# computeOriginalCDI() reproduces every CMS step, so the 2019 ranks should
# match closely. Output: dev/validate_original_cdi_2019.txt and an .rds cache.
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr)
  library(readr)
  devtools::load_all(here::here(), quiet = TRUE)
})

out_txt <- here::here("dev", "validate_original_cdi_2019.txt")
out_rds <- here::here("dev", "cdi_original_2019.rds")

# --- 1. Reference: CMS published ranks --------------------------------------
cms_url <- "https://raw.githubusercontent.com/CMSgov/community-deprivation-index/main/CDI_NATRANK.csv"
cms <- readr::read_csv(cms_url, col_types = "cc") |>
  dplyr::rename(GEOID = FIPS_BG, cms_label = CDI_NATRANK) |>
  dplyr::mutate(
    GEOID = stringr::str_pad(GEOID, 12, pad = "0"),
    cms_suppressed = cms_label %in% c("GQ", "PH", "GQ-PH"),
    cms_rank = suppressWarnings(as.integer(cms_label))
  )

# --- 2. Ours: computeOriginalCDI(2019, "all") -------------------------------
if (!file.exists(out_rds)) {
  t0 <- Sys.time()
  res <- computeOriginalCDI(year = 2019, states = "all", quiet = FALSE)
  elapsed <- Sys.time() - t0
  saveRDS(list(df_cdi = res$df_cdi, pca = res$pca[c("coefficients", "loadings", "n_obs")],
               imp_src = res$data_imputed[, c("GEOID", grep("_imp_src$", names(res$data_imputed), value = TRUE))],
               elapsed = elapsed),
          out_rds)
}
ours <- readRDS(out_rds)

sink(out_txt, split = TRUE)
cat("computeOriginalCDI(2019, 'all') vs CMS CDI_NATRANK.csv\n")
cat("Run time:", format(ours$elapsed), "\n\n")

cmp <- ours$df_cdi |>
  dplyr::select(GEOID, suppression_flag, CDI_NATRANK, cdi_rank, annotation) |>
  dplyr::full_join(cms, by = "GEOID")

cat("Block groups - ours:", nrow(ours$df_cdi), " CMS:", nrow(cms),
    " in both:", sum(!is.na(cmp$CDI_NATRANK) & !is.na(cmp$cms_label)), "\n\n")

# --- 3. Suppression agreement -----------------------------------------------
cat("Suppression label (rows = ours, cols = CMS):\n")
lab <- function(x) ifelse(x %in% c("GQ", "PH", "GQ-PH"), x, ifelse(is.na(x), "<absent>", "ranked"))
print(table(ours = lab(cmp$CDI_NATRANK), cms = lab(cmp$cms_label), useNA = "ifany"))

# --- 4. Rank agreement (both ranked) ----------------------------------------
both <- cmp |> dplyr::filter(!is.na(cdi_rank), !is.na(cms_rank))
d <- both$cdi_rank - both$cms_rank
cat("\nBoth ranked:", nrow(both), "\n")
cat(sprintf("Exact match:   %.2f%%\n", 100 * mean(d == 0)))
cat(sprintf("Within +/- 1:  %.2f%%\n", 100 * mean(abs(d) <= 1)))
cat(sprintf("Within +/- 5:  %.2f%%\n", 100 * mean(abs(d) <= 5)))
cat(sprintf("Mean |diff|:   %.3f\n", mean(abs(d))))
cat(sprintf("Spearman:      %.5f\n", cor(both$cdi_rank, both$cms_rank, method = "spearman")))
cat("\nRank difference (ours - CMS), distribution:\n")
print(table(cut(d, c(-Inf, -10, -5, -2, -1, 0, 1, 2, 5, 10, Inf), right = TRUE)))

cat("\nExact match by annotation:\n")
print(both |> dplyr::mutate(exact = cdi_rank == cms_rank) |>
        dplyr::group_by(annotation) |>
        dplyr::summarise(n = dplyr::n(), exact_pct = round(100 * mean(exact), 2), .groups = "drop"))

cat("\nWorst states by exact match:\n")
print(both |> dplyr::mutate(st = substr(GEOID, 1, 2), exact = cdi_rank == cms_rank) |>
        dplyr::group_by(st) |>
        dplyr::summarise(n = dplyr::n(), exact_pct = round(100 * mean(exact), 2),
                         mean_abs = round(mean(abs(cdi_rank - cms_rank)), 3), .groups = "drop") |>
        dplyr::arrange(exact_pct) |> head(10))

cat("\nPCA: n_obs =", ours$pca$n_obs, "\nScoring coefficients:\n")
print(round(ours$pca$coefficients, 5))

# --- 5. Previous package version (no stabilization / imputation) ------------
old_rds <- here::here("dev", "cdi_national_noImp.rds")
if (file.exists(old_rds)) {
  old <- readRDS(old_rds) |>
    dplyr::inner_join(cms, by = "GEOID") |>
    dplyr::filter(!is.na(our_CDI), !is.na(cms_rank))
  cat(sprintf("\nFor reference, previous computeCDI (no shrinkage/imputation/suppression): exact %.2f%%, within +/-1 %.2f%%, Spearman %.4f (n = %d)\n",
              100 * mean(old$our_CDI == old$cms_rank), 100 * mean(abs(old$our_CDI - old$cms_rank) <= 1),
              cor(old$our_CDI, old$cms_rank, method = "spearman"), nrow(old)))
}
sink()
