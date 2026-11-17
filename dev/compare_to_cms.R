# ---------------------------------------------------------------------------
# Validate ReproducibleCDI against the official CMS published index.
#
# Reference: https://github.com/CMSgov/community-deprivation-index
#   CDI_NATRANK.csv = block-group national rank (1-100) produced by the CMS SAS,
#   plus GQ / PH / GQ-PH labels for suppressed block groups.
#
# The CDI is built on 2019 ACS 5-year (2015-2019) block-group data, so we
# compute our index with year = 2019 nationally and compare the 1-100 ranks.
#
# Expect a very high - but not exact - agreement: this package deliberately
# omits CMS's MOE shrinkage and suppression and uses spatial (rather than
# hierarchical) imputation, so some block groups will move a few ranks.
# ---------------------------------------------------------------------------

suppressMessages({
  library(dplyr)
  library(readr)
  devtools::load_all(".")   # or library(ReproducibleCDI)
})

# --- 1. Reference: CMS published ranks --------------------------------------
cms_url <- "https://raw.githubusercontent.com/CMSgov/community-deprivation-index/main/CDI_NATRANK.csv"
cms <- readr::read_csv(cms_url, col_types = "cc") |>          # both as character
  dplyr::rename(GEOID = FIPS_BG, cms_rank_chr = CDI_NATRANK) |>
  dplyr::mutate(
    GEOID      = stringr::str_pad(GEOID, 12, pad = "0"),
    suppressed = cms_rank_chr %in% c("GQ", "PH", "GQ-PH"),
    cms_rank   = suppressWarnings(as.integer(cms_rank_chr))
  )

# --- 2. Our index: national, 2019, block group ------------------------------
# NOTE: this is a large run (~220k block groups, ~51 ACS calls).
#   * impute = FALSE  -> complete cases only, no geometry download (recommended
#                        first pass; much faster / lighter).
#   * impute = TRUE   -> queen-contiguity spatial imputation; downloads geometry
#                        for the whole US (very heavy).

if(!file.exists(file.path(here::here('dev'),'cdi_national_noImp.rds'))){
  ours <- computeCDI(
    geo    = "block group",
    year   = 2019,
    states = "all",
    impute = FALSE,
    scope  = "national",
    return_format = "minimal"   # GEOID, CDI_std, our_CDI, annotation
  ) |>
    dplyr::mutate(GEOID = stringr::str_pad(GEOID, 12, pad = "0"))
  
  saveRDS(
    ours,
    file = file.path(here::here('dev'),'cdi_national_noImp.rds')
  )
} else {
  ours <- readRDS(file = file.path(here::here('dev'),'cdi_national_noImp.rds'))
}

# --- 3. Join and compare (numeric ranks only) -------------------------------
cmp <- ours |>
  dplyr::inner_join(cms, by = "GEOID") |>
  dplyr::filter(!suppressed, !is.na(our_CDI), !is.na(cms_rank))

cat("\nBlock groups compared:", nrow(cmp), "\n")

# Rank agreement
sp <- cor(cmp$our_CDI, cmp$cms_rank, method = "spearman")
pe <- cor(cmp$our_CDI, cmp$cms_rank, method = "pearson")
cat(sprintf("Spearman rank correlation: %.4f\n", sp))
cat(sprintf("Pearson  correlation:      %.4f\n", pe))
cat(sprintf("Mean |rank difference|:    %.2f (of 100)\n", mean(abs(cmp$our_CDI - cmp$cms_rank))))
cat(sprintf("Within +/- 1:  %.1f%%\n", 100 * mean(abs(cmp$our_CDI - cmp$cms_rank) <= 1)))
cat(sprintf("Within +/- 5:  %.1f%%\n", 100 * mean(abs(cmp$our_CDI - cmp$cms_rank) <= 5)))

# Decile agreement matrix (like the paper's Table S1/S2)
cmp <- cmp |>
  dplyr::mutate(
    our_dec = ntile(our_CDI, 10),
    cms_dec = ntile(cms_rank, 10)
  )
cat("\nDecile cross-tab (rows = ours, cols = CMS):\n")
print(table(ours = cmp$our_dec, cms = cmp$cms_dec))

# Optional: save the joined comparison for inspection / plotting
# readr::write_csv(cmp, "dev/cdi_vs_cms_2019.csv")


library(ggplot2)
# install.packages("hexbin")   # geom_hex needs this

sp <- cor(cmp$our_CDI, cmp$cms_rank, method = "spearman")

ggplot(cmp, aes(x = our_CDI, y = cms_rank)) +
  geom_hex(binwidth = c(2, 2)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey30") +
  scale_fill_viridis_c(
    trans = "log10",              # counts pile up on the diagonal
    name = "Block groups",
    labels = scales::comma
  ) +
  coord_equal(xlim = c(1, 100), ylim = c(1, 100)) +
  labs(
    title = "ReproducibleCDI vs. CMS published CDI",
    subtitle = paste0("Block-group national rank, 2019 ACS5  (Spearman = ",
                      sprintf("%.3f", sp), ")"),
    x = "our_CDI (this package)",
    y = "CDI_NATRANK (CMS)"
  ) +
  theme_minimal(base_size = 12)
Key choices, and why:
  

