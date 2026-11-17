test_that("CDI cells match the CMS SAS &censusvar list cell-for-cell", {
  # Copied from cdi_construction_s1.sas (%let censusvar = ...).
  cms <- strsplit(paste(
    "B15003_002 B15003_003 B15003_004 B15003_005 B15003_006 B15003_007 B15003_008 B15003_009 B15003_010 B15003_011",
    "B15003_012 B15003_013 B15003_014 B15003_015 B15003_016 B15003_001 B15003_017 B15003_018 B15003_019 B15003_020 B15003_021",
    "B15003_022 B15003_023 B15003_024 B15003_025 C24010_003 C24010_039 C24010_027 C24010_063 C24010_001 B19013_001 B19001_002",
    "B19001_003 B19001_004 B19001_014 B19001_015 B19001_016 B19001_017 B25077_001 B25064_001 B25088_002 B25003_002 B25003_001",
    "B23025_005 B23025_002 B17010_002 B17010_001 C17002_002 C17002_003 C17002_004 C17002_005 C17002_001 B11003_010 B11003_016",
    "B11003_001 B25044_003 B25044_010 B25044_001 B28002_013 B28002_001 B25016_007 B25016_016 B25016_001 B25014_005 B25014_006",
    "B25014_007 B25014_011 B25014_012 B25014_013 B25014_001 B27010_017 B27010_033 B27010_050 B27010_066 B27010_001"
  ), " ")[[1]]
  expect_setequal(ReproducibleCDI:::cdiAcsCells(), cms)
})

test_that("18 PCA variables plus educ_hs_pl_perc", {
  expect_length(ReproducibleCDI:::cdiVarNames(TRUE), 18)
  expect_length(ReproducibleCDI:::cdiVarNames(FALSE), 19)
  expect_false("educ_hs_pl_perc" %in% ReproducibleCDI:::cdiVarNames(TRUE))
})

test_that("group-quarters cell depends on the year", {
  expect_equal(ReproducibleCDI:::cdiSuppressionCells(2018)[["gq"]], "B09019_038")
  expect_equal(ReproducibleCDI:::cdiSuppressionCells(2019)[["gq"]], "B09019_026")
})
