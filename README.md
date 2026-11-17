
<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/HanLabCollaboration/ReproducibleCDI/graph/badge.svg)](https://app.codecov.io/gh/HanLabCollaboration/ReproducibleCDI)
<!-- badges: end -->

# ReproducibleCDI

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The goal of `ReproducibleCDI` is to provide a simple, tidy interface for
calculating the **Community Deprivation Index (CDI)**, a composite
measure of neighborhood-level socioeconomic deprivation, using data from
the US Census Bureau’s American Community Survey (ACS).

## What is the Community Deprivation Index?

The Community Deprivation Index (Robst et al., 2024) follows the Singh /
Area Deprivation Index (ADI) methodology and updates it: it re-estimates
the factor weights on recent ACS data, modernizes several inputs (e.g.,
household internet in place of telephone service), and **adds a
percentage-uninsured variable**, for a total of **18 standardized
variables**. It is designed as a **census block group** index.

The variables in this package are **exactly** those of the official CMS
reference implementation
([CMSgov/community-deprivation-index](https://github.com/CMSgov/community-deprivation-index),
SAS): same ACS cells, same formulas. No other table is ever substituted;
if a CMS table is not published at a geography for a given year, the
package stops and says so.

The 18 component variables (CMS names, lowercase) are:

- **`educ_12less_perc`** — population aged 25+ with less than 12 years
  of schooling (no diploma)
- **`educ_ba_pl_perc`** — population aged 25+ with a bachelor’s degree
  or higher
- **`emp_white_col_perc`** — civilian employed 16+ in white-collar
  occupations (management/business/science/arts plus sales & office)
- **`fam_bel_pov_perc`** — families with income below the poverty line
- **`hhld_1pl_room_perc`** — occupied units with more than one occupant
  per room (crowding)
- **`hhld_no_int_perc`** — households with no internet access
- **`hhld_no_veh_perc`** — households with no vehicle available
- **`hous_no_plumb_perc`** — housing units lacking complete plumbing
- **`inc_dis_upd`** — income disparity: natural log of 100 × (households
  under \$20k) / (households over \$100k)
- **`med_gross_rent`** — median gross rent
- **`med_hhld_inc`** — median household income
- **`med_home_val`** — median home value
- **`med_month_mort`** — median monthly owner cost with a mortgage
- **`no_ins_perc`** — population with no health insurance coverage
- **`one_par_hhld_perc`** — single-parent family households with own
  children under 18
- **`own_occ_hous_perc`** — owner-occupied housing units
- **`pop_bel_150_pov_perc`** — population below 150% of the poverty line
- **`unemp_perc`** — unemployed share of the civilian labor force 16+

CMS also computes **`educ_hs_pl_perc`** (high-school diploma or more)
but does not use it in the index; the package returns it for
completeness. Values are used exactly as `tidycensus` returns them
(Census jam codes are already missing; top-coded medians such as
\$250,001 stay numeric).

Higher CDI scores indicate greater deprivation.

## Two ways to compute the CDI

|  | `computeOriginalCDI()` | `computeCDI()` |
|----|----|----|
| Purpose | Reproduce **every** CMS step | CDI built the ReproducibleX way |
| Geography | block group | county, tract, block group |
| Years (ACS5) | 2017–2022 | 2017 and later |
| Rate with a zero denominator | 0 (CMS) | missing (as in the sibling packages) |
| Stabilization | CMS shrinkage toward the tract | optional (`stabilize = TRUE`), same function |
| Imputation | hierarchical: tract → county → state | optional spatial queen-contiguity (`impute = TRUE`) |
| Order of the two | stabilize, then impute | your choice: `order = "stabilize_first"` or `"impute_first"` |
| Suppression | `GQ` / `PH` / `GQ-PH` block groups | only zero-population geographies dropped |
| Ranking | SAS `proc rank groups=100` | `dplyr::ntile(…, 100)` |

### `computeOriginalCDI()`: the CMS method, step by step

1.  Fetch ACS5 data at block group (with margins of error), tract,
    county, and state. Every variable is computed at each level from
    that level’s own tables.
2.  Compute the 19 CMS variables; block-group margins of error follow
    the CMS formulas.
3.  Flag block groups for suppression: `PH` (population \< 100 or
    housing units \< 30), `GQ` (≥ 1/3 of residents in group quarters),
    or `GQ-PH`.
4.  **Stabilize** each block-group value toward its tract with the CMS
    empirical Bayes shrinkage (`stabilizeCDI()`):
    $w = \frac{1/s^2}{1/s^2 + 1/t^2}$, where $s^2 = (MOE/1.645)^2$ and
    $t^2$ is the spread of the tract’s block groups around the tract
    value.
5.  **Impute** remaining gaps hierarchically from the tract, county,
    then state value, and take the log of the income-disparity ratio.
6.  Standardize, fit a one-component PCA on the non-suppressed block
    groups, score with the scoring coefficients, rescale to mean 100 /
    SD 20.
7.  Rank into 1–100 exactly as SAS `proc rank groups=100`; suppressed
    block groups get their flag instead of a rank (`CDI_NATRANK`).

The national ranking (`states = "all"`) includes the 50 states, DC, and
Puerto Rico, as in the CMS file. Only ACS5 2017–2022 are supported: the
household-internet table (B28002) starts in 2017, and the group-quarters
table (B09019) needed for `GQ` suppression is not published at block
group from ACS5 2023 on. Asking for a later year is an error.

### `computeCDI()`: ReproducibleX conventions, with your choice of order

`computeCDI()` uses the same 18 variables, sets a rate with a zero
denominator to missing, drops zero-population geographies, and fits the
PCA on complete cases, as `ReproducibleYostIndex`, `ReproducibleNDI`,
and `ReproducibleSDI` do. Two optional steps:

- **`stabilize = TRUE`** — the CMS shrinkage (`stabilizeCDI()`, the same
  function `computeOriginalCDI()` uses) toward the parent level: block
  group → tract (as CMS), tract → county, or county → state (as in
  `ReproducibleYostIndex`). Only observed values with a margin of error
  are stabilized.
- **`impute = TRUE`** — spatial queen-contiguity imputation
  (`imputeMissing()`), population-weighted by default
  (`weight_var = "tot_pop"`, or `"none"`).

When both are on, **`order`** chooses the sequence:

- **`"stabilize_first"`** (default) — stabilize the observed values,
  then impute the remaining gaps from the stabilized neighbors (the same
  order as `computeYostIndex()`).
- **`"impute_first"`** — impute from the raw neighbors, then stabilize.
  Imputed cells have no margin of error, so they are not stabilized
  themselves, but they count in the tract spread $t^2$ used to stabilize
  their observed neighbors.

The PCA, scoring, and rescaling are as above; `our_CDI` ranks into
percentiles 1–100. The `scope` argument controls the reference
population: `"national"` (default, one PCA) or `"state"` (one PCA per
state).

## Usage

``` r
library(ReproducibleCDI)

# Requires a Census API key:
# tidycensus::census_api_key("YOUR_KEY_HERE", install = TRUE)

# The CMS index, reproduced exactly (2019 block groups, national)
cdi_cms <- computeOriginalCDI(year = 2019, states = "all")
head(cdi_cms$df_cdi)   # ..., CDI_std, cdi_rank, CDI_NATRANK, annotation

# ReproducibleX-style CDI for Delaware block groups (2022 ACS5),
# stabilized toward the tract, then spatially imputed
cdi_de <- computeCDI(
  geo       = "block group",
  year      = 2022,
  states    = "DE",
  stabilize = TRUE,
  impute    = TRUE,
  order     = "stabilize_first"
)

head(cdi_de$df_cdi)   # GEOID, the components, CDI, CDI_std, our_CDI, annotation
```

Key `computeCDI()` options (mirroring `ReproducibleNDI` and
`ReproducibleYostIndex`):

- **`geo`**: `"county"`, `"tract"`, or `"block group"` (`"cbg"`).
- **`states`**: a vector of state abbreviations (e.g. `c("CA", "NY")`),
  or `"all"` for the entire US.
- **`stabilize`**, **`impute`**, **`order`**: see above.
- **`scope`**: `"national"` or `"state"`.
- **`return_format`**: `"detailed"` (default; a list with the data, the
  stabilization diagnostics, and the `psych::principal` object) or
  `"minimal"` (`GEOID`, `CDI_std`, `our_CDI`, `annotation`).

The 18-variable PCA needs more complete geographies than variables; very
small selections (e.g., a handful of counties) are flagged
`"PCA failed"`.

## References

1.  Robst J, Forlines G, Kautter J, et al. The development of the
    Community Deprivation Index and its application to Accountable Care
    Organizations. *Health Affairs Scholar.* 2024.
    [doi:10.1093/haschl/qxae161](https://doi.org/10.1093/haschl/qxae161)
    ([PubMed](https://pubmed.ncbi.nlm.nih.gov/39664484/))
2.  Singh GK. Area deprivation and widening inequalities in US
    mortality, 1969–1998. *American Journal of Public Health.*
    2003;93(7):1137-1143.
    [doi:10.2105/ajph.93.7.1137](https://doi.org/10.2105/ajph.93.7.1137)
3.  Kind AJH, Buckingham WR. Making Neighborhood-Disadvantage Metrics
    Accessible — The Neighborhood Atlas. *New England Journal of
    Medicine.* 2018;378(26):2456-2458.
    [doi:10.1056/NEJMp1802313](https://doi.org/10.1056/NEJMp1802313)

See also the Massachusetts Health Data Tool [Community Deprivation
Index](https://healthdatatool.mass.gov/ss_whatsnewitem/community-deprivation-index/).

## Status

This package is in early development.

## Installation

You can install the development version of `ReproducibleCDI` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("HanLabCollaboration/ReproducibleCDI")
```
