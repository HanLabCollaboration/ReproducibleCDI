# Plan: CMS-faithful CDI + ReproducibleX-style CDI

Date: 2026-10-26

> **Status (2026-10-26): implemented.** `computeOriginalCDI(2019, "all")` reproduces
> `CDI_NATRANK.csv` exactly: all 220,333 block groups, identical suppression labels
> (2,649 GQ, 671 GQ-PH, 1,393 PH), and 100 % identical ranks for the 215,620 ranked
> block groups (see `dev/validate_original_cdi_2019.txt`). The CMS national file
> includes Puerto Rico, so `computeOriginalCDI(states = "all")` does too.
> 71 unit tests pass; `R CMD check` has no new issues.
Reference: https://github.com/CMSgov/community-deprivation-index
(`cdi_construction_s1.sas` = variables, MOEs, shrinkage, imputation, suppression;
`cdi_construction_s2.sas` = PCA, scoring, ranking; `CDI_NATRANK.csv` = published 2019 BG ranks)

## Goals

The package will expose **two** entry points that share the same variable definitions
and the same stabilization function:

| | `computeOriginalCDI()` | `computeCDI()` |
|---|---|---|
| Purpose | Reproduce **every** CMS step exactly | CDI built the ReproducibleX way |
| Variables | The 18 CMS variables | The same 18 CMS variables |
| Zero denominator | → 0 (CMS) | → `NA` (as in Yost / NDI / SDI) |
| Stabilization | CMS shrinkage BG → tract | Same function (`stabilizeCDI()`), optional; BG → tract, tract → county, county → state (as Yost) |
| Imputation | Hierarchical BG → tract → county → state | Spatial queen-contiguity (`imputeMissing()`), optional |
| Order | Fixed: stabilize, then impute | User choice: `order = "stabilize_first"` or `"impute_first"` |
| Suppression | PH / GQ / GQ-PH | `tot_pop == 0` only (as siblings) |
| PCA weights | score coefficients | score coefficients (same ranks as loadings) |
| Ranking | SAS `proc rank groups=100` | `dplyr::ntile(…, 100)` (as siblings) |
| Geographies | block group only | county / tract / block group |

Variable rule for both: use exactly the CMS ACS cells. If a CMS table is not published
at block group for a year, that year is tract-only (no substitute tables).

---

## 1. What CMS actually does (summary of the SAS)

### 1.1 Inputs
ACS 5-year data at **four levels**: block group, tract, county, state. Every variable is
computed at all four levels from that level's own ACS table (tract values are *not*
aggregated from BGs). MOEs are only computed at BG.

### 1.2 The 18 PCA variables (+1 computed but unused)

| CMS name | Formula (CMS: denom 0 → 0) | MOE approach |
|---|---|---|
| Educ_12less_perc | sum(B15003_002:016) / B15003_001 | moe3 |
| Educ_ba_pl_perc | sum(B15003_022:025) / B15003_001 | moe3 |
| Emp_white_col_perc | (C24010_003+039+027+063) / C24010_001 | moe3 |
| Fam_bel_pov_perc | B17010_002 / B17010_001 | moe2 |
| Hhld_1pl_room_perc | (B25014_005:007 + 011:013) / B25014_001 | moe3 |
| Hhld_no_int_perc | B28002_013 / B28002_001 | moe2 |
| Hhld_no_veh_perc | (B25044_003+010) / B25044_001 | moe3 |
| Hous_no_plumb_perc | (B25016_007+016) / B25016_001 | moe3 |
| Inc_dis_upd | 100 · sum(B19001_002:004) / sum(B19001_014:017); missing if either sum is 0. **Kept on the ratio scale through shrinkage and imputation, then `log()`** | moe3 on the ratio, then ×100 |
| Med_gross_rent | B25064_001 | moe1 (published MOE) |
| Med_hhld_inc | B19013_001 | moe1 |
| Med_home_val | B25077_001 | moe1 |
| Med_month_mort | B25088_002 | moe1 |
| No_ins_perc | (B27010_017+033+050+066) / B27010_001 | moe3 |
| One_par_hhld_perc | (B11003_010+016) / B11003_001 | moe3 |
| Own_occ_hous_perc | B25003_002 / B25003_001 | moe2 |
| Pop_bel_150_pov_perc | sum(C17002_002:005) / C17002_001 | moe3 |
| Unemp_perc | B23025_005 / B23025_002 | moe2 |
| *Educ_hs_pl_perc* | *sum(B15003_017:025) / B15003_001 — computed, not in PCA* | *moe3* |

MOE formulas (BG only):
- **moe2 / moe3**: `n_moe = sqrt(sum(num_moe_i^2))`, `d_moe = sqrt(sum(den_moe_j^2))`
  (both numerator **and** denominator may have several cells — needed for `Inc_dis_upd`),
  `p = num/den`; `x = n_moe^2 − p^2·d_moe^2`; if `x > 0` use `sqrt(x)/den`, else
  `sqrt(n_moe^2 + p^2·d_moe^2)/den`. Only when `den ≠ 0` (else MOE missing).
  No special zero-estimate rule (all cell MOEs are summed in quadrature), so
  **do not use `tidycensus::moe_sum(estimate=...)`**.
- The "subtract" proportion formula is also applied to `Inc_dis_upd`, though it is a
  ratio. Replicate as-is and document.

### 1.3 Shrinkage (BG toward its tract)
Per variable `v`, for each BG in tract `T`:
- `bg_count` = number of BG rows in tract `T` (**all** rows, incl. zero-pop, suppressed,
  and BGs where `v` is missing).
- `sqdiff_T` = Σ over BGs in `T` with both `v_bg` and `v_T` non-missing of `(v_bg − v_T)^2`.
- `t2 = sqdiff_T / (bg_count − 1)`; `s2 = (moe / 1.645)^2`;
  `w = (1/s2) / (1/s2 + 1/t2)`; `v_shrunk = w·v + (1 − w)·v_T`.
- `v_shrunk = v` (not shrunk) when `v` or `v_T` is missing, MOE is missing,
  `bg_count = 1`, or `s2 ≤ 0` or `t2 ≤ 0`.

### 1.4 Imputation (hierarchical, after shrinkage)
`v_final = first non-missing of (v_shrunk, v_tract, v_county, v_state)`.
For `Inc_dis_upd`: `log()` of the first value that is non-missing **and non-zero**.
Because rates with a zero denominator are 0 (not missing), imputation in practice hits
mainly the four medians and `Inc_dis_upd`. CMS records nothing about which cells were
imputed.

### 1.5 Suppression (BG only)
- `PH`: `B01001_001 < 100` or `B25032_001 < 30`.
- `GQ`: `B09019_<gq> / B09019_001 ≥ 0.3333`.
- `GQ-PH`: both. Suppressed BGs are excluded from standardization, PCA and ranking, and
  their rank is the flag string. They **are** still used in the tract `t2`.

### 1.6 PCA, score, rank
- `proc factor method=principal nfactors=1` (correlation matrix) on non-suppressed BGs,
  listwise complete.
- `proc standard` to mean 0 / SD 1 over non-suppressed BGs (per-variable available cases).
- `raw_score = Σ z_j · scorecoef_j` (score coefficients = loadings / eigenvalue).
- Rescale to mean 100 / SD 20 (reported, not used for rank).
- `proc rank groups=100` on `raw_score` (ties = mean), then +1 → 1–100:
  `floor(rank · 100 / (n + 1)) + 1`. This is **not** `dplyr::ntile()`.
- No explicit sign alignment.

---

## 2. Shared building blocks

### 2.1 One variable spec — `R/cdiVarSpec.R` (new)
Internal tibble, one row per CMS variable:

```
var | type ("rate","median","ratio_log") | num (cells) | den (cells) | moe_type (1/2/3) | in_pca
```

Drives fetching, variable calculation, MOE computation, stabilization and the PCA list
in both entry points. A unit test checks it cell-for-cell against the CMS `&censusvar`
list.

Names: CMS names in lowercase snake case (`educ_12less_perc`, …, `inc_dis_upd`,
`med_hhld_inc`, `one_par_hhld_perc`). The current `_imp` suffixes go away. No other
tidySDOH package references the current names (checked ReproducibleSDOH, -app, -data,
sherlock_tidySDOH).

### 2.2 Availability guard — `R/checkCDIAvailability.R` (new)
After fetching BG data, if any component's cells are all `NA` for the whole request,
stop with: `"<table> is not published at block group for ACS5 <year>; use geo = 'tract'."`
Never substitute another table.

Probe so far (RI, Bristol County, `_001` cells, 2026-10-26):

| Year | 18 component tables at BG? | Note |
|---|---|---|
| 2016 | **No** — B28002 does not exist at any geography | CDI impossible → `year >= 2017` (current `>= 2016` is a bug) |
| 2017–2023 | Yes | — |

Suppression GQ cell (only `computeOriginalCDI()` needs it):

| Year | GQ cell | At BG? |
|---|---|---|
| ≤ 2018 | `B09019_038` (`_026` is a household-type cell!) | yes |
| 2019–2022 | `B09019_026` | yes |
| 2023–2024 | `B09019_026` | **no** (B26001 neither) → `computeOriginalCDI()` limited to 2017–2022 |

### 2.3 Fetching — `R/fetchCDIAcs.R`
- Keep estimates **and** MOEs (`E` / `M` columns) when stabilization is requested.
- Optional suppression cells (`B01001_001`, `B25032_001`, year-specific B09019 cells).
- New `fetchCDIHierarchy(geo, year, states, levels)` fetching the target level plus
  whichever parents are needed (tract for stabilization; tract + county + state for
  hierarchical imputation).
- Estimates and MOEs are used exactly as tidycensus returns them: jam codes already come
  back as `NA` (checked 2019 BG, CA/NY/MA/TX: no negative values), and top/bottom-coded
  medians come back as numeric caps (e.g. 250,001, 2,000,001, 3,501, 4,001) with `NA`
  MOE, so they are never stabilized. No extra handling (decision D3).

### 2.4 Variable calculation — `R/calculateCDIVars.R`
`calculateCDIVars(data, spec, zero_denominator = c("zero", "na"), with_moe)`:
- `"zero"` → CMS rule (used by `computeOriginalCDI()`).
- `"na"` → sibling-package rule: plain division, `0/0 = NaN → NA` (used by `computeCDI()`).
  Confirmed in Yost (`calculateYostVars`), SDI (`calculateSDIVars`), NDI Powell-Wiley/ICE.
- Always returns `inc_dis_upd` on the **ratio** scale; `log()` is applied later by
  `finalizeIncDis()`.
- With `with_moe`, adds `<var>_moe` via `R/cdiMoe.R` (new): ports of `moe1/moe2/moe3`,
  vectorised (the Yost `moe3` is `rowwise()` and supports only one denominator cell;
  CDI needs four for `Inc_dis_upd`).

### 2.5 Stabilization — `R/stabilizeCDI.R` (new, used by both entry points)
`stabilizeCDI(data, data_parent, varlist, geo_parent_key_len)`: a literal port of
§1.3, same interface and output columns as `ReproducibleYostIndex::stabilizeYostIndex()`
(`<var>_parent`, `<var>_sqdiff`, `<var>_s2`, `<var>_wgt`, `<var>_stabilized`,
`stabilizedflg_<var>`), plus `nvar_stabilized` per row.
- **Row count (`bg_count`) = `dplyr::n()` within the parent**, exactly as Yost: every row
  passed in counts, including rows where the variable is missing. The function has no
  option for this; the two entry points differ only in which rows they pass:
  - `computeOriginalCDI()` passes **all** BGs (zero-pop and suppressed included) → CMS.
  - `computeCDI()` drops `tot_pop == 0` first, then stabilizes → Yost
    (`processYostScope()` does the same).
- Parent level follows the Yost hierarchy: block group → tract (`key_len = 11`),
  tract → county (`5`), county → state (`2`). `computeOriginalCDI()` uses only
  BG → tract (CMS); tract and county stabilization in `computeCDI()` are extensions
  beyond CMS, kept for consistency with Yost.
- Cells without an MOE (missing values, imputed values) are never stabilized themselves.
- `stabilizedflg_<var>` = 1 only when shrinkage was actually applied (value, parent value
  and MOE present, `count > 1`, `s2 > 0`, `t2 > 0`), else 0 — never `NA`.
  ReproducibleYostIndex is not changed by this plan.

### 2.6 Imputation
- `imputeHierarchical(data, parents, varlist)` (new): literal port of §1.4. Used only by
  `computeOriginalCDI()`.
- `imputeMissing()` (spatial queen): unchanged, identical to the sibling packages. Used
  by `computeCDI()`.

### 2.7 PCA and score — `R/fitCDIPCA.R` (new, factored out of `processCDIScope()`)
Standardize, `psych::principal(nfactors = 1, rotate = "none")`, composite
`z %*% fit$weights`, sign-aligned so `fam_bel_pov_perc` loads positive, rescale to
mean 100 / SD 20. Returns score, coefficients, fit object. The ranking function is
chosen by the caller.

---

## 3. `computeOriginalCDI()` — exact CMS replication (new)

```r
computeOriginalCDI(year, states = "all", quiet = FALSE, return_format = "detailed")
```

Steps, in CMS order, with no options to change them:
1. Fetch BG (estimates + MOEs + suppression cells), tract, county, state.
2. `calculateCDIVars(zero_denominator = "zero", with_moe = TRUE)` at BG; same without
   MOE at tract / county / state. Includes `educ_hs_pl_perc` (not used in PCA).
3. Suppression flags (§1.5) using the year-specific GQ cell.
4. `stabilizeCDI()` over **all** BG rows (suppressed included).
5. `imputeHierarchical()` (BG → tract → county → state), then `log()` on `inc_dis_upd`
   with the non-zero rule.
6. Drop suppressed BGs; standardize; PCA; score; rescale 100 / 20.
7. `rankGroupsSAS()` (new `R/rankGroupsSAS.R`) → `CDI_NATRANK` 1–100; suppressed BGs get
   `"GQ"`, `"PH"` or `"GQ-PH"`.

Output: `minimal` = `GEOID`, `CDI_NATRANK` (character, as CMS) and numeric `cdi_rank`;
`detailed` adds raw / `_moe` / `_stabilized` / final values per variable,
`suppression_flag`, the imputation source level per variable (`<var>_imp_src` ∈
BG/tract/county/state — extra information, does not change results), stabilization
diagnostics, PCA coefficients.

Partial `states`: allowed, with a message that ranks are relative to the requested states
(CMS ranks nationally).

Years: `year` must be in 2017–2022 (decision D2). The check runs **first**, before any
Census API call, and uses `stop()` (an error, not a warning — nothing is computed):
- `year > 2022`:
  `"computeOriginalCDI() supports ACS5 2017-2022 only. For <year>, the Census Bureau does not publish table B09019 (group quarters) at the block-group level, so the CMS GQ suppression step cannot be reproduced. Use computeCDI() for <year>."`
- `year < 2017`:
  `"computeOriginalCDI() supports ACS5 2017-2022 only. Table B28002 (household internet), one of the 18 CDI variables, is not published before ACS5 2017."`
- Unit tests assert both errors (e.g. `expect_error(computeOriginalCDI(2023), "2017-2022")`)
  without hitting the API.

GQ cell: `B09019_038` for 2017–2018, `B09019_026` for 2019–2022.

---

## 4. `computeCDI()` — ReproducibleX style with user-selected order

```r
computeCDI(
  geo, year, states,
  stabilize = FALSE,                                 # CMS shrinkage via stabilizeCDI(); toward parent level (as Yost)
  impute    = FALSE,                                 # spatial queen via imputeMissing()
  order     = c("stabilize_first", "impute_first"),  # used only when stabilize & impute
  quiet = FALSE, return_format = "detailed", scope = "national",
  weight_var = "tot_pop"
)
```

Common steps (as today and as the siblings):
1. Fetch (+ MOEs and the parent level if `stabilize`; geometry if `impute`).
2. `calculateCDIVars(zero_denominator = "na", with_moe = stabilize)`.
3. Drop `tot_pop == 0` (returned later as "No population").
4. Stabilize and/or impute according to `order` (below).
5. `log()` on `inc_dis_upd` (ratio ≤ 0 → `NA`).
6. Complete cases → PCA → score → 100 / 20 → `ntile(…, 100)` → annotation.

### 4.1 `order = "stabilize_first"` (default)
1. `stabilizeCDI()` on the observed data only. Missing cells stay missing and contribute
   nothing to `sqdiff`; observed cells with an MOE are pulled toward the tract value.
2. `imputeMissing()` (queen) fills the remaining gaps using the **stabilized** neighbour
   values.

This is the same order as `computeYostIndex(stabilize = TRUE, impute = TRUE)`.

### 4.2 `order = "impute_first"`
1. `imputeMissing()` (queen) fills gaps from the **raw** neighbour values.
2. `stabilizeCDI()` on the imputed data. Imputed cells have no MOE, so they are not
   stabilized themselves, but they are now non-missing and enter `sqdiff_T`, which
   changes `t2` and therefore the weights of their observed tract siblings.

Effects to document: the two orders differ in (a) whether neighbours are averaged
before or after stabilization and (b) whether imputed values feed the tract variance
`t2`. (With CMS's hierarchical imputation the two orders would be identical, because a
BG filled with its tract value adds 0 to `sqdiff`. That is why the order choice lives
in `computeCDI()` only.)

### 4.3 `inc_dis_upd` scale
Stabilize **and** impute on the ratio scale (the MOE is on that scale), then `log()`,
mirroring CMS's "log last" rule.

### 4.4 Outputs
As today, plus when `stabilize = TRUE`: `<var>_raw`, `<var>_raw_moe`,
`stabilizedflg_<var>`, `nvar_stabilized` (Yost naming convention), and
`settings = list(stabilize, impute, order)`.
`annotation` adds "Stabilized" / "Stabilized + imputed" levels.
Stabilization and imputation both work at all three levels (§2.5).

---

## 5. File-by-file change list

| File | Change |
|---|---|
| `R/cdiVarSpec.R` | **new** — CMS variable spec |
| `R/checkCDIAvailability.R` | **new** — BG availability guard |
| `R/fetchCDIAcs.R` | spec-driven; MOEs; suppression cells; `fetchCDIHierarchy()` |
| `R/cleanAcsNamesCDI.R` | handle `state` level |
| `R/calculateCDIVars.R` | spec-driven; `zero_denominator`; ratio-scale `inc_dis_upd`; `educ_hs_pl_perc`; `_moe`; CMS names |
| `R/cdiMoe.R` | **new** — vectorised `moe1/2/3` |
| `R/stabilizeCDI.R` | **new** — shared stabilization |
| `R/imputeHierarchical.R` | **new** — CMS imputation |
| `R/imputeMissing.R` | unchanged |
| `R/flagCDISuppression.R` | **new** — PH / GQ |
| `R/rankGroupsSAS.R` | **new** |
| `R/fitCDIPCA.R` | **new** — PCA + score, shared |
| `R/computeOriginalCDI.R` | **new** — exported |
| `R/computeCDI.R` | `stabilize`, `order`, `year >= 2017`, docs |
| `R/processCDIScope.R` | orchestrates §4 (stabilize / impute / order) |
| `R/ReproducibleCDI-package.R` | `globalVariables` |
| `tests/testthat/` | **new** |
| `DESCRIPTION` | fix description (it claims hierarchical imputation today); `testthat` in Suggests |
| `README.Rmd` | two entry points, order option |
| `dev/compare_to_cms.qmd` | validation (§6) |

---

## 6. Implementation order and verification

1. **Availability probe** (`dev/`): every component cell and the suppression cells at BG,
   ACS5 2017–2023, 3–4 states. Freeze into `checkCDIAvailability()`.
2. **Spec, variables, MOE.** Tests: spec equals CMS `&censusvar`; `zero_denominator`
   both ways; MOE ports vs hand-computed cases (`x_sub ≤ 0` branch, zero denominators,
   multi-cell denominator, ratio ×100).
3. **`stabilizeCDI()`.** Tests on a toy tract: `bg_count` counts rows with missing values;
   `t2`, `s2`, `w`, stabilized value; every "not stabilized" branch; output identical to
   `stabilizeYostIndex()` on shared inputs.
4. **`imputeHierarchical()` + `rankGroupsSAS()` + suppression.** Tests: fallback chain;
   `inc_dis_upd` non-zero rule; SAS rank formula on small vectors with ties; PH/GQ
   thresholds and year-specific GQ cell.
5. **`computeOriginalCDI(2019)` vs `CDI_NATRANK.csv`.** Targets: suppression flags ≈ 100 %
   agreement; rank exact-match rate as high as possible. Investigate residuals (e.g. sign,
   data vintage) and document them.
6. **`computeCDI()` orders.** Tests: `stabilize_first` with no missing data equals
   `impute_first` (nothing to impute); with missing data they differ only through the
   two channels in §4.2. In `dev/compare_to_cms.qmd`, compare to CMS:
   no stabilize/no impute, stabilize only, impute only, `stabilize_first`, `impute_first`.
7. Docs, README, `roxygen2::roxygenise()`, `R CMD check`.

---

## 7. Decisions (all settled 2026-10-26)

- ~~D1 – `bg_count` in `computeCDI()`~~ **Decided 2026-10-26: match Yost** (§2.5).
  `n()` over the rows passed in; `computeOriginalCDI()` passes all BGs, `computeCDI()`
  drops `tot_pop == 0` first. `computeCDI()` also stabilizes tract → county and
  county → state, like Yost.
- ~~D2 – GQ suppression for ACS5 ≥ 2023~~ **Decided 2026-10-26: `computeOriginalCDI()`
  only accepts ACS5 2017–2022** and stops otherwise. B09019 is 100 % `NA` at block group in
  2023 and 2024 (confirmed RI Bristol and CA Los Angeles; still published at tract), so the
  GQ step cannot run. `computeCDI()` is unaffected (no CMS suppression, never requests
  B09019) and accepts `year >= 2017`.
- ~~D3 – Top-coded medians~~ **Decided 2026-10-26: rely on tidycensus.** Use values as
  tidycensus returns them (caps stay numeric, jam codes are already `NA`); no switch,
  no special handling (§2.3).
