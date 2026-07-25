# RAILS Procedure Guide

This folder contains the core estimation functions (`AoU_Fun.R`), the main analysis script (`04_Global_RAILS.R`), and the downstream analyses that consume its weights. It assumes the aggregated datasets and the individual-level AoU dataset have already been prepared — see the [Data Processing guide](../Data%20Processing/README.md) first.

> **Run inside the AoU Researcher Workbench.** Upload `AoU_Fun.R` alongside the analysis scripts on the Workbench disk — they load it via `source("AoU_Fun.R")`.

---

## Required Packages

```r
library(tidyverse); library(survey); library(Matrix)
library(ggplot2); library(ggrepel); library(xtable)
```

---

## Files

| File | Description |
|---|---|
| `AoU_Fun.R` | Helper and estimation functions (documented below) |
| `04_Global_RAILS.R` | Main analysis — loads data, precomputes margins, calls `fun.rails.threeway`, joins weights to participants, saves results |
| `05_Prevalence_Analysis.R` | Flags GBD disease categories from ICD9/10 codes, estimates weighted prevalence under every method, produces national and by-region figures |
| `06_Phecode_Categories.R` | Maps diseases to phecode groupings; prevalence-ratio, relative-bias, and volcano figures colored by phecode category |
| `07_ASCVD_Analysis.R` | ASCVD risk distributions before/after weighting, plus the risk-category shift table |
| `08_NHIS_Comparison.R` | AoU vs NHIS 2020 for care place, self-reported health, and insurance under three weighting schemes; emits the LaTeX table |

---

## Functions in `AoU_Fun.R`

### `fun.out(w)`

Diagnostic summary of a weight vector: sum, variance, proportion non-positive, count below 1, min, max.

---

### `fun.lkd(x, names_univar, LKHD, dt_s, dt_aou, m_s)`

Tests a single candidate three-way interaction `x` for inclusion via a forward likelihood-ratio test.

| Argument | Description |
|---|---|
| `x` | Candidate term, e.g. `"agegroup:sex:income"` |
| `names_univar` | Terms already in the current model |
| `LKHD` | Log-likelihood of the current model |
| `dt_s` / `dt_aou` | Aggregated reference and AoU cell counts (column `weight`) |
| `m_s` | Model matrix for the current model on `dt_s` |

Returns `(log-likelihood, LRT, p-value, LRT/df, df)`, or all `NA` if any cell count is zero or the solver fails (capped at 1000 iterations).

---

### `fun.nps(cat_temp, dt_aou, dt_s, nsiz = sum(dt_s$weight), theta_init = NULL)`

Pseudo-likelihood propensity score weights via Newton-Raphson on the aggregated cells. Returns weights scaled to sum to `nsiz`, with the converged `theta` attached as an attribute for warm-starting the next fit.

---

### `fun.rails.threeway(dt_agg_aou, dt_agg_pums, pop_totals, names_univar, alpha, nsiz)`

Global RAILS estimator over the **two-way base + selected three-way interactions** model space:

1. Fits a two-way NPS model to obtain the baseline log-likelihood.
2. Runs forward LRT selection over all three-way candidates — adds the largest LRT/df among those reaching `alpha`, stops when none qualify.
3. Applies **LIFO stepwise** raking — adds selected terms one at a time in selection order, refitting the NPS model and raking to the population totals at each step; stops at the first model that fails to converge and keeps the weights from the last successful step.
4. Computes benchmark methods on the same cells.

All computation is on the aggregated cross-tabulation using sparse model matrices and warm-started Newton steps. Population totals are **not** computed inside the function: all margins are precomputed once in `04_Global_RAILS.R` and passed in as `pop_totals`; each LIFO step subsets the margins it needs via `create_v3`.

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `dt_agg_aou` | — | Aggregated AoU cell counts |
| `dt_agg_pums` | — | Aggregated reference-sample cell counts |
| `pop_totals` | — | Named vector of population margins (one-way + two-way + three-way) |
| `names_univar` | the seven shared covariates | Main-effect variables; all two- and three-way combinations are considered |
| `alpha` | `0.05` | Significance threshold for forward LRT selection |
| `nsiz` | `sum(dt_agg_pums$weight)` | Target population total the weights are scaled to. Override when the reference sample differs from the source of `pop_totals`, or when targeting a subgroup total |

When a step fails to converge, a **warning** names the covariate combination that failed to pile up to the higher order, and the last successful model's weights are returned.

#### Output

Returns `dt_agg_aou` with appended columns (one row per covariate cell):

| Column | Description |
|---|---|
| `d_unweighted` | Equal weights `nsiz / n` (naive benchmark) |
| `d_cal1` / `d_cal2` | Raking from equal weights to one-way / one-way + two-way margins (no PS) |
| `d_nps1` / `d_nps2` | NPS weights, one-way / two-way model (no raking) |
| `d_nps1_rake` / `d_nps2_rake` | NPS weights as starting point, raked to the same margins |
| `d_rails` | RAILS weights after LIFO stepwise raking |
| `selected_terms` | All interaction terms selected by forward LRT |
| `calibrated_terms` | Terms of the final converged model — shorter than `selected_terms` when the walk stopped early, `NA` if no step converged |

All weight columns are **cell-level totals**; divide by the cell count (`weight`) to obtain per-individual weights. `d_cal2` may be `NA` with a warning when full two-way raking does not converge — itself an informative benchmark result.

---

## Main Analysis (`04_Global_RAILS.R`)

**1. Load data**, **2. apply harmonized factor levels**, then:

**3. Precompute all PUMS margins** once via a single sparse model matrix:

```r
twovars   <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
threevars <- combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))

max_formula <- formula(paste0("~", paste(c(names_univar, twovars, threevars), collapse = "+")))
mat_max     <- sparse.model.matrix(max_formula, data = dt_agg_pums, keep.order = TRUE)
pop_totals  <- as.numeric(Matrix::crossprod(mat_max, dt_agg_pums$weight))
names(pop_totals) <- colnames(mat_max)
```

> Sanity check: `sum(pop_totals == 0)` should be `0`. A zero margin means some covariate combination has no PUMS support, and raking on a term involving it will fail.

**4. Run Global RAILS**

```r
result_rails <- fun.rails.threeway(
  dt_agg_aou   = dt_agg_aou,
  dt_agg_pums  = dt_agg_pums,
  pop_totals   = pop_totals,
  names_univar = names_univar,
  alpha        = 0.05
)
```

**5. Join weights to participants and save** — cell totals divided by the cell count give per-individual weights:

```r
dt_aou_weighted <- dt_aou %>%
  left_join(
    result_rails %>%
      mutate(
        w_unweighted = d_unweighted / weight,  w_cal1      = d_cal1      / weight,
        w_cal2       = d_cal2       / weight,  w_nps1      = d_nps1      / weight,
        w_nps2       = d_nps2       / weight,  w_nps1_rake = d_nps1_rake / weight,
        w_nps2_rake  = d_nps2_rake  / weight,  w_rails     = d_rails     / weight
      ) %>%
      select(all_of(names_univar), starts_with("w_"),
             selected_terms, calibrated_terms),
    by = names_univar
  )
```

Verify: `sum(dt_aou_weighted$w_rails)` ≈ `nsiz`.

---

## Optional: Hybrid Design (NHIS Reference)

The script ends with one flag-gated run: **`run_hybrid <- TRUE`** uses NHIS as the propensity reference sample while raking to the PUMS `pop_totals`. Requires `dt_agg_nhis_v2.csv` (see `03_NHIS_Prep.R`) and an explicit `nsiz`, since the NHIS weight sum does not match the PUMS-based totals scale. Output: `hybrid_rails_weights.csv`.

---

## Relation to the Subgroup Version

`fun.sub.rails.threeway` (in `../../Subgroup RAILS/Sub_AoU_Fun.R`) wraps the same procedure in a loop over the levels of a subgroup variable, fitting a separate model within each stratum.

| | Global RAILS | Subgroup RAILS |
|---|---|---|
| Function | `fun.rails.threeway` | `fun.sub.rails.threeway` |
| Stratification | None | By `subgroup_var` |
| `nsiz` | Total PUMS weight sum | Per-subgroup PUMS weight sum |
| Population totals | Precomputed `pop_totals` vector | Recomputed per subgroup from the PUMS cell subset |

---

## Downstream Analyses

### `05_Prevalence_Analysis.R`

Flags 15 GBD disease categories from ICD9/10 code ranges (a participant is a case with >1 distinct condition date), then estimates weighted prevalence under all eight methods, nationally and by region, against US reference prevalence. Produces the national comparison figures (all methods, US-vs-naive-vs-RAILS, ratio-to-US) and the by-region grids.

### `06_Phecode_Categories.R`

Maps AoU disease categories to phecodes and phecode groups (majority vote across matched ICD codes, with an explicit list resolving ties), then produces the manuscript scatter figures: weighted/unweighted prevalence ratio, relative bias reduction, relative difference in bias, and the volcano plots.

### `07_ASCVD_Analysis.R`

Compares the 10-year ASCVD risk distribution before and after RAILS weighting — histograms, densities, and a risk-category table (low <5%, intermediate 5–10%, high >10%) reporting percentage-point shifts. Figures export as 600 dpi TIFF, vector PDF, and PNG.

### `08_NHIS_Comparison.R`

Compares three dynamic health outcomes — care place, self-reported health, insurance coverage — between NHIS 2020 and AoU under three schemes:

| Column | Weights |
|---|---|
| Unweighted | equal |
| Recalibrate | RAILS |
| Double-Weighting | RAILS ÷ item-response propensity (logistic on the demographic covariates) |

Weights are rescaled to `nsiz` per outcome, after filtering to that outcome's respondents. Confidence intervals are logit-scale, with the variance size taken as each estimate's own sum of weights. Emits the comparison table as LaTeX and CSV.
