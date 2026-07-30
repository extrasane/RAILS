# RAILS

**Raking-Assisted Integration of Linked Surveys** — an R package for combining
a non-probability sample (e.g. a volunteer biobank) with a probability
reference sample to produce population-representative weights.

The core estimator fits a two-way pseudo-likelihood propensity model, selects
three-way interactions by a forward likelihood-ratio test, and applies LIFO
stepwise raking to reference population margins. All computation runs on
aggregated covariate cells. A suite of benchmark weighting methods is returned
alongside the RAILS weights.

## Installation

```r
# install.packages("remotes")
remotes::install_github("extrasane/RAILS", subdir = "RAILS")
```

## Functions

| Function | Purpose |
|---|---|
| `fun.rails.threeway()` | Global RAILS estimator: two-way base + selected three-way interactions, with benchmark methods |
| `fun.sub.rails.threeway()` | Runs the estimator within each level of a subgroup variable |
| `fun.nps()` | Pseudo-likelihood propensity score weights (Newton-Raphson) |
| `fun.lkd()` | Forward likelihood-ratio test for one candidate interaction |
| `fun.out()` | Weight diagnostics (sum, variance, min, max, ...) |

## Inputs

Both `dt_agg_aou` (non-probability) and `dt_agg_pums` (reference) are
**aggregated cell tables**: one row per unique covariate combination, with a
column `weight` giving the cell's summed weight. `pop_totals` is a named vector
of the reference population's one-, two-, and three-way margins, precomputed
once by the caller:

```r
library(Matrix)
vars      <- c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region")
twovars   <- combn(vars, 2, FUN = function(x) paste(x, collapse = ":"))
threevars <- combn(vars, 3, FUN = function(x) paste(x, collapse = ":"))
f         <- formula(paste0("~", paste(c(vars, twovars, threevars), collapse = "+")))
mm        <- sparse.model.matrix(f, data = dt_agg_pums, keep.order = TRUE)
pop_totals <- setNames(as.numeric(crossprod(mm, dt_agg_pums$weight)), colnames(mm))
```

## Usage

```r
library(RAILS)

result <- fun.rails.threeway(
  dt_agg_aou   = dt_agg_aou,
  dt_agg_pums  = dt_agg_pums,
  pop_totals   = pop_totals,
  names_univar = vars,
  alpha        = 0.05
)

# result has one row per covariate cell with per-individual weight columns:
#   d_unweighted, d_cal1, d_cal2, d_nps1, d_nps2, d_nps1_rake, d_nps2_rake, d_rails
# plus selected_terms and calibrated_terms.
```

Subgroup analysis (one model per level of `region`):

```r
result_by_region <- fun.sub.rails.threeway(
  dt_agg_aou   = dt_agg_aou,
  dt_agg_pums  = dt_agg_pums,
  subgroup_var = "region"
)
```

## Notes

- Weight columns are **per-individual** (cell totals divided by the cell
  count). Sum over participants approximates `nsiz`.
- `calibrated_terms` reports the final converged model, which may be shorter
  than `selected_terms` if a raking step failed to converge.
- The AoU Researcher Workbench application scripts that use these functions live
  in the parent repository under `Application/`.

## Development

Documentation (`man/`) is generated from the roxygen2 comments in `R/`:

```r
devtools::document()   # regenerates NAMESPACE and man/
devtools::check()      # R CMD check
```
