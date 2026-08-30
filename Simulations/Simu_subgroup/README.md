# Simu_subgroup — Subgroup Calibration Simulation

Simulation study for the **subgroup** RAILS paper (*Subgroupwise Calibration for
Large-Scale National Biobanks*). It compares **global** and **subgroup**
calibration strategies for integrating a non-probability sample with a
probability reference sample, when the population is partitioned by a
categorical separator `z`.

This is a self-contained companion to the main-paper simulations in
[`../Simu_methods`](../Simu_methods) and [`../Simu_variance`](../Simu_variance).
Unlike those, every function this study needs is bundled in one source file
(`subgroup_functions.R`) rather than sourced from the repo-level `Functions/`.

---

## Files

This folder (`Simulations/Simu_subgroup`) holds the **run** side; the **reporting**
side lives in [`Results/Res_subgroup`](../../Results/Res_subgroup), mirroring the
`Simu_methods` / `Res_methods` split.

| File | Location | Role |
|---|---|---|
| `subgroup_functions.R` | `Simu_subgroup` | All shared code: data generation, the estimators, and the diagnostic helpers (see [Functions](#related-functions) below). |
| `subgroup_simfun.R` | `Simu_subgroup` | The master driver `sim.fun`: runs all eight methods for **one** replicate and returns the list `fun.rep` consumes. |
| `subgroup_simulation.Rmd` | `Simu_subgroup` | Per-setting **run** driver. Set `case` (the three-bit code), then knit to generate `X5_V12_<case>.RData` and the truth CSVs. Sources the two files above. |
| `subgroup_summary.Rmd` | `Res_subgroup` | Loads one setting's `sim_result`, aggregates it with `fun.rep`, and renders the result tables, box/violin plots, coverage heatmap, and selected-variable histograms. |
| `Heat.R` | `Res_subgroup` | Builds the **cross-setting** subgroup-coverage heatmap (all eight settings side by side). |

`subgroup_simulation.Rmd` assembles each setting from the three-bit `case` code
`(covariate dependence, outcome heterogeneity, sampling heterogeneity)` — the
same code that suffixes the result file `X5_V12_<case>.RData` — by toggling
the covariate distributions (`p1`–`p4`), the outcome coefficients `alpha`, and
the non-probability coefficients `gamma` across the four `z` subgroups. The
homogeneous/heterogeneous parameter blocks are hard-coded once and shared across
settings.

---

## Data-generating process

`sample.function()` and `outcome.function()` in `subgroup_functions.R`.

A finite population of size `nsiz` is generated with a categorical **separator**
`z ∈ {1,2,3,4}` (`catx5`) drawn from `z_prob`, and four auxiliary covariates
generated **conditionally on `z`**:

- `x1, x2` — Bernoulli, level-specific probabilities (`p1_list`, `p2_list`);
- `x3` — 4-category multinomial (`p3_list`);
- `x4` — Bernoulli whose probability depends on the `z` level (`p4_list`, a
  base + slope), inducing covariate–separator dependence.

`outcome.function()` then draws three Bernoulli indicators from logistic models
with up to two-way interactions (`x1x2`, `x2x4`), fit **per `z` stratum**:

| Indicator | Coefficients | Meaning |
|---|---|---|
| `y`   (outcome)                 | `alpha` (list by `z`) | the prevalence estimand is `mean(y)` |
| `s`   (probability sample)      | `beta`  (shared)      | reference-sample membership |
| `aou` (non-probability sample)  | `gamma` (list by `z`) | biobank membership |

Survey weights for the probability sample are the inverse of its true inclusion
probability. The estimand is the population prevalence `mean(y)` overall and
within each `z` subgroup.

## Experimental settings (the 3-bit code)

Each setting toggles three factors on/off; the eight combinations are named by a
**three-bit code** `abc` (used in the `X5_V12_<code>` file names):

| Bit | Factor | On = 1 |
|---|---|---|
| 1st | **Covariate dependence** — covariate parameters differ across `z` levels | `p*_list` varies by level |
| 2nd | **Outcome heterogeneity** — the `y` model (`alpha`) differs across `z` | `alpha` varies by level |
| 3rd | **Sampling heterogeneity** — the non-prob model (`gamma`) differs across `z` | `gamma` varies by level |

So `000` is fully homogeneous and `111` toggles all three. The manuscript reports
six of the eight (`000, 100, 010, 001, 101, 111` → paper Settings 1–6); all eight
are produced here. Sampling heterogeneity (3rd bit) is the dominant factor that
decides whether subgroup calibration beats global calibration.

## Methods compared

Eight estimators, in the fixed order `fun.rep` reports them:

| # | Method | Strategy | Implemented by |
|---|---|---|---|
| 1 | **Naive**    | none (unweighted biobank mean) | — |
| 2 | **Oracle**   | true `paou` inclusion probabilities | — |
| 3 | **G-Raking** | global raking to one-way population margins | `fun.oneway.rake` |
| 4 | **NPS**      | global nested propensity score, main effects only | `fun.nps` |
| 5 | **GVS**      | global NPS + greedy forward interaction selection | `fun.gvs` (`d_vs`) |
| 6 | **G-RAILS**  | GVS + LIFO stepwise raking (global) | `fun.gvs` → `fun.lifo` |
| 7 | **S-Raking** | raking to marginal totals **within each `z` subgroup** | `fun.subgroup.rake` |
| 8 | **S-RAILS**  | full RAILS **within each `z` subgroup**; falls back to the G-RAILS weights when a subgroup fails to converge | per-subgroup `fun.gvs` → `fun.lifo` (assembled in `sim.fun`) |

## <a name="related-functions"></a>Related functions (`subgroup_functions.R`)

Grouped as in the file:

- **Utilities** — `expit`, `.log1pexp`, `.solve_ridge_mat`.
- **Margin string helpers** — `transform_term`, `transform_term_aou`,
  `transform_term_sub` (and `.transform_single` / `.transform_interaction`):
  turn model-matrix column names into the `sum(catxV == v & …)` count expressions
  that build population totals for `survey::calibrate`. The `_sub` variant adds a
  subgroup indicator so totals are computed within a `z` level.
- **Propensity score** — `fun.nps`: pseudo-likelihood PS by Newton–Raphson,
  returning normalized design weights `d` (scaled to `nsiz`) and the
  log-likelihood `LKHD`.
- **Selection** — `fun.lkd` (one LRT step for a candidate interaction) and
  `fun.gvs` (forward greedy selection: main-effects start, adds the most
  informative significant two-way term until none remain).
- **Calibration** — `fun.oneway.rake` (global one-way raking) and `fun.lifo`
  (LIFO: drop the last-selected term until raking converges — the RAILS
  safeguard). `fun.subgroup.rake` runs marginal raking independently per `z`
  level and flags non-convergent subgroups (`index_nonconv`).
- **Data generation** — `sample.function`, `outcome.function` (above).
- **Diagnostics** — `fun.out` (weighted mean, WR-variance CI, bias, coverage,
  weight summaries; optional `svymean` design variance), `svymean.out`,
  `fun.pre0`/`fun.pre` (per-subgroup size/mean/variance/CI/coverage), and
  `fun.ee` (estimating-equation residuals per margin).

## <a name="pipeline"></a>Pipeline

```
sample.function() + outcome.function()      # one population per replicate
        │
        ▼
   sim.fun(seed, …)   ── runs all 8 methods ─▶  per-replicate list:
        │                   $<method>$pop  (mean, var, CI, …)  via fun.out
        │                   $<method>$sub  (size/mean/var/LB/UB/coverage × subgroup) via fun.pre0
        │                   $<method>$time
        │                   $rails$model / $rails$n_stepwise      (G-RAILS selection)
        │                   $rails_subgroup$model_ori/model_lifo  (S-RAILS selection)
        │                   $avg_size, $error
        ▼
   sim_result <- lapply(1:n, sim.fun)  ─saved as─▶  X5_V12_<code>.RData     (n = 1000)
        │
        ▼
   fun.rep(sim_result, n, pre_pop, truey, …)   # in subgroup_summary.Rmd
        │   → population: relative bias, Avg/Emp/Median SD, MAD, nominal &
        │     oracle coverage, NA ratio
        │   → subgroup: bias, |bias|, bias², coverage, relative bias
        │   → selected-variable rankings/histograms, timing, average size
        ▼
   subgroup_summary.Rmd  → per-setting HTML report (cached as X5_V12_summary_<code>.RData)
   Heat.R                → coverage heatmap across all 8 settings
```

`pre_pop` (subgroup true prevalences) and `truey` (population true prevalence)
are read from the `pre_pop_*` / `pre_all_*` CSVs written alongside each result.
`var_method = c("rails","rails_subgroup")` names the G-RAILS and S-RAILS objects;
`var_stratified = "catx5"` is the separator `z`.

## Running

Set `case` at the top of `subgroup_simulation.Rmd` and knit it (once per 3-bit
code); it sources `subgroup_functions.R` + `subgroup_simfun.R`, builds the
truths, runs `n = 1000` replicates, and writes `X5_V12_<case>.RData` into
`Results/Res_subgroup`. The core loop is simply:

```r
sim_result <- lapply(seq_len(n), function(seed)
  suppressWarnings(
    sim.fun(seed, nsiz, alpha, beta, gamma, pars, names_var,
            pre_pop, true_y, tolerance = tolerance)))
```

Then, on the reporting side (`Results/Res_subgroup`): set the matching `case` in
`subgroup_summary.Rmd` and knit for the per-setting tables/figures, and run
`Heat.R` for the cross-setting coverage heatmap.

The `.RData`, `.csv`, and `.html` outputs are intentionally **git-ignored**
(see the repo-root `.gitignore`); only the source (`.R`, `.Rmd`, `.md`) is tracked.

## Provenance (ACCRE / HPC)

The published `V12` results were produced on the ACCRE cluster. The repo files
correspond to the HPC ones as: `subgroup_functions.R` = `sup.fun6.R`,
`subgroup_simfun.R` = `sim.fun3.R`. The HPC driver parallelizes the replicate
loop with `Rmpi` / `doMPI` (`foreach(seed = 1:1e3) %dopar% sim.fun(...)`); the
repo driver is the serial `lapply` equivalent (same results, no cluster needed).
The canonical `sim.fun` signature takes `pre_pop`/`true_y` (no `pre_all` — it was
an unused argument and has been dropped here); `z_prob = c(0.15, 0.35, 0.49, 0.01)`.
