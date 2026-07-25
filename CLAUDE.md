# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

This is an R-based simulation study for **RAILS** (Raking-Assisted Integration of Linked Surveys), a statistical methodology for combining non-probability samples (e.g., "All of Us" / AoU cohort) with probability samples (e.g., NHANES) to produce unbiased population prevalence estimates. The project compares multiple weighting methods across simulation scenarios.

## How to Run

Everything is executed through R Markdown (`.Rmd`) files — there is no build system or test runner.

**Run a simulation:**
```r
# In RStudio or R, knit the relevant Rmd:
rmarkdown::render("Simulations/Simu_methods/S1_simulation.Rmd")
```

**Run a results summary:**
```r
rmarkdown::render("Results/Res_methods/S1_summary.Rmd")
```

Simulations cache results as `.RData` files in `Results/`. If the `.RData` already exists, the simulation block loads from cache rather than re-running.

**Restore packages (Simu_variance only):**
```r
# From within Simulations/Simu_variance/
renv::restore()
```

## Architecture

### Directory Structure

- `Functions/` — shared R source files sourced by all simulation Rmds
- `Simulations/Simu_methods/` — simulation Rmds for scenarios S1–S5 (comparing weighting methods)
- `Simulations/Simu_variance/` — simulation Rmds S3\_1 through S3\_6 (variance estimation comparisons)
- `Results/Res_methods/` and `Results/Res_variance/` — summary Rmds that load `.RData` outputs and produce tables/plots

### Function Files and Their Roles

| File | Purpose |
|------|---------|
| `sample_functions.R` | Generates synthetic population data (`sample.function`): age, sex, income, and two extra covariates; creates AoU (`aou`) and probability sample (`s`) membership via logistic models |
| `supplementary_functions.R` | `expit`, string helpers for calibration formula construction (`transform_term`, `transform_term_aou`), and the core PS estimation solver `fun.vs` (pseudo-likelihood propensity score via Newton-Raphson) |
| `output_functions.R` | `fun.out`: computes weighted mean, variance, CI, bias, nominal coverage, weight diagnostics; optionally wraps `svymean` for survey-design variance |
| `simfun.R` | `sim.fun`: the master simulation function; runs all ~10 competing methods (M1–M10) for one replication and returns a named list |
| `report_function.R` | `fun.rep`: aggregates a list of `sim.fun` outputs across replications into summary tables (bias, empirical SD, average SD, coverage, oracle coverage, variable selection counts) and ggplot-ready long-format data |
| `apply_row_colors.R` | kableExtra helper for coloring table rows in summary Rmds |

### Simulation Methods Compared

`sim.fun` evaluates these methods in order:

1. **unweighted** — naive AoU mean
2. **trueweight** — oracle IPW using known inclusion probabilities
3. **svy\_oneway\_linear / rake** — survey calibration on marginals only
4. **svy\_twoway\_linear / rake** — calibration on marginals + all two-way interactions
5. **svy\_twoway\_rake\_hybrid** — two-way raking with one-way fallback on non-convergence
6. **ps\_oneway\_raw/linear/rake** — pseudo-likelihood PS weighting, one-way model
7. **ps\_twoway\_raw/linear/rake** — pseudo-likelihood PS weighting, two-way model
8. **ps\_vs\_twoway\_\*** — PS with forward likelihood-ratio variable selection
9. **ps\_vs\_stepwise\_twoway\_rake** (G-RAILS) — the proposed method: stepwise VS + raking calibration
10. **oracle\_rails** — RAILS with true model (no VS needed)

Variance correction outputs (`var_correction`, `var_correction_oracle_rails`, `var_correction_svy_oneway_rake`) use stacked estimating equations computed inside `fun.rails.var` and `fun.svy.var` (defined in `supplementary_functions.R`).

### Data Flow

```
sample.function()  →  dt (full population)
                       ├─ dt_aou  (non-probability sample)
                       └─ dt_s    (probability sample)
                              ↓
                        sim.fun() [calls fun.vs, survey::calibrate, fun.out]
                              ↓
                        sim_result list  →  saved as .RData
                              ↓
                        fun.rep()  →  summary tables, plots
```

### Scenario Naming

- `S1`–`S5` vary simulation parameters (`alpha`, `beta`, `gamma`) controlling outcome prevalence, probability-sample selection, and AoU selection bias.
- `S3_1`–`S3_6` (variance series) fix the S3 data-generating process while varying the number of calibration variables or variance estimation approaches.

## Key R Packages

`survey` (calibration, svydesign, svymean), `dplyr`/`tidyverse`, `VGAM`/`distributionsrd` (truncated Pareto for income), `kableExtra`/`DT` (tables), `ggplot2`/`plotly` (plots), `stringr`, `glmnet`.
