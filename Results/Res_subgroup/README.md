# Res_subgroup — Subgroup Calibration Reporting

Reporting side of the subgroup calibration simulation. It consumes the
per-replicate result objects produced by the run side in
[`Simulations/Simu_subgroup`](../../Simulations/Simu_subgroup) and renders the
paper's tables and figures. Mirrors the `Res_methods` / `Res_variance` layout.

## Files

| File | Role |
|---|---|
| `subgroup_summary.Rmd` | Per-setting report. Defines `fun.rep` (the aggregator), loads one setting's results, and renders population + subgroup tables, box/violin plots, the coverage heatmap, timing, average size, and selected-variable histograms. |
| `Heat.R` | Cross-setting subgroup-coverage heatmap: loads (or builds) every setting's `fun.rep` summary and tiles their subgroup coverage side by side. |

## Inputs (per setting, named by the 3-bit code `<code>`)

Written by the run side alongside each `X5_V12_<code>.RData`:

- `X5_V12_<code>.RData` — the length-`n` list of `sim.fun` outputs (`sim_result`).
- `pre_pop_X5_V12_<code>.csv` — true prevalence: row 1 is the population truth
  `truey`; the rest are the subgroup truths `pre_pop`.
- `pre_all_X5_V12_<code>.csv` — the calibration-margin truths (`pre_all`).

The 3-bit code is `(covariate dependence, outcome heterogeneity, sampling
heterogeneity)`; see the [run-side README](../../Simulations/Simu_subgroup/README.md)
for the factor definitions and the eight-method list.

## Running

```r
# One setting: set `case` (e.g. "000") near the top of the Rmd, then knit.
rmarkdown::render("subgroup_summary.Rmd")
```

`fun.rep` output is cached to `X5_V12_summary_<code>.RData`, so re-knitting a
setting reloads rather than recomputes. `var_method = c("rails","rails_subgroup")`
names the G-RAILS and S-RAILS objects; `var_stratified = "catx5"` is the
separator `z`.

`Heat.R` expects `fun.rep`, `sim_result`, `n`, `pre_pop`, `truey`, `var_method`,
and `var_stratified` to be in scope (i.e. run it after a `subgroup_summary.Rmd`
session); it reuses the cached `*_summary_*.RData` where present.

## Reported quantities

- **Population** — relative bias (mean/median), average / empirical / median SD,
  MAD, nominal coverage, oracle coverage, NA ratio.
- **Subgroup** (per `z` level) — bias, |bias|, bias², coverage, relative bias,
  and a size/mean/variance/coverage summary table.
- **Diagnostics** — elapsed time, average subgroup size, average number of
  variables selected, and RAILS / S-RAILS selected-term rankings and histograms.

Outputs (`.RData`, `.html`, `.csv`) are **git-ignored**; only source is tracked.
