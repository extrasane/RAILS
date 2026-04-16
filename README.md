# README

## Overview

This repository contains the code used to produce all simulation results in:

> **"Synthetic Sampling Weights for Volunteer-Based National Biobanks: A Case Study with the All of Us Research Program"**  


The repository is organized into three main folders: **Functions**, **Simulations**, and **Results**.

---

## Repository Structure

```
├── Functions/
│   ├── simfun.R
│   ├── sample_functions.R
│   ├── output_functions.R
│   ├── report_function.R
│   ├── supplementary_functions.R
│   └── apply_row_colors.R
│
├── Simulations/
│   ├── Sim_methods/
│   │   ├── S1_simulation.Rmd
│   │   ├── S2_simulation.Rmd
│   │   ├── S3_simulation.Rmd
│   │   ├── S4_simulation.Rmd
│   │   └── S5_simulation.Rmd
│   └── Simu_variance/
│       └── S3_simulation_1.Rmd  … S3_simulation_6.Rmd
│
└── Results/
    ├── Res_methods/
    │   ├── S1.RData, S1_summary.RData
    │   ├── S2.RData, S2_summary.RData
    │   ├── S3.RData, S3_summary.RData
    │   ├── S4.RData, S4_summary.RData
    │   ├── S5.RData, S5_summary.RData
    │   ├── S1_summary.Rmd
    │   ├── S2_summary.Rmd
    │   ├── S3_summary.Rmd
    │   ├── S4_summary.Rmd
    │   └── S5_summary.Rmd
    └── Res_variance/
        ├── S3_1.RData … S3_6.RData
        ├── S3_summary_1.RData … S3_summary_6.RData
        ├── S3_full_summary_1.RData … S3_full_summary_6.RData
        ├── S3_summary_1.Rmd … S3_summary_6.Rmd
        └── var_summary.Rmd
```

---

## Functions

All helper functions are stored in the `Functions/` folder and sourced at the top of each simulation file. **No function files need to be run independently.**

| File | Description |
|---|---|
| `simfun.R` | Core simulation function `sim.fun()`. Generates synthetic population, draws biobank and reference samples, fits propensity score and calibration weights (RAILS and benchmark methods), and returns estimates and standard errors for each replicate. |
| `sample_functions.R` | Helper functions for drawing stratified biobank and probability reference samples from the synthetic population. |
| `output_functions.R` | Functions to organize raw simulation replicates into summary metrics (bias, variance, coverage, elapsed time, etc.) across methods. |
| `report_function.R` | `fun.rep()` — post-processes a list of simulation results into formatted output tables for display in the summary R Markdown files. |
| `supplementary_functions.R` | Auxiliary statistical functions, including the RAILS variance estimation function `fun.rails.var()` and the likelihood ratio test helper used for variable selection. |
| `apply_row_colors.R` | Formatting utility for `kableExtra` tables; applies conditional row coloring to summary output tables. |

---

## Simulations

### Simulation Scenarios

Each scenario is defined by a unique combination of population parameters (`alpha`, `beta`, `gamma`) controlling selection into the biobank, participation in the reference survey, and the outcome model, respectively. All scenarios use a synthetic population of size N = 3,340,000, a reference survey of n = 5,500, and a biobank of n = 35,000.

| Simulation | Description | Methods Compared | Results Saved In |
|---|---|---|---|
| **S1** | Baseline scenario. True propensity score model uses main effects only (`catx1`, `catx2`, `catx3`). No extra covariates. | All weighting methods including RAILS, raking, PS-based, and true-weight benchmarks. | `Results/Res_methods/S1.RData` |
| **S2** | Includes additional covariates (`catx4`, `catx5`) relevant to the outcome but not biobank selection. True model is additive. | Same as S1. | `Results/Res_methods/S2.RData` |
| **S3** (×6 sub-scenarios) | Variance estimation study. Reference survey sample size varies across six levels. True model includes interaction terms (`catx1:catx2`, `catx2:catx3`). | Focuses on RAILS variance estimation and coverage properties across sample sizes. | `Results/Res_variance/S3_1.RData` … `S3_6.RData` |
| **S4** | True model includes both extra covariates and interaction terms. Both biobank selection and outcome depend on `catx4`, `catx5`. | Same as S1. | `Results/Res_methods/S4.RData` |
| **S5** | Same covariate structure as S4 but with a different biobank selection mechanism (non-zero `alpha` for extra covariates). | Same as S1. | `Results/Res_methods/S5.RData` |

### How to Run

Each `*_simulation.Rmd` file is self-contained. To reproduce results:

1. Source files will load automatically via `source("../../Functions/...")` calls at the top of each file. Both `Sim_methods/` and `Simu_variance/` sit two levels below the project root, so the relative paths to `Functions/` are consistent across all files.
2. Open the desired `*_simulation.Rmd` in RStudio and click **Knit**, or run `rmarkdown::render()`.
3. If the corresponding `.RData` result file already exists in the `Results/` folder, the simulation will be skipped and the existing file will be loaded. To re-run from scratch, delete or rename the `.RData` file.
4. Each simulation runs **1,000 replicates** via `lapply(1:n, function(seed) sim.fun(...))`. Running a single analysis (Section 4 of the paper) takes approximately **3 minutes** on a standard desktop.

> **Note:** All of Us data is not publicly available. To run the application code, access must be obtained through the All of Us Researcher Workbench. NHIS and PUMS 2022 data are publicly available from the CDC and U.S. Census Bureau, respectively.

---

## Summary Files

After running simulations, the corresponding `*_summary.Rmd` files aggregate results and produce tables and figures.

| Summary File | Loads From | Produces |
|---|---|---|
| `S1_summary.Rmd` | `Results/Res_methods/S1.RData` | Population prevalence table, bias boxplot, variance correction coverage table (Table/Figure in main text or supplement). |
| `S2_summary.Rmd` | `Results/Res_methods/S2.RData` | Same outputs as S1 for Scenario 2. |
| `S3_summary_1.Rmd` … `S3_summary_6.Rmd` | `Results/Res_variance/S3_1.RData` … `S3_6.RData` | Per-sub-scenario variance and coverage summary tables. |
| `S3_summary.Rmd` | All six `S3_summary_*.RData` files | Combined summary across all S3 sub-scenarios. |
| `S4_summary.Rmd` | `Results/Res_methods/S4.RData` | Same outputs as S1 for Scenario 4. |
| `S5_summary.Rmd` | `Results/Res_methods/S5.RData` | Same outputs as S1 for Scenario 5. |
| `var_summary.Rmd` | `Results/Res_variance/S3_full_summary_1.RData` … `S3_full_summary_6.RData` | Coverage plots (Nominal vs. Oracle) across RAILS, Cal_1, True_RAILS methods; confidence interval figures by SD type (Empirical, Simplified, Stacked). Saves `Nominal-Oracle Coverage.png`, `CI Emp SD.png`, `CI Avg SD.png`. |

---

## R Package Dependencies

All simulation and summary files require the following R packages. Please ensure these are installed before running any file.

| Package | Version Used |
|---|---|
| `knitr` | 1.49 |
| `dplyr` | 1.1.4 |
| `tidyverse` | 2.0.0 |
| `cat` | 0.0-7 |
| `survey` | 4.2-1 |
| `kableExtra` | 1.4.0 |
| `VGAM` | 1.1-8 |
| `distributionsrd` | 0.0.6 |
| `ggplot2` | 3.5.0 |
| `truncnorm` | 1.0-9 |
| `reshape2` | 1.4.4 |
| `paletteer` | 1.6.0 |
| `plotly` | 4.10.4 |
| `stringr` | 1.5.0 |
| `purrr` | 1.0.1 |
| `RVAideMemoire` | 0.9-83-7 |
| `glmnet` | 4.1-7 |
| `DT` | 0.33 |
| `cowplot` | 1.1.1 |

>R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

---

