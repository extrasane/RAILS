# Data Processing Guide

This folder covers the data preparation steps required before running Global RAILS: building the PUMS reference dataset, extracting the AoU biobank dataset, and (for the hybrid design) preparing NHIS.

> Steps 1 and 3 can be run locally or in any R environment. **Step 2 must be run inside the AoU Researcher Workbench.**

---

## Required Packages

```r
# Step 1 (PUMS)
library(httr); library(jsonlite); library(dplyr); library(tidyverse)

# Step 2 (AoU — Workbench only)
library(tidyverse); library(bigrquery); library(lubridate)

# Step 3 (NHIS)
library(tidyverse); library(survey)
```

---

## Step 1 — PUMS Reference Data (`01_PUMS_Prep.R`)

### 1a. Download from Census API

Downloads 2022 ACS 1-year PUMS microdata for all four Census regions. Results are cached locally as `PUMS_2022.csv` so the API is called only once.

```r
# Request your own free key at: https://api.census.gov/data/key_signup.html
CENSUS_KEY <- "YOUR_CENSUS_API_KEY"
```

`REGION` is requested explicitly in the `get=` list so that region is a per-record variable, guaranteed row-aligned by the API contract.

### 1b. Recode and Harmonize Variables

| Variable | Source field(s) | Categories |
|---|---|---|
| `agegroup` | `AGEP` | 18–24, 25–44, 45–64, 65–74, 75+ |
| `sex` | `SEX` | Female, Male |
| `race_eth` | `RACBLK`, `RACASN`, `RACWHT`, `HISP`, … | Hispanic, NH Asian, NH Black, NH White, Others |
| `income` | `HINCP` | \<35k, 35k–50k, 50k–75k, 75k–100k, \>100k |
| `edu` | `SCHL` | Less than highschool → College graduate or advanced (5 levels) |
| `homeown` | `TEN` | Own, Rent, Others |
| `region` | `REGION` (Census region code 1–4) | Northeast, Midwest, South, West |

Two details matter for correctness:

- **All coded fields are coerced with `as.numeric()` before any comparison.** The Census API returns every field as a character string with zero-padded codes (e.g. `HISP` = `"01"`), so a comparison like `hispanic != 1` would silently misclassify on a fresh download while behaving differently after a cached `read.csv()`.
- **Records are restricted to adults (age > 17) at the top of the recode chain**, and weights are rescaled to the pre-`na.omit` adult total. That anchor is also the `nsiz` used downstream.

> The AoU side maps state of residence to region via a state lookup; PUMS uses the official `REGION` variable. Both produce the same four regions.

### 1c. Aggregate and Save

```r
names_univar <- c("agegroup", "sex", "edu", "homeown", "income", "race_eth", "region")
cat_formula  <- formula(paste0("weight ~", paste(names_univar, collapse = "+")))
dt_agg_pums  <- aggregate(cat_formula, data = dt_pums_analysis, sum)
```

**Output:** `dt_agg_pums_v2.csv` — one row per unique covariate cell, `weight` = sum of PUMS person weights in that cell.

When run inside the Workbench the script uploads to the bucket automatically; when run locally it skips the upload and you copy the file across manually.

---

## Step 2 — AoU Biobank Data (`02_AoU_Prep.R`)

> **Run inside the AoU Researcher Workbench.** Queries controlled-tier BigQuery tables; access requires Workbench registration and dataset approval.

The active code uses the **Workbench 2.0** grammar. Set the CDR before running, either by resource id (from the Resources panel) or directly:

```r
Sys.setenv(AOU_CDR_RESOURCE_ID = "my-controlled-tier-cdr")
# or
Sys.setenv(AOU_CDR_DATASET = "project.dataset")
# optional, to copy outputs to a bucket resource
Sys.setenv(OUTPUT_BUCKET_RESOURCE_ID = "my-output-bucket")
```

Queries run straight into R (`bq_project_query` + `bq_table_download`), with tables referenced fully qualified as `` `project.dataset.table` ``. A Workbench 1.0 implementation (`bq_table_save` → `gsutil` export round-trip) is retained commented out at the top of the file for environments that still require it.

### 2a. BigQuery Queries

Three tables are queried and joined by `person_id`. **All queries restrict the cohort to participants with EHR data** (`has_ehr_data = 1` in `cb_search_person`). If that table is absent from the CDR, the script warns and falls back to all persons — check for the warning, since it changes the cohort and every downstream estimate.

| Query | Table | Fields |
|---|---|---|
| Demographics | `person` | Gender, date of birth, race, ethnicity, sex at birth |
| Survey | `ds_survey` | Income, home ownership, education, usual care place |
| State | `person_ext` | State of residence |

Survey concept IDs: `1585370`, `1585375`, `1585899`, `1585940`, `43530593`.

### 2b. Recode and Harmonize Variables

Variables are recoded to match the PUMS categories exactly. Age is computed as of `AGE_REFERENCE_DATE` (default 2024-08-01; keep fixed so runs are comparable).

- `race`: "None Indicated", "None of these", "Middle Eastern or North African", "Native Hawaiian or Other Pacific Islander", "More than one population" → **Others**; "I prefer not to answer" / "PMI: Skip" → `NA`
- `ethnicity`: Hispanic or Latino → **Yes**; all non-Hispanic responses → **No**; prefer-not-to-answer / skip → `NA`
- `sex`: derived from `sex_at_birth`; Intersex, None, PMI responses → `NA`
- `homeown`: "Other Arrangement" / "PMI: Don't Know" → **Other** (recoded to "Others" in the analysis step)
- `income`, `edu`: mapped from AoU's verbose survey labels to the five-level PUMS-aligned categories

State of residence is extracted from the `PII State: XX` string format. `pivot_wider(values_fn = first)` guards against participants with multiple answers to the same question.

**Outputs:** `aou_raking_dt.csv` (individual level) and `dt_agg_aou_v3.csv` (aggregated cells).

---

## Step 3 — NHIS Reference Data (`03_NHIS_Prep.R`)

> Needed only for the **hybrid design**, where NHIS serves as the propensity reference sample while raking calibrates to PUMS totals.

Downloads the public 2020 NHIS adult file from NCHS, recodes to the same categories, and produces:

- `NHIS_2020.csv` — recoded individual-level data (`weight` = `WTFA_A`)
- `dt_agg_nhis_v2.csv` — weighted cell counts, ready to pass to `fun.rails.threeway` in the reference-sample slot

```r
result_hybrid <- fun.rails.threeway(
  dt_agg_aou  = dt_agg_aou,
  dt_agg_pums = dt_agg_nhis,            # NHIS as propensity reference
  pop_totals  = pop_totals,             # still computed from PUMS
  nsiz        = sum(dt_agg_pums$weight) # PUMS adult total
)
```

---

## Outputs Summary

| File | Created in | Used in |
|---|---|---|
| `dt_agg_pums_v2.csv` | Step 1 | RAILS Procedure |
| `aou_raking_dt.csv` | Step 2 | RAILS Procedure (individual-level join) |
| `dt_agg_aou_v3.csv` | Step 2 | RAILS Procedure |
| `NHIS_2020.csv` | Step 3 | Hybrid-design analyses |
| `dt_agg_nhis_v2.csv` | Step 3 | RAILS Procedure (hybrid design only) |
