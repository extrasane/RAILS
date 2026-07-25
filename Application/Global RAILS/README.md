# Global RAILS

Application of the RAILS (Raking-Assisted Integration of Linked Surveys) methodology to the [All of Us (AoU) Research Program](https://allofus.nih.gov/), using 2022 ACS PUMS as the probability reference sample to produce nationally representative prevalence estimates from the AoU biobank.

> AoU participant data are not publicly available. Steps involving AoU data must be run inside the [All of Us Researcher Workbench](https://workbench.researchallofus.org/).

---

## Structure

```
Global RAILS/
├── Data Processing/
│   ├── README.md                ← Guide for PUMS, AoU, and NHIS data preparation
│   ├── 01_PUMS_Prep.R           ← Download and recode PUMS reference data
│   ├── 02_AoU_Prep.R            ← Query and harmonize AoU data (Workbench only)
│   └── 03_NHIS_Prep.R           ← Download and recode NHIS data (hybrid design)
│
└── RAILS Procedure/
    ├── README.md                ← Function reference and analysis guide
    ├── AoU_Fun.R                ← Core functions (fun.nps, fun.lkd, fun.rails.threeway)
    ├── 04_Global_RAILS.R        ← Main analysis script
    ├── 05_Prevalence_Analysis.R ← Weighted disease prevalence + national/region figures
    ├── 06_Phecode_Categories.R  ← Phecode-grouped prevalence figures
    ├── 07_ASCVD_Analysis.R      ← ASCVD risk distribution figures
    └── 08_NHIS_Comparison.R     ← AoU vs NHIS 2020 health-outcome table
```

---

## Workflow Overview

Scripts are numbered in run order. Steps 1–8 are the Global RAILS pipeline; the
Subgroup RAILS analysis continues the sequence at step 9 (see
[Subgroup RAILS](../Subgroup%20RAILS/README.md)).

| Step | Where to run | Script | Output |
|---|---|---|---|
| 1. Prepare PUMS data | Local or any R environment | `01_PUMS_Prep.R` | `dt_agg_pums_v2.csv` |
| 2. Prepare AoU data | AoU Workbench | `02_AoU_Prep.R` | `aou_raking_dt.csv`, `dt_agg_aou_v3.csv` |
| 3. Prepare NHIS data *(hybrid design only)* | Local or any R environment | `03_NHIS_Prep.R` | `dt_agg_nhis_v2.csv` |
| 4. Run Global RAILS | AoU Workbench | `04_Global_RAILS.R` | `global_rails_weights.csv` — weights per participant |
| 5. Prevalence analysis | AoU Workbench | `05_Prevalence_Analysis.R` | Prevalence tables + national/region figures |
| 6. Phecode category figures | AoU Workbench | `06_Phecode_Categories.R` | Phecode-grouped figures |
| 7. ASCVD risk analysis | AoU Workbench | `07_ASCVD_Analysis.R` | ASCVD risk distribution figures |
| 8. NHIS comparison table | AoU Workbench | `08_NHIS_Comparison.R` | AoU vs NHIS table (LaTeX + CSV) |
| 9. Subgroup RAILS | AoU Workbench | `../Subgroup RAILS/09_Sub_RAILS.R` | Subgroup weights |

Step 9 is **reusable for any stratification** — sex, region, race group — by changing the subgroup variable.

See the [Data Processing guide](Data%20Processing/README.md) for Steps 1–3 and the [RAILS Procedure guide](RAILS%20Procedure/README.md) for Steps 4–8 and function documentation.

---

## Covariates

All datasets are harmonized to seven categorical covariates:

| Variable | Levels |
|---|---|
| `agegroup` | 18–24, 25–44, 45–64, 65–74, 75+ |
| `sex` | Female, Male |
| `race_eth` | Hispanic, NH Asian, NH Black, NH White, Others |
| `income` | <35k, 35k–50k, 50k–75k, 75k–100k, >100k |
| `edu` | Less than highschool → College graduate or advanced (5 levels) |
| `homeown` | Own, Rent, Others |
| `region` | Northeast, Midwest, South, West |
