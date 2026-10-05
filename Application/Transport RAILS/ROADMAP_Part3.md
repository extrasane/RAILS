# Transport RAILS Part 3: NHIS / NHANES Benchmarks

**Status:** design agreed; no code yet. Builds on Part 1 (`T01`, VEHSS benchmark) and Part 2 (`T02`, phenotype panel without a benchmark). Notation as in [DEFINITIONS.md](DEFINITIONS.md).

---

## 1. Aim

Re-assess the transportability of RAILS weights (S-, P- and T-RAILS) against a **survey benchmark with design-based uncertainty**, for both **prevalences** and **demographic associations**. The design also directly tests the hypothesis raised by Parts 1–2: that the region-wide transport contrasts arise from **regional differences in EHR recording**.

## 2. Decisions

| # | Decision | Choice |
|---|---|---|
| D1 | AoU outcome definition | **Both**: EHR (PhecodeX code list per concept) **and** survey self-report (Personal and Family Health History, PFHH) |
| D2 | Region map | **Official Census regions** for AoU, PUMS and NHIS: DE, MD and DC in the South; OK in the South. S-RAILS therefore differs slightly from the manuscript's (only these states move) |
| D3 | NHIS years | **Pooled 2021–2023** adult files, centred on PUMS 2022. Weights divided by 3; strata / PSU kept per year |
| D4 | NHANES | **National definition-gap check only**. Public files have no geography, so regional NHANES would require an NCHS Research Data Center |
| — | Population margins | PUMS 2022, adults 18+, as in Part 2 (original S-RAILS age groups), re-mapped to the official regions |

## 3. Roles of the data sources

| Source | Quantity | Level | Uncertainty |
|---|---|---|---|
| PUMS 2022 | $\boldsymbol T^{pop,a}$ (margins of $\boldsymbol X$) | region | treated as fixed |
| NHIS 2021–23 | $\tilde\mu^a_j$ (self-reported diagnosis); $\tilde\beta^a_{j}$ (demographic associations) | region | design-based (strata, PSU, weights) |
| NHANES 2017–Mar 2020 (and/or 2021–23) | measured / diagnosed / self-reported prevalence | national | design-based |
| AoU | S-, P-, T-RAILS with EHR and survey outcomes | region | naive now; bootstrap later |

## 4. Outcome crosswalk (to be built and graded)

| Concept | NHIS item* | AoU EHR (PhecodeX)* | AoU survey (PFHH)* | NHANES measured* | Alignment |
|---|---|---|---|---|---|
| Hypertension | `HYPEV_A` | hypertension | "ever told" hypertension | BP ≥ 130/80 or medication | to grade |
| High cholesterol | `CHLEV_A` | hyperlipidemia | high cholesterol | total / LDL cholesterol | to grade |
| Diabetes | `DIBEV_A` | type 1/2 diabetes | diabetes | HbA1c ≥ 6.5% / glucose | to grade |
| Coronary heart disease | `CHDEV_A` | ischemic heart disease | CHD | — | to grade |
| Myocardial infarction | `MIEV_A` | MI | heart attack | — | to grade |
| Angina | `ANGEV_A` | angina | angina | — | to grade |
| Stroke | `STREV_A` | cerebrovascular disease | stroke | — | to grade |
| Asthma (ever) | `ASEV_A` | asthma | asthma | — | to grade |
| COPD / emphysema / chronic bronchitis | `COPDEV_A` | COPD | COPD | — | to grade |
| Any cancer (excl. skin?) | `CANEV_A` | malignant neoplasms | cancer | — | to grade |
| Arthritis | `ARTHEV_A` | arthropathies | arthritis | — | to grade |
| Depression | `DEPEV_A` | depressive disorders | depression | PHQ-9 ≥ 10 | to grade |
| Anxiety | `ANXEV_A` | anxiety disorders | anxiety | — | to grade |
| Chronic kidney disease | `KIDWEAKEV_A` | CKD | kidney disease | eGFR < 60 / ACR ≥ 30 | to grade |
| Obesity | BMI (self-report) | — | — | measured BMI ≥ 30 | to grade |

\*Variable names and code lists are **provisional**: NHIS names must be verified against each year's codebook, and PFHH question / answer concepts against the AoU CDR (step T03-0 below).

**Grading.**
- **close:** same construct and wording;
- **approximate:** same construct, different scope, e.g. NHIS "COPD, emphysema or chronic bronchitis" vs a narrower code list;
- **poor:** excluded from the main analysis.

## 5. Analyses

1. **Prevalence vs NHIS, by region.** For each source (S / P / T), method and AoU definition (EHR / survey):
   - difference $D = \widehat\mu - \tilde\mu$, with CI $D \pm 1.96\sqrt{\hat V_{AoU} + \hat V_{NHIS}}$;
   - prevalence ratio and benchmark coverage;
   - transport contrast $\Delta^a_j$ (benchmark-free).
2. **Ascertainment test.** Compare the regional pattern of $\Delta^a_j$ and of $D$ under the EHR vs the survey definition. EHR-driven recording differences predict:
   - a persistent region-wide sign pattern for EHR;
   - an attenuated pattern for survey self-report;
   - S-RAILS (EHR) deviating from NHIS in the Parts 1–2 direction (Northeast / Midwest high, South / West low).
3. **Associations vs NHIS.** Survey-weighted logistic regressions in NHIS by region (conditional and marginal, same exposures and reference levels as Part 2), against $\beta$ from S-, P- and T-RAILS.
4. **Benchmark population check.** NHIS weighted distribution of $\boldsymbol X$ (age, sex, race/ethnicity, education, family income, tenure) against the PUMS margins by region. NHIS income is *family* income, not household income; this must be documented.
5. **National definition gap (NHANES).** For diabetes, hypertension, CKD, high cholesterol, depression and obesity: measured vs self-reported vs AoU (G-RAILS, national). This quantifies the expected EHR / self-report gap independently of weighting.
6. **PFHH non-response.** Survey outcomes are available only for PFHH respondents. Use an item-response adjustment (as the double weighting in `08_NHIS_Comparison.R`) and report the response rate by region.

## 6. Build order

| Step | Script | Where | Content |
|---|---|---|---|
| T03-0 | `T03_0_AoU_PFHH_explore.R` | Workbench | List PFHH question / answer concepts and counts, so the survey column of the crosswalk can be fixed |
| T03-1 | `T03_1_NHIS_prep.R` | local (internet) | Download NHIS adult 2021–2023; verify variables; harmonize $\boldsymbol X$ and outcomes; official regions; pooled weights; regional prevalences and $\beta$ with design SEs |
| T03-2 | `T03_2_NHANES_prep.R` | local (internet) | National measured / self-reported prevalences |
| T03-3 | `T03_3_AoU_outcomes.R` | Workbench | EHR (PhecodeX lists) and PFHH outcomes for the crosswalk concepts; official region map |
| T03-4 | `T03_4_Transport_NHIS.R` | Workbench | Calibration (official regions), estimates, comparisons, ascertainment test, figures |
| T03-5 | report | local | Rmd / HTML as for Parts 1–2 |

## 7. Known limitations

- **NHIS outcomes are self-reported.** Agreement with AoU-survey is expected to be closer than with AoU-EHR; neither equals true prevalence (see the NHANES check).
- **Covariate definitions differ:** family income in NHIS vs household income in PUMS and AoU.
- **NHIS regional sample sizes** limit precision for rare conditions and sparse covariate levels.
- **Region map change** breaks exact reproduction of the manuscript's S-RAILS.
