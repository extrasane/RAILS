# Transport RAILS: Roadmap

**Status:** application only (no simulation). A **preliminary** run is implemented in [`T01_Transport_VEHSS_region.R`](T01_Transport_VEHSS_region.R); see §0. Items marked **[Decision]** are still open for the full study.

---

## 0. Preliminary implementation (T01)

`T01_Transport_VEHSS_region.R` runs the whole design on data that already exists in the Workbench. It needs no new queries. Formal definitions of every dataset, population, source sample and preprocessing step, in the manuscript's notation, are in [DEFINITIONS.md](DEFINITIONS.md).

> **Naming update.** This roadmap was drafted with cohorts $A_r$ / $B$ / $C_r$. They are now:
> - $A_r$ → **within-region** source $\mathscr U_A^a$ (**S-RAILS**);
> - $B$ → **all-region** source $\mathscr U_A$ (**P-RAILS**);
> - $C_r$ → **out-of-region** source $\mathscr U_A^{-a}$ (**T-RAILS**).
>
> Read the sections below with that mapping; the target region $r$ is $a$ in the manuscript's notation.

| Roadmap element | Preliminary choice |
|---|---|
| Targets | 4 Census regions (pipeline state map + OK) |
| Population | AoU and PUMS 40+, the stage-1 / stage-2 files of the `VEHSS_comparison*.R` scripts |
| Outcomes | AMD, glaucoma, cataract (PhecodeX, ≥ 2 code dates) |
| Benchmark | VEHSS, states aggregated to regions by implied population; domains 40+, 40-64, 65+; cataract only 65+ (Medicare denominator) |
| Calibration set | the 6 pipeline covariates without region |
| Methods | unweighted, one-way raking, two-way raking, two-way NPS + raking, RAILS (three-way, two-way fallback), RAILS trimmed at the 99th percentile |
| Overlap | joint-cell support, two-way support, membership c-statistic |
| Weight diagnostics | ESS, CV, max/median, P99/P1, balance on 5-year age and the full 6-way joint (not calibrated) |
| Heterogeneity | marginal (VEHSS region vs rest) and conditional (AoU in/out LRT, standardized outcome-model gap) |
| Variance | naive linearized (as in the existing application scripts). **Not yet:** bootstrap |
| Size-matched B / C | optional, `N_SIZE_MATCH` |
| Not yet | NHIS / BRFSS benchmarks, survey-based outcomes, alternative calibration sets, states |

---

## 1. Study questions

| # | Question | Main comparison |
|---|---|---|
| Q1 | Does weighting participants from **outside** a region to that region's margins (transport) perform worse than weighting the region's **own** participants (generalization)? | Source cohort C vs A, same method and target |
| Q2 | Does adding outside participants to the local ones help or hurt? | Source cohort B vs A |
| Q3 | Which measurable features predict whether transport succeeds? | Error regressed on overlap, ESS, weight variability, prevalence, heterogeneity, covariate–outcome strength |
| Q4 | How sensitive are Q1–Q3 to analytic choices? | Methods, calibration sets, trimming, benchmarks, harmonization rules |

### Notation

| Symbol | Meaning |
|---|---|
| $r$ | target region |
| $k$ | validation outcome |
| $c \in \{A, B, C\}$ | source cohort |
| $m$ | weighting method |
| $\hat\theta_{rkcm}$ | weighted AoU prevalence |
| $\tilde\theta_{rk}$ | external benchmark prevalence, standard error $s_{rk}$ |
| $X$ | calibration covariates (`region` is excluded: it is constant in the target) |

---

## 2. Design overview

```
for each target region r
├── Target inputs
│   ├── population margins of X       ← ACS / PUMS (region r)
│   └── benchmark prevalence θ̃_rk    ← NHIS / BRFSS / VEHSS (region r)
├── Source cohorts
│   ├── A_r : AoU participants in r                 (generalization)
│   ├── B   : all AoU participants                  (pooled)
│   └── C_r : AoU participants outside r            (transport)
└── for each cohort c
    ├── overlap of c with target r
    └── for each method m
        ├── calibrate c → margins of r
        ├── diagnostics: balance, weight variability, ESS, positivity
        └── for each outcome k
            ├── θ̂_rkcm
            └── compare with θ̃_rk: |D|, ARD, PR, CI(D), coverage
across r, k
├── Q1/Q2: paired cohort contrasts
├── Q3: error ~ predictors
└── Q4: sensitivity re-runs
```

---

## 3. Phase 0 — Definitions and harmonization

This phase must be settled first, because it fixes what "the same target" means across the three data sources.

### 3.1 Target population

| Item | Proposal | Why it matters |
|---|---|---|
| Age | adults 18+ (or 20+ to match NHIS adult file conventions) | AoU, PUMS and the benchmark must cover the same ages |
| Residence | household population only (drop PUMS group quarters) | NHIS and BRFSS exclude institutionalized adults |
| Year | PUMS year aligned with the benchmark year(s) | prevalence and margins drift over time |
| Geography | **[Decision]** 4 Census regions, 9 divisions, or states | sets how many target units Q3 has (see §8.2) |

### 3.2 Calibration covariates

- **Primary:** `agegroup, sex, edu, homeown, income, race_eth` — the Global RAILS set without `region`.
- Category coding must be identical in AoU, PUMS and (for balance checks) the benchmark survey.

### 3.3 Validation outcomes

**[Decision]** Choose 6–12 outcomes. For each, record the definition in all three sources in one harmonization table: `outcome_harmonization.csv`.

| Outcome type | Examples | AoU source | Benchmark |
|---|---|---|---|
| Self-reported diagnosis | diabetes, hypertension, asthma, COPD, depression, CHD, stroke | AoU surveys (Personal / Family Health History) | NHIS, BRFSS ("ever told") |
| Behaviour / status | current smoking, insurance, fair/poor health | AoU surveys (Lifestyle, Overall Health, Basics) | NHIS, BRFSS |
| Measured | obesity (BMI ≥ 30) | AoU physical measurements | BRFSS (self-report), NHANES (measured; no region identifiers in public files) |
| Clinically modeled | AMD, glaucoma, cataract | EHR / PhecodeX | VEHSS (already in `Sub_VEHSS/VEHSS_*`) |

**Primary harmonization rule:** prefer AoU **survey** items when the benchmark is self-reported. EHR depth varies by region, so EHR-based outcomes mix measurement differences into the transport error. EHR definitions become a sensitivity analysis (§7).

**Item non-response:** handle the same way as `08_NHIS_Comparison.R`, where weights are rescaled per outcome among respondents. Double-weighting is a sensitivity analysis.

---

## 4. Phase 1 — Data assembly

| Step | Input | Output |
|---|---|---|
| 1.1 Target margins | PUMS microdata by state → region | one-, two- and three-way margins of $X$ per region, from the **by-state individual** PUMS (avoids the aggregated-PUMS region-joint issue) |
| 1.2 Target microdata | same PUMS records | individual target sample for overlap modeling (§5) |
| 1.3 Benchmarks | NHIS public use (has `REGION`); BRFSS (state → region, BRFSS design weights); VEHSS | $\tilde\theta_{rk}$, $s_{rk}$ per region and outcome, design-based SEs |
| 1.4 Source cohorts | AoU prepped data (`02_AoU_Prep.R`) | $A_r$, $B$, $C_r$ with harmonized $X$ and outcomes |
| 1.5 Benchmark covariates | NHIS / BRFSS covariates mapped to $X$ | used only to check that the benchmark population matches the PUMS margins |

Step 1.5 is a check on the benchmark itself. Large gaps between the benchmark's weighted $X$ distribution and PUMS mean the "truth" refers to a different population.

---

## 5. Phase 2 — Covariate overlap (per region × cohort, before weighting)

| Metric | Definition | Reads as |
|---|---|---|
| Cell support | share of target population (PUMS weight) in two-way cells of $X$ with 0 or < 5 source participants | direct positivity check for calibration |
| Membership c-statistic | AUC of a logistic model for "in target (PUMS) vs in source (AoU)" on $X$ with two-way terms | 0.5 = identical $X$ distributions; → 1 = poor overlap |
| Density-ratio tail | share of target weight with estimated ratio $p_t(x)/p_c(x)$ above the 99th percentile of $A_r$'s ratios | how much of the target relies on few source people |
| Standardized differences | unweighted SMD for each level of $X$ | descriptive |

$A_r$ is the natural reference for $C_r$. Same-region AoU participants already differ from the region's population, so transport should be judged by how much **worse** $C_r$'s overlap is than $A_r$'s.

---

## 6. Phase 3 — Weighting and diagnostics

### 6.1 Methods

The same target margins are used for every cohort.

| Code | Method | Existing function |
|---|---|---|
| UW | unweighted | — |
| R1 | one-way raking | `survey::calibrate` |
| R2 | two-way raking | `survey::calibrate` |
| PS | pseudo-likelihood PS (main effects) | `fun.nps` / `AoU_Fun.R` |
| RAILS-2 | two-way RAILS (forward selection + LIFO raking) | `fun.rails.twoway` |
| RAILS-3 | three-way RAILS | `fun.rails.threeway` |

**Methodological point.**
- The pseudo-likelihood PS step assumes the source is a subsample of the reference population.
- That holds for $A_r$. It holds partly for $B$, and not at all for $C_r$, which is disjoint from the target.
- For $C_r$ the fitted "propensity" is really a source-to-target density ratio, up to a constant. The point estimate is still well defined, because calibration enforces the target margins whatever the starting weights.
- The stacked-estimating-equation variance (`fun.rails.var`) must be **re-derived** for $C_r$ and $B$, or replaced by a bootstrap that resamples AoU participants (§6.3).

### 6.2 Weight diagnostics

| Diagnostic | Definition |
|---|---|
| Convergence | raking converged; LIFO terms dropped; empty-margin warnings |
| ESS | Kish $n_\text{eff} = (\sum w)^2 / \sum w^2$; also as a share of $n_c$ |
| Weight variability | CV of $w$; max/median; 99th/1st percentile ratio |
| Balance on **calibrated** margins | exact by construction; report only as a convergence check |
| Balance on **non-calibrated** features | weighted SMD vs PUMS for (a) joints above the calibrated order (three-way cells for RAILS-2), and (b) auxiliary PUMS variables not used in calibration (e.g. marital status, employment, nativity, rurality where available in AoU) |

Balance has to be judged on features left out of calibration, because calibration matches the margins it used exactly.

### 6.3 Variance

- **Primary:** bootstrap over AoU participants within cohort, re-running the whole weighting procedure, selection included. 200–500 replicates.
- **Secondary:** analytic variance where it is valid ($A_r$).

---

## 7. Phase 4 — Benchmark comparison

For each $(r, k, c, m)$, with $D = \hat\theta - \tilde\theta$:

| Metric | Formula |
|---|---|
| Absolute difference | $\lvert D \rvert$ |
| Absolute relative difference | $\lvert D \rvert / \tilde\theta$ |
| Prevalence ratio | $\text{PR} = \hat\theta / \tilde\theta$; CI on the log scale by the delta method |
| CI for $D$ | $D \pm 1.96\sqrt{\hat V_\text{AoU} + s_{rk}^2}$ (the two samples are independent) |
| Benchmark coverage | **(i)** $\tilde\theta \in \text{CI}_\text{AoU}$ (AoU uncertainty only); **(ii)** $0 \in \text{CI}_D$ (both uncertainties). Report both. |

For Q1, the **signed** contrast $\hat\theta_C - \hat\theta_A$ does not involve the benchmark, so benchmark error cancels. Absolute errors do not cancel. Report both views.

---

## 8. Phase 5 — Cross-region synthesis

### 8.1 Q1 / Q2: cohort comparisons

- Paired contrasts over $(r, k)$, with the same method and the same target:
  - $\Delta^{CA} = \lvert D_C \rvert - \lvert D_A \rvert$
  - $\Delta^{BA} = \lvert D_B \rvert - \lvert D_A \rvert$
- Summaries: mean and median of the contrasts, share of pairs where C beats A, and a bootstrap CI.
- **Size confound:** $C_r$ is roughly 3× larger than $A_r$. Add a **size-matched** run that subsamples $C_r$ (and $B$) to $n_{A_r}$, 100 draws. This separates "outside participants are different" from "outside participants are more numerous".

### 8.2 Q3: predictors of transport success

Analysis unit: $(r, k, c)$, or the contrast $\Delta^{CA}_{rk}$. Fit a linear mixed model with random intercepts for region and outcome.

| Predictor | Measure |
|---|---|
| Covariate overlap | membership c-statistic, cell-support share (§5) |
| Effective sample size | $\log n_\text{eff}$ |
| Weight variability | CV of weights |
| Outcome prevalence | $\tilde\theta_{rk}$ (on the logit scale) |
| Regional outcome heterogeneity — marginal | $\lvert \tilde\theta_{rk} - \tilde\theta_{-r,k} \rvert$ from benchmarks |
| Regional outcome heterogeneity — **conditional** | within AoU: LRT for `region × X` terms in the outcome model, and the **internal transport-bias proxy** $\hat E_{t}[\hat m_{-r}(X) - \hat m_r(X)]$, averaged over target-$r$ margins |
| Covariate–outcome strength | AUC / pseudo-$R^2$ of outcome on $X$ in the source; and $\lvert \hat\theta_\text{RAILS} - \hat\theta_\text{UW} \rvert$, which measures how far calibration moves the estimate |

The conditional heterogeneity proxy measures outcome heterogeneity given $X$ (the subgroup simulations' bit 2) in real data. It is the quantity that drives transport bias under good selection modeling. If it predicts $\Delta^{CA}$, it can serve as an a-priori warning that needs no benchmark.

**Power caveat:** with 4 regions and about 10 outcomes there are only about 40 units, and they are correlated. With 4 regions, Q3 is exploratory. **[Decision]** Divisions (9), or states with BRFSS benchmarks (50 + DC), would make Q3 testable. States would also reuse the state-level work in `10_State_Prevalence_Maps.R`.

---

## 9. Phase 6 — Sensitivity analyses (Q4)

| Factor | Variants |
|---|---|
| Weighting method | all of §6.1; plus the S-RAILS fallback rule from the VEHSS scripts |
| Calibration set | (a) primary 6; (b) demographics only (`agegroup, sex, race_eth`); (c) primary + an extra SES or health-access variable if it exists in PUMS and AoU; (d) model order two-way vs three-way |
| Weight trimming | cap at the 95th / 99th percentile, or at 5 × median; re-rake after trimming |
| Benchmarks | NHIS vs BRFSS vs VEHSS where outcomes overlap; pooled benchmark years |
| Harmonization | AoU survey vs EHR/phecode definitions; strict vs broad code lists (the `PHENO_DEF` switch pattern) |
| Cohort size | size-matched $B$ / $C$ (§8.1) |
| Population definition | 18+ vs 20+; with and without group quarters |

---

## 10. Expected pattern (guides interpretation; no simulation planned)

From the exact population limits worked out in the design discussion:
- If the outcome model given $X$ is the same across regions, $C$ is unbiased for any covariate shift.
- $C$'s bias grows with outcome heterogeneity given $X$. The conditional heterogeneity proxy (§8.2) estimates exactly this quantity.
- $A$ stays approximately unbiased, but its variance grows as the local cohort shrinks.

---

## 11. Proposed file layout

```
Application/Transport RAILS/
├── ROADMAP.md                      ← this file
├── T01_Transport_VEHSS_region.R    ← PRELIMINARY: whole design on the VEHSS eye outcomes (§0)
├── README.md                       ← written once the full study code exists
├── Trans_Fun.R                     ← cohort builders, overlap, diagnostics, metrics, bootstrap
├── T01_Target_Margins.R            ← PUMS by state → region margins + target microdata
├── T02_Benchmarks.R                ← NHIS / BRFSS / VEHSS region prevalences + SEs
├── T03_Source_Cohorts.R            ← A_r, B, C_r with harmonized X and outcomes
├── T04_Overlap.R                   ← §5
├── T05_Weighting.R                 ← §6.1–6.3 for every (r, c, m)
├── T06_Benchmark_Comparison.R      ← §7
├── T07_Synthesis.R                 ← §8 tables, contrasts, mixed models, figures
├── T08_Sensitivity.R               ← §9 driver over variant grid
└── outcome_harmonization.csv       ← §3.3
```

Every T-script that touches AoU data runs in the Researcher Workbench. It sources `AoU_Fun.R` and `Sub_AoU_Fun.R` unchanged.

---

## 12. Planned outputs

| Output | Content |
|---|---|
| Table 1 | Target populations and benchmarks: region sizes, benchmark sample sizes, outcome prevalences |
| Table 2 | Overlap and weight diagnostics by region × cohort (primary method) |
| Table 3 | Benchmark comparison: $\lvert D \rvert$, ARD, PR, coverage by cohort × method, averaged over region and outcome |
| Figure 1 | Forest plot of $D$ (with CI) by outcome, faceted by region, colored by cohort |
| Figure 2 | $\Delta^{CA}$ against the conditional heterogeneity proxy |
| Figure 3 | Error against ESS and overlap |
| Supplement | Sensitivity grid |

---

## 13. Open decisions (to settle before coding)

| # | Decision | Options | Default proposal |
|---|---|---|---|
| D1 | Target geography | 4 regions / 9 divisions / states | regions (primary) + states via BRFSS (Q3) |
| D2 | Outcome list | §3.3 | 8–10 self-reported outcomes + the 3 VEHSS eye outcomes |
| D3 | Primary benchmark | NHIS / BRFSS | NHIS for regions, BRFSS for states |
| D4 | Years | PUMS and benchmark alignment | PUMS 2022 with NHIS 2022 (or pooled 2021–2023) |
| D5 | Target population | age floor, group quarters | 18+, household population |
| D6 | Primary method | §6.1 | RAILS-3 (manuscript method); RAILS-2 when the cohort is small |
| D7 | Variance | bootstrap / analytic | bootstrap primary |
| D8 | Order of work | — | **Decided:** application only; preliminary VEHSS run first (T01), then Phase 0 for the full outcome set |
