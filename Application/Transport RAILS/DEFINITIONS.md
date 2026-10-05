# Transport RAILS: Data, Populations and Preprocessing

Formal definitions for the preliminary analysis (`T01_Transport_VEHSS_region.R`). Everything here describes what the code **actually does**, traced to the stage-1 / stage-2 code of `../Subgroup RAILS/VEHSS_comparison*.R` and to T01.

**Notation follows the Subgroup RAILS manuscript** (`full_text.tex`, §2.1 Notations and §4 Application). Symbols that the manuscript does not have are marked **(new)**. §6 defines every diagnostic; §7 lists preprocessing issues found while writing this down.

> **Names.** In the manuscript, subscript $A$ is the non-probability sample (AoU) and $B$ the probability sample (PUMS). The earlier cohort labels `A_local`, `B_pooled`, `C_external` clashed with this and are **retired**. The three source samples are now **within-region**, **all-region** and **out-of-region**, and their RAILS estimates are **S-RAILS**, **P-RAILS** and **T-RAILS** (§3.3). Code, tables and figures use these names. Cached files from earlier runs are translated on load.

---

## 1. Notation

**From the manuscript:**

| Symbol | Meaning in this analysis |
|---|---|
| $\mathscr U$, $N$ | Target population: US residents aged 40+ in the 50 states and DC (§4.4), of size $N$ |
| $\mathscr U_A$, $n_A$ | Non-probability sample: the AoU analytic sample (§3.1) |
| $\mathscr U_B$, $n_B$ | Probability sample: ACS 2022 PUMS (§2), with survey weights $d_i^B$ |
| $\boldsymbol X$, $\boldsymbol x_i$ | Auxiliary (calibration) variables: age group, sex, education, home ownership, household income, race/ethnicity (§4.1) |
| $z$, $a \in \mathscr Z$, $L = \lvert \mathscr Z \rvert$ | Partitioning variable $z$ = Census region of residence; levels $\mathscr Z = \{\text{NE}, \text{MW}, \text{S}, \text{W}\}$, so $L = 4$ |
| $\mathscr U^a = \{ i \in \mathscr U : z_i = a \}$, $N^a$ | **Target population for region $a$** and its size |
| $\mathscr U_A^a$, $n_A^a$; $\mathscr U_B^a$ | AoU participants / PUMS records in region $a$ |
| $\rho_a = N^a / N$ | Population share of region $a$ |
| $y_{ij} \in \{0,1\}$ | Outcome $j$ (AMD, glaucoma, cataract) for subject $i \in \mathscr U_A$, under the AoU EHR case definition (§3.2) |
| $\boldsymbol T^{pop}$ | Calibration targets (population totals of $\boldsymbol X$ and its interactions) |
| $\widehat d_i^A$, $\widehat w_i^A$ | NPS base weight; calibrated (raked) weight |
| $\widehat\mu^a_{j,m}$ | Hájek-ratio prevalence of outcome $j$ in region $a$ under method $m$ |

**New for transport:**

| Symbol | Meaning |
|---|---|
| $\mathscr U_A^{-a} = \mathscr U_A \setminus \mathscr U_A^a$, $n_A^{-a}$ | AoU participants **outside** region $a$ |
| $\mathscr S$ | **Source sample**: the AoU set whose weights are calibrated to region $a$; $\mathscr S \in \{\mathscr U_A^a,\ \mathscr U_A,\ \mathscr U_A^{-a}\}$ (§3.3) |
| $\boldsymbol T^{pop,a}$ | Calibration targets of region $a$: totals of $\boldsymbol X$ (one-, two-, three-way) over $\mathscr U^a$, estimated from $\mathscr U_B^a$ with $d_i^B$ |
| $\widehat w^{A}_{i}(\mathscr S \to a)$ | Weight of $i \in \mathscr S$ after calibrating $\mathscr S$ to $\boldsymbol T^{pop,a}$ |
| $g$ | Age domain: 40+, 40-64, 65+ |
| $\mu^{a}_{j}$, $\tilde\mu^{a}_{j}$ | Estimand (§5) and external benchmark (VEHSS) for region $a$ |

---

## 2. Probability sample $\mathscr U_B$ (PUMS) and the targets $\boldsymbol T^{pop,a}$

**Source.** ACS 2022 1-year Public Use Microdata Sample, person records, downloaded from the Census API state by state (`vehss_pums_2022_bystate_40plus.csv`, stage 1 of the VEHSS scripts).

**Preprocessing, in order:**

1. Keep persons with `AGEP >= 40`.
2. Map the state (`ST`) to $z$ with the official Census region map (§4.4); DC falls in the South.
3. **Population size.** $N = \sum$ `PWGTP` over **all** 40+ records of the mapped states. This includes group-quarters residents and records with missing covariates.
4. Recode $\boldsymbol X$ (§4.1) and drop records with any missing component. HINCP and TEN are undefined for group-quarters residents, so **group quarters drop out here** (§7).
5. **Survey weights.** For the remaining records, $d_i^B = \text{PWGTP}_i \cdot N \big/ \sum_{\text{complete}} \text{PWGTP}$, so that $\sum_{i \in \mathscr U_B} d_i^B = N$.
6. $\mathscr U_B^a = \{ i \in \mathscr U_B : z_i = a \}$, with $N^a = \sum_{i \in \mathscr U_B^a} d_i^B$ and $\rho_a = N^a / N$.

**Targets.** $\boldsymbol T^{pop,a} = \sum_{i \in \mathscr U_B^a} d_i^B\, \boldsymbol x_i^{(3)}$, where $\boldsymbol x^{(3)}$ is the design vector of all one-, two- and three-way terms of $\boldsymbol X$ (`target_totals()` in T01). In the manuscript's terms, this is the coarse information $\mathscr F^a = \{ \mathscr{CO}(\boldsymbol x_i), i \in \mathscr U^a \}$. The PUMS individual records themselves are used only for the overlap and balance diagnostics.

---

## 3. Non-probability sample $\mathscr U_A$ (All of Us)

### 3.1 Analytic sample $\mathscr U_A$

All conditions are applied in this order (stage 2 of `VEHSS_comparison.R`):

| # | Condition | Source in the CDR |
|---|---|---|
| 1 | **Has EHR data** | `cb_search_person.has_ehr_data = 1` |
| 2 | **Age ≥ 40** on the reference date 2024-08-01 | `person.birth_datetime`; age is **rounded**, see §7 |
| 3 | **State of residence** maps to a region | `person_ext.state_of_residence_source_value` ("PII State: XX") |
| 4 | **All components of $\boldsymbol x_i$ non-missing** after recoding (§4.1) | `person` (race, ethnicity, sex at birth); The Basics survey (education, income, home ownership; concept ids 1585370, 1585375, 1585899, 1585940, 43530593) |

T01 reads the AMD stage-2 file as $\mathscr U_A$. It then joins the glaucoma and cataract outcomes by `person_id` and reports any participant missing from those files.

### 3.2 Outcomes $y_{ij}$

$y_{ij} = 1$ if participant $i$ has **at least 2 distinct dates** with a code of the PhecodeX phecode for outcome $j$ (child phecodes included) **at any time in the EHR**. There is no look-back window, so this is lifetime recorded diagnosis.

| $j$ | PhecodeX | Codes | Tables searched |
|---|---|---|---|
| AMD | SO_374.51 | 3 ICD-9-CM + 46 ICD-10-CM (no H35.30) | `condition_occurrence` |
| Glaucoma | SO_375.1 | glaucoma, suspects excluded | `condition_occurrence` + `observation` |
| Cataract | SO_371 | includes pseudophakia / aphakia / congenital | `condition_occurrence` + `observation` |

A code is matched on the source concept's ICD code or on the raw source value.

### 3.3 Source samples $\mathscr S$ for target region $a$

| Source $\mathscr S$ | Name (code) | Size | Relation to the target $\mathscr U^a$ | RAILS estimate | Old label |
|---|---|---|---|---|---|
| $\mathscr U_A^a$ | **Within-region** (`within`) | $n_A^a$ | $\mathscr U_A^a \subsetneq \mathscr U^a$: drawn from the target, the manuscript's own assumption $\mathscr U_A \subsetneq \mathscr U$ applied within $a$ | **S-RAILS** (generalization) | `A_local` |
| $\mathscr U_A$ | **All-region** (`all`) | $n_A$ | $\mathscr U_A \cap \mathscr U^a = \mathscr U_A^a$: contains the target's participants plus the other regions' | **P-RAILS** (pooled) | `B_pooled` |
| $\mathscr U_A^{-a}$ | **Out-of-region** (`out`) | $n_A - n_A^a$ | $\mathscr U_A^{-a} \cap \mathscr U^a = \varnothing$: no target member. **The manuscript's assumption $\mathscr U_A \subsetneq \mathscr U$ fails for target $\mathscr U^a$** | **T-RAILS** (transport) | `C_external` |

The estimator names apply to the **RAILS** weights only. Unweighted and raking estimates from the same source are named by the source, e.g. "unweighted, out-of-region".

By construction, $\mathscr U_A = \mathscr U_A^a \cup \mathscr U_A^{-a}$.

**Relation to the manuscript's estimators:**
- **$\mathscr S = \mathscr U_A^a$ is S-RAILS:** the source is calibrated to its own region's targets $\boldsymbol T^{pop,a}$, as in the manuscript ($\widehat w^{A,a}_{i,\mathrm S}$).
- **$\mathscr S = \mathscr U_A$ is not G-RAILS.** G-RAILS calibrates $\mathscr U_A$ once to the **national** targets $\boldsymbol T^{pop}$ (region included) and evaluates on $\mathscr U_A^a$. Here $\mathscr U_A$ is calibrated to $\boldsymbol T^{pop,a}$ and **all** of it contributes to the estimate for $a$.

---

## 4. Shared definitions

### 4.1 Auxiliary variables $\boldsymbol X$ (identical coding in $\mathscr U_A$ and $\mathscr U_B$)

| Component | Levels | $\mathscr U_A$ (AoU) | $\mathscr U_B$ (PUMS) |
|---|---|---|---|
| Age group | 40-64, 65-84, 85+ | age on 2024-08-01 | `AGEP` |
| Sex | Female, Male | sex at birth (intersex / skip / prefer not → missing) | `SEX` |
| Education | < high school, some high school, high-school graduate, some college, college graduate+ | The Basics "highest grade" | `SCHL` (<12 / 12-15 / 16-17 / 18-20 / 21+) |
| Home ownership | Own, Rent, Others | The Basics (other arrangement / don't know → Others) | `TEN` (1-2 / 3 / 4) |
| Household income | <35k, 35-50k, 50-75k, 75-100k, >100k | The Basics annual income | `HINCP` |
| Race/ethnicity | Hispanic, NH Asian, NH Black, NH White, Others | `person` race + ethnicity: Hispanic first, then single race; MENA, NHPI, "more than one", "none of these" → Others | `HISP`, then `RACWHT` / `RACBLK` / `RACASN`, in that order |

The partitioning variable $z$ (region) is **not** in $\boldsymbol X$: it is constant within $\mathscr U^a$. This mirrors S-RAILS in the manuscript.

### 4.2 Calibration

Each source $\mathscr S$ is calibrated to $\boldsymbol T^{pop,a}$ with the RAILS procedure of the manuscript (`fun.rails.threeway`):
1. NPS base weights $\widehat d_i^A$ from the pseudo-likelihood with $\mathscr U_B^a$;
2. GVS forward selection of three-way terms;
3. LIFO raking, with constraints $\sum_{i \in \mathscr S} \widehat w^A_i(\mathscr S \to a)\, \boldsymbol x_i = \boldsymbol T^{pop,a}$ (selected terms) and $\sum_{i \in \mathscr S} \widehat w^A_i(\mathscr S \to a) = N^a$.

The fallback is `fun.rails.twoway`. The selected terms may differ across the three $\mathscr S$; they are logged in `transport_fit_log.csv`.

For $\mathscr S = \mathscr U_A^{-a}$ the pseudo-likelihood's premise $\mathscr S \subset \mathscr U^a$ does not hold. The fitted "propensity" is then a source-to-target density ratio up to a constant. Raking still enforces $\boldsymbol T^{pop,a}$, but the analytic RAILS variance does not apply.

### 4.3 Age domains $g$

Domain estimates restrict the region-level weights to the domain's age groups; nothing is re-calibrated. Age group is a component of $\boldsymbol X$, so the domain's weight total equals its PUMS total $N^{a,g}$.

### 4.4 Region map (defines $z$)

The **official Census regions**, identical to the PUMS `REGION` variable, used for $\mathscr U_A$, $\mathscr U_B$ and the VEHSS state benchmark:
- DE, MD, DC and OK are in the **South**;
- the Northeast is CT, ME, MA, NH, NJ, NY, PA, RI, VT;
- AK and HI are in the West.

(Earlier runs used the older pipeline map: DE and MD in the Northeast, OK added to the South, DC excluded. Outputs from before this change are not comparable.)

---

## 5. Estimands, estimators and benchmarks

**Estimand.** The prevalence of the **AoU-defined** outcome $j$ in region $a$ (domain $g$ analogous):

$$\mu^{a}_{j} = \frac{1}{N^a} \sum_{i \in \mathscr U^a} y_{ij}.$$

**Estimator from source $\mathscr S$** (Hájek ratio, as in the manuscript, but summing over the **source** rather than over $\mathscr U_A^a$):

$$\widehat\mu^{a}_{j}(\mathscr S) = \frac{\sum_{i \in \mathscr S} \widehat w^A_i(\mathscr S \to a)\, y_{ij}}{\sum_{i \in \mathscr S} \widehat w^A_i(\mathscr S \to a)}.$$

The unweighted analogue $\widehat\mu^{a}_{j,\mathrm N}(\mathscr S)$ is the sample mean over $\mathscr S$.

**Condition for consistency** (beyond the manuscript's ignorability given $\boldsymbol X$ within the source): the outcome model given $\boldsymbol X$ must be the same in $\mathscr S$ as in $\mathscr U^a$.
- For $\mathscr S = \mathscr U_A^a$ this is the manuscript's own assumption within region $a$.
- For $\mathscr S = \mathscr U_A^{-a}$ it additionally requires $E[y_{ij} \mid \boldsymbol x, z = a] = E[y_{ij} \mid \boldsymbol x, z \ne a]$: no **concept shift** across regions.

**Transport contrast (new).** $\Delta^{a}_{j} = \widehat\mu^{a}_{j}(\mathscr U_A^{-a}) - \widehat\mu^{a}_{j}(\mathscr U_A^{a})$ = T-RAILS − S-RAILS (`delta_out_within` in the code). It does not involve any benchmark. The analogous P-RAILS − S-RAILS is `delta_all_within`. The concept-shift proxy (`het_proxy`) estimates the part of $\Delta^{a}_{j}$ due to concept shift.

**Benchmark comparisons in the code.** `abs_diff_<source>` $= \lvert \widehat\mu^{a}_{j}(\mathscr S) - \tilde\mu^{a}_{j} \rvert$. `d_absdiff_out_within` $= \lvert \text{T-RAILS} - \tilde\mu \rvert - \lvert \text{S-RAILS} - \tilde\mu \rvert$, where a value above 0 means T-RAILS is further from VEHSS.

**External benchmark $\tilde\mu^{a}_{j}$ (VEHSS).** VEHSS measures a different quantity, so $\tilde\mu^{a}_{j} \ne \mu^{a}_{j}$ in general:

| $j$ | VEHSS quantity | Year | Denominator | Comparable domains |
|---|---|---|---|---|
| AMD | **Modeled** prevalence of any AMD, including undiagnosed | 2019 | residents | 40+, 40-64, 65+ |
| Glaucoma | **Modeled** prevalence (40+ derived from age bands) | 2022 | residents | 40+, 40-64, 65+ |
| Cataract | **Diagnosed** in Medicare FFS + Advantage claims | 2022 | Medicare beneficiaries | 65+ only |

The regional benchmark is $\tilde\mu^{a}_{j} = \sum_{s \in a} \text{cases}_s \big/ \sum_{s \in a} \text{pop}_s$ over states $s$ in region $a$, where pop = cases / prevalence (implied population). A state enters only if VEHSS reports every age band of the domain. No benchmark SE is available.

**Consequence.** $\widehat\mu^{a}_{j}(\mathscr S) - \tilde\mu^{a}_{j}$ mixes the weighting error with the definitional gap $\mu^{a}_{j} - \tilde\mu^{a}_{j}$. The gap is common to all three sources, so **$\Delta^{a}_{j}$ is free of it**; absolute agreement with VEHSS is not.

---

## 6. Diagnostics (fig 5, fig 7, `transport_overlap.csv`, `transport_weight_diag.csv`)

All are computed per target region $a$ and source $\mathscr S$. An **X-cell** is one of the $3 \cdot 2 \cdot 5 \cdot 3 \cdot 5 \cdot 5 = 2{,}250$ combinations of the six components of $\boldsymbol X$.

**Before weighting (overlap):**

| Column | Definition | Reading |
|---|---|---|
| `membership_auc` (c-statistic) | Stack $\mathscr U_B^a$ (label 1, weights $d_i^B$) and $\mathscr S$ (label 0, one per participant), each class rescaled to total weight 1. Fit a logistic regression of the label on $\boldsymbol X$ (main effects + all two-way terms). The c-statistic is the AUC: $\Pr(\text{score of a random target member} > \text{score of a random source participant})$ | 0.5 = same covariate mix; 1 = fully separable |
| `support_joint0` | $\sum_{c:\, n_{\mathscr S}(c) = 0} N^a(c) \,/\, N^a$: the share of the target population in X-cells $c$ that occur in $\mathscr U_B^a$ ($N^a(c) > 0$) but have **no participant in the AoU source** $\mathscr S$. Cells empty in PUMS do not count | These people cannot be represented directly; raking spreads their share over neighbouring cells |
| `support_2way_lt5_mean` / `_worst` | The same idea for each pair of covariates: the share of $N^a$ in two-way cells with fewer than 5 source participants; the mean over the 15 pairs, and the worst pair | Sparse two-way margins make raking fragile |

**After weighting:**

| Column | Definition | Reading |
|---|---|---|
| `ess` | Kish $(\sum \widehat w_i)^2 / \sum \widehat w_i^2$ over $\mathscr S$ | Precision left after weighting |
| `cv_w` | SD / mean of the per-person weights | Weight variability |
| `tvd_joint6` | $\tfrac12 \sum_c \lvert p_{\mathscr S}(c) - p_{\mathscr U_B^a}(c) \rvert$ over all X-cells, with the weighted source and target shares | Balance on the full joint distribution, beyond the calibrated terms. 0 = identical |
| `smd_age` | (weighted mean age in $\mathscr S$ − target mean age) / target SD of age | Balance on detailed age, which is not calibrated (only the age group is). 0.1 is the usual threshold |

**Concept shift (`transport_heterogeneity.csv`):** binomial fits of $y_{ij}$ on $\boldsymbol X$ in AoU cells, unweighted.

| Column | Definition |
|---|---|
| `lrt`, `lrt_p` | Likelihood-ratio test of the region-in/out $\times \boldsymbol X$ interactions |
| `het_proxy` | $\sum_{c} N^a(c)\,[\widehat m_{-a}(c) - \widehat m_{a}(c)] \,/\, N^a$, where $\widehat m_{a}$ and $\widehat m_{-a}$ are the outcome models fitted in $\mathscr U_A^a$ and $\mathscr U_A^{-a}$ |

---

## 7. Preprocessing issues found while writing this down

| # | Issue | Where | Effect | Suggested fix |
|---|---|---|---|---|
| 1 | **Age is rounded, not floored** (`round(interval / 1 year)`) | stage 2, all VEHSS scripts | Participants aged 39.5–39.99 are counted as 40; age-group boundaries shift by half a year | `floor()` |
| 2 | **Multiracial coding differs.** PUMS assigns "White alone or in combination" to NH White first, then Black, then Asian; AoU puts "more than one population" in Others | stage 1 vs stage 2 | $\boldsymbol T^{pop,a}$ has fewer Others and more NH White/Black/Asian than the AoU coding implies; affects every race margin | Use PUMS `RAC1P` (single race alone; two or more → Others) |
| 3 | **Group quarters handled implicitly.** Dropped through missing HINCP/TEN, while $N$ includes them | stage 1 + T01 | $\boldsymbol T^{pop,a}$ describes households; $N$ includes institutions | Drop GQ explicitly before computing $N$, or state it |
| 4 | **Income not inflation-adjusted** (`HINCP` without `ADJINC`) | stage 1 | Small shift at band edges | Multiply by `ADJINC / 1e6` |
| 5 | **Reference dates differ.** AoU age as of 2024-08-01, PUMS 2022, VEHSS 2019 / 2022 | — | Minor for $\boldsymbol T^{pop,a}$; adds to the definitional gap | Document |
| 6 | **AMD searches `condition_occurrence` only**, while glaucoma and cataract also search `observation` | stage 2 | Outcome capture differs by $j$ | Align, or document |
| 7 | **EHR-only cohort** (`has_ehr_data = 1`) | stage 2 | $\mathscr U_A$ is AoU *with EHR*, whose completeness may vary by region and health system: a candidate cause of the regional pattern in the results | Add an EHR-depth measure as a covariate or diagnostic |

Issues 1–4 change the inputs of the existing VEHSS pipeline as well, so any fix should be made there (stages 1–2) and the stage files regenerated.
