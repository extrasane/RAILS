# Subgroup RAILS

Applies the RAILS methodology **within levels of a subgroup variable** (e.g. sex, Census region, race group), reusing the Global RAILS machinery unchanged.

Script numbering continues the Global RAILS sequence: steps 1–8 live in [Global RAILS](../Global%20RAILS/README.md), and this folder is **step 9**, run after the Global RAILS pipeline has produced the aggregated cell tables.

> **Run inside the AoU Researcher Workbench.** Copy **both** `Sub_AoU_Fun.R` and `AoU_Fun.R` (from `../Global RAILS/RAILS Procedure/`) into the same directory — the scripts here `source("Sub_AoU_Fun.R")`, which in turn sources `AoU_Fun.R`.

---

## Files

| File | Description |
|---|---|
| `Sub_AoU_Fun.R` | `fun.sub.rails.threeway` — runs the RAILS procedure within **each level** of a subgroup variable and row-binds the results. Also `fun.rails.twoway` / `fun.sub.rails.twoway`: the same procedure one order lower (one-way base model, forward selection of **two-way** interactions, LIFO stepwise raking), for smaller non-probability samples such as the eye-care-restricted AoU cohort; same output columns |
| `09_Sub_RAILS_sex.R` | Subgroup by **sex**, restricted to the **Female** stratum (breast cancer study) |
| `09_Sub_RAILS_region.R` | Subgroup by **official Census region** (as PUMS `REGION`: DE, MD, DC and OK in the South), all four regions compared side by side |
| `10_State_Prevalence_Maps.R` | State-level unweighted / G-RAILS / S-RAILS prevalence per phenotype, the D (population-weighted Jensen-Shannon divergence) and E (dispersion entropy) metrics, the five phenotype-set choropleths and the D-vs-E scatter for the manuscript. **Sex-specific phenotypes** (conditions anatomically restricted to one sex, since every prevalence uses the full-population denominator) are dropped when `EXCLUDE_SEX_SPECIFIC <- TRUE`; `SEX_SPECIFIC_RULE = "icd"` flags a cause whose ICD-9/10 ranges all fall in the male-genital, female-genital or pregnancy blocks (from `all_icd_codes.csv`), `"name"` uses cause-name patterns (automatic fallback if the lookup is absent); `EXCLUDE_BREAST_CANCER` (default `FALSE`) adds breast cancer. The ICD rule compares codes at full sub-code granularity and requires **all** of a cause's *diagnostic* ranges to be in blocks of the **same** sex; supplementary codes (ICD-10 Z, ICD-9 V/E: screening, history, family history) are ignored (`IGNORE_SUPPLEMENTARY_CODES`), otherwise e.g. prostate cancer's Z12.5 / Z85.46 / Z80.42 would make it "mixed". Partly sex-specific or male+female aggregates are "mixed" and kept; every cause's class (M / F / mixed / none) and its per-range sexes are written to `sex_specific_icd_classification*.csv` for review, and hand decisions go in `SEX_SPECIFIC_MANUAL_EXCLUDE` (default: "Maternal and neonatal disorders", whose neonatal codes cannot occur in adults) / `SEX_SPECIFIC_MANUAL_KEEP`. Outputs of the exclusion run get the `_nosex` suffix (maps, `plot_D_vs_E_nosex.png`, `state_*_nosex.csv`, plus `sex_specific_phenotypes_excluded_nosex.csv` listing what was dropped), so the original all-phenotype outputs from a `FALSE` run are kept side by side. Also builds the manuscript's **D-E regime examples** table (`DE_regime_examples*.csv`, `plot_D_vs_E_regimes*.png`): the four (high/low D) x (high/low E) regimes of the interpretation table cut at quantiles (`REGIME_*`), `N_EXAMPLES` phenotypes closest to each corner, and each phenotype's **leverage subgroup**: the state with the largest share q_a = rho_a d_a / D of its divergence, reported as the leverage state (`leverage_state`) when that share is a majority (`LEVERAGE_SHARE_MIN = 0.5`) and as `spread` otherwise; the share is the only criterion. `map_DE_regimes_{SG,SU,GU}*.png` show the same examples as choropleths of the S-RAILS / G-RAILS, S-RAILS / unweighted and G-RAILS / unweighted state ratios, one row per regime (high D high E, high D low E, low D high E, low D low E), the examples across the row, each panel titled with the phenotype (row 1) and the leverage subgroup with its share and S/G ratio (row 2, e.g. `[CT ~95% D, S/G = 3.16]`, or `[spread; NY ~33% D, S/G = 0.94]` below the majority cut). Colour is a symmetric log scale squished to [1/`MAP_RATIO_CAP`, `MAP_RATIO_CAP`]. D and E themselves compare S-RAILS with G-RAILS symmetrically (Jensen-Shannon); the unweighted estimate enters only the SU / GU maps |
| *VEHSS studies* | The VEHSS benchmark downloads and all `VEHSS_comparison*.R` scripts (AMD, glaucoma, cataract, eye-care cohort, age group) now live in [`../Sub_VEHSS/`](../Sub_VEHSS/README.md), separate from the original subgroup pipeline. They use `Sub_AoU_Fun.R` from this folder |

Both scripts source `Sub_AoU_Fun.R`, drop the subgroup variable from the model (it is constant within a stratum), and scale each subgroup's weights to that subgroup's **own** PUMS total. They differ only in the subgroup variable and whether the data is pre-filtered:

| | `09_Sub_RAILS_sex.R` | `09_Sub_RAILS_region.R` |
|---|---|---|
| `subgroup_var` | `"sex"` | `"region"` |
| `names_univar` | 6 vars incl. `region`, excl. `sex` | 6 vars incl. `sex`, excl. `region` |
| Pre-filter | yes → `sex == "Female"` (one level) | no → runs all four regions |
| Output | `dt_sub_aou_sex_female.csv` | `dt_sub_aou_region.csv` |

To run **both** sexes instead of female-only, drop the `filter(sex == "Female")` lines in `09_Sub_RAILS_sex.R` and pass the full tables — the wrapper then loops both levels, exactly like the region script.

---

## `fun.sub.rails.threeway(dt_agg_aou, dt_agg_pums, subgroup_var, names_univar, alpha)`

For each level of `subgroup_var`, the wrapper:

1. Subsets both aggregated cell tables to the level (factor levels kept, so model matrix columns match between the AoU and PUMS subsets).
2. Builds that subgroup's population margins — all one-way, two-way, and three-way margins of `names_univar` — from the PUMS cell subset.
3. Calls `fun.rails.threeway` with `nsiz` = the subgroup's PUMS total.
4. Tags results with `subgroup_run` and row-binds all levels.

| Argument | Default | Description |
|---|---|---|
| `dt_agg_aou` / `dt_agg_pums` | — | Aggregated cell tables; `subgroup_var` **must** be one of their aggregation dimensions |
| `subgroup_var` | `"region"` | Stratification variable; automatically excluded from `names_univar` |
| `names_univar` | the seven shared covariates minus `subgroup_var` | Main-effect variables within each subgroup |
| `alpha` | `0.05` | Forward LRT threshold, passed through |

Output has one row per covariate cell per subgroup, with all of `fun.rails.threeway`'s columns (`d_unweighted`, `d_cal1`, `d_cal2`, `d_nps1`, `d_nps2`, `d_nps1_rake`, `d_nps2_rake`, `d_rails`, `selected_terms`, `calibrated_terms`) plus `subgroup_run`. Selection and the LIFO walk run independently within each subgroup, so the selected and calibrated terms can differ across levels.

### Example

```r
result_sub <- fun.sub.rails.threeway(
  dt_agg_aou   = dt_agg_aou,      # aggregated WITH the subgroup variable as a dimension
  dt_agg_pums  = dt_agg_pums,
  subgroup_var = "region"
)

result_sub %>%
  group_by(subgroup_run) %>%
  summarise(n_cells = n(),
            total_weight = sum(d_rails, na.rm = TRUE),
            calibrated = first(calibrated_terms))
```

---

## The two analysis scripts

Both load the aggregated tables + individual-level AoU data, harmonize factors, call `fun.sub.rails.threeway`, join the per-individual weights back to participants (by the model covariates, plus `region` for the region run), and save. The individual join key differs:

- **`09_Sub_RAILS_sex.R`** joins on the 6 covariates only (data already female).
- **`09_Sub_RAILS_region.R`** joins on the 6 covariates **plus `region`**, since each participant matches their own region's cell.

**Small-subgroup caveat:** with fewer cells per stratum, empty AoU margins are more common than in the global run — expect more skipped selection candidates and earlier LIFO stopping, all reported via warnings and `calibrated_terms`.

**Region-joints caveat:** the region run calibrates within each region using the **aggregated** `dt_agg_pums`, so its correctness depends on that file's *within-region* joint distributions being right. If you suspect residual region-joint distortion in the aggregated PUMS, calibrate region from **by-state** individual PUMS instead (state → region map) — that was the motivation behind the earlier by-state region variant.
