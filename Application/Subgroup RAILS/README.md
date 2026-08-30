# Subgroup RAILS

Applies the RAILS methodology **within levels of a subgroup variable** (e.g. sex, Census region, race group), reusing the Global RAILS machinery unchanged.

Script numbering continues the Global RAILS sequence: steps 1–8 live in [Global RAILS](../Global%20RAILS/README.md), and this folder is **step 9**, run after the Global RAILS pipeline has produced the aggregated cell tables.

> **Run inside the AoU Researcher Workbench.** Copy **both** `Sub_AoU_Fun.R` and `AoU_Fun.R` (from `../Global RAILS/RAILS Procedure/`) into the same directory — the scripts here `source("Sub_AoU_Fun.R")`, which in turn sources `AoU_Fun.R`.

---

## Files

| File | Description |
|---|---|
| `Sub_AoU_Fun.R` | `fun.sub.rails.threeway` — runs the RAILS procedure within **each level** of a subgroup variable and row-binds the results |
| `09_Sub_RAILS_sex.R` | Subgroup by **sex**, restricted to the **Female** stratum (breast cancer study) |
| `09_Sub_RAILS_region.R` | Subgroup by **Census region**, all four regions compared side by side |

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
