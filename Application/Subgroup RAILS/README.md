# Subgroup RAILS

Applies the RAILS methodology **within levels of a subgroup variable** (e.g. sex, Census region, race group), reusing the Global RAILS machinery unchanged.

Script numbering continues the Global RAILS sequence: steps 1–8 live in [Global RAILS](../Global%20RAILS/README.md), and this folder is **step 9**, run after the Global RAILS pipeline has produced the aggregated cell tables.

> **Run inside the AoU Researcher Workbench.** Copy `AoU_Fun.R` from `../Global RAILS/RAILS Procedure/` into the same directory — the scripts here load it via `source("AoU_Fun.R")`.

---

## Files

| File | Description |
|---|---|
| `Sub_AoU_Fun.R` | `fun.sub.rails.threeway` — runs **every level** of a subgroup variable in one call |
| `09_Sub_RAILS.R` | Runs **one stratum**, set via `SUB_VAR` / `SUB_LEVEL` at the top of the script |

### Which one to use

- **One stratum of interest** (e.g. females only, for a breast cancer analysis) → `09_Sub_RAILS.R`. Change the three settings lines and re-run; nothing else in the script is subgroup-specific.
- **Every level at once** (e.g. all four Census regions, compared side by side) → `fun.sub.rails.threeway` from `Sub_AoU_Fun.R`, which loops the levels and row-binds the results with a `subgroup_run` tag.

Both scale weights to the **subgroup's own** PUMS total and drop the subgroup variable from the model, since it is constant within the stratum.

```r
## 09_Sub_RAILS.R — the only lines to change between analyses
SUB_VAR   <- "sex"        # or "region", "race_eth", ...
SUB_LEVEL <- "Female"     # or "South", "NH Black", ...
OUT_FILE  <- "dt_sub_aou_femaleonly.csv"
```

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

## Single-Stratum Script (`09_Sub_RAILS.R`)

Loads the aggregated tables and the individual-level AoU data, filters all three to `SUB_VAR == SUB_LEVEL`, drops the subgroup variable from `names_univar`, builds the subgroup's margins, calls `fun.rails.threeway` directly, joins the weights back to participants, and saves.

The script reports the subgroup's cell count, participant count, and `nsiz` before the fit — worth checking before a long selection run, since small strata are where sparse cells cause early LIFO stops.

**Small-subgroup caveat:** with fewer cells per stratum, empty AoU margins are more common than in the global run — expect more skipped selection candidates and earlier stopping, all reported via warnings and `calibrated_terms`.
