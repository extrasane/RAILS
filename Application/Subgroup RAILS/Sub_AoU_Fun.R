#### Sub_AoU_Fun.R — Subgroup RAILS wrapper
#### Run inside the AoU Researcher Workbench
####
#### Thin wrapper around fun.rails.threeway: runs the full Global RAILS
#### procedure (two-way base + selected three-way interactions) separately
#### within each level of a subgroup variable. All fixes, speedups, and
#### benchmark methods of fun.rails.threeway apply to every subgroup
#### automatically.
####
#### Requires AoU_Fun.R (fun.nps, fun.lkd, create_v3, fun.rails.threeway)
#### in the same directory — copy it from "../Global RAILS/RAILS Procedure/".

library(dplyr)
library(Matrix)
library(survey)
source("AoU_Fun.R")

##### Subgroup RAILS estimator
### For each level of subgroup_var:
###   1. Subsets both aggregated cell tables to the level.
###   2. Builds the subgroup's PUMS population margins (one-way + two-way +
###      three-way of names_univar) from the PUMS cell subset.
###   3. Calls fun.rails.threeway with nsiz = the subgroup's PUMS total.
### Returns the row-bound cell-level results, tagged with subgroup_run.
### All d_* weight columns are PER-INDIVIDUAL (see fun.rails.threeway).
###
### Notes:
###   - subgroup_var must be an aggregation dimension of BOTH cell tables
###     (i.e. included when the tables were aggregated), and is excluded
###     from names_univar automatically.
###   - Factor levels are kept (no droplevels) so model matrix columns stay
###     identical between the AoU and PUMS subsets.
fun.sub.rails.threeway <- function(
    dt_agg_aou,
    dt_agg_pums,
    subgroup_var = "region",
    names_univar = setdiff(
      c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
      subgroup_var
    ),
    alpha        = 0.05
) {
  if (subgroup_var %in% names_univar) {
    stop("fun.sub.rails.threeway: subgroup_var must not appear in names_univar")
  }
  if (!subgroup_var %in% names(dt_agg_aou) || !subgroup_var %in% names(dt_agg_pums)) {
    stop("fun.sub.rails.threeway: subgroup_var must be a column of both cell tables")
  }

  ## Max cross-tabulation formula for the subgroup-level population margins
  twovars     <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  threevars   <- combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))
  max_formula <- formula(paste0("~", paste(c(names_univar, twovars, threevars), collapse = "+")))

  subgroup_levels <- sort(unique(dt_agg_aou[[subgroup_var]]))

  out <- lapply(subgroup_levels, function(lev) {

    message("=== Subgroup ", subgroup_var, " = ", lev, " ===")

    agg_aou_sub  <- dt_agg_aou  %>% filter(.data[[subgroup_var]] == lev)
    agg_pums_sub <- dt_agg_pums %>% filter(.data[[subgroup_var]] == lev)

    if (nrow(agg_pums_sub) == 0) {
      warning("fun.sub.rails.threeway: no PUMS cells for ", subgroup_var,
              " = ", lev, " — skipped.")
      return(NULL)
    }

    ## Subgroup population margins, computed once from the PUMS cell subset
    mat_max    <- sparse.model.matrix(max_formula, data = agg_pums_sub, keep.order = TRUE)
    pop_totals <- setNames(
      as.numeric(Matrix::crossprod(mat_max, agg_pums_sub$weight)),
      colnames(mat_max)
    )

    fun.rails.threeway(
      dt_agg_aou   = agg_aou_sub,
      dt_agg_pums  = agg_pums_sub,
      pop_totals   = pop_totals,
      names_univar = names_univar,
      alpha        = alpha,
      nsiz         = sum(agg_pums_sub$weight)
    ) %>%
      mutate(subgroup_run = lev)
  })

  bind_rows(out)
}
