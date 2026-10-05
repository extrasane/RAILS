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


########################################################################
## TWO-WAY variants — for a SMALLER non-probability sample (e.g. the AoU
## cohort restricted to eye-care contact). Same procedure as the three-way
## functions above, one order lower: the base model is the ONE-WAY (main
## effects) model, the forward LRT selection runs over the TWO-WAY
## interactions, and the LIFO stepwise raking adds the selected two-way terms
## one at a time. The existing functions are not changed.
########################################################################

##### RAILS estimator, one-way base + selected TWO-WAY interactions
### Output columns are the same as fun.rails.threeway so downstream joins do
### not change:
###   d_unweighted, d_cal1, d_nps1, d_nps1_rake   — as in fun.rails.threeway
###   d_cal2, d_nps2, d_nps2_rake                 — FULL two-way benchmarks, NA
###                                                 (with a warning) when the
###                                                 full two-way model cannot be
###                                                 fitted on the small sample
###   d_rails                                     — RAILS weights (one-way base +
###                                                 selected two-way terms)
###   selected_terms / calibrated_terms           — one-way terms + selected /
###                                                 calibrated two-way terms
### pop_totals must contain at least the one-way and two-way margins.
fun.rails.twoway <- function(
    dt_agg_aou,
    dt_agg_pums,
    pop_totals,
    names_univar = c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
    alpha        = 0.05,
    nsiz         = sum(dt_agg_pums$weight)
) {
  ## One-way base model and its likelihood
  cat_oneway <- formula(paste0("~", paste(names_univar, collapse = "+")))

  d_init <- fun.nps(cat_oneway, dt_agg_aou, dt_agg_pums, nsiz = nsiz)
  m_pums <- sparse.model.matrix(cat_oneway, dt_agg_pums)
  m_aou  <- sparse.model.matrix(cat_oneway, dt_agg_aou)
  theta1 <- attr(d_init, "theta")
  LKHD1  <- sum(as.numeric(crossprod(m_aou, dt_agg_aou$weight)) * theta1) -
    sum(dt_agg_pums$weight * log1p(exp(as.numeric(m_pums %*% theta1))))

  ## Forward two-way selection
  vars_new         <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  names_twoway_new <- names_univar
  LKHD_new         <- LKHD1
  m_s_new          <- m_pums

  index_add <- 1
  while (index_add != 0) {
    if (length(vars_new) == 0) break
    temp_p <- apply(
      as.data.frame(vars_new), 1, fun.lkd,
      names_univar = names_twoway_new, LKHD = LKHD_new,
      dt_s = dt_agg_pums, dt_aou = dt_agg_aou, m_s = m_s_new
    )
    colnames(temp_p) <- vars_new
    temp_sign <- temp_p[, !is.na(temp_p[3, ]) & temp_p[3, ] < alpha, drop = FALSE]
    if (ncol(temp_sign) == 0 || all(is.na(temp_sign))) {
      index_add <- 0
    } else {
      selected_var     <- colnames(temp_sign)[which.max(temp_sign[4, ])]
      message(format(Sys.time(), "%H:%M:%S"), "  GVS (two-way): selected '", selected_var, "'")
      names_twoway_new <- c(names_twoway_new, selected_var)
      LKHD_new         <- temp_sign[1, selected_var]
      m_s_new <- sparse.model.matrix(
        as.formula(paste0("~", paste(names_twoway_new, collapse = "+"))),
        dt_agg_pums
      )
      vars_new <- vars_new[vars_new != selected_var]
    }
  }

  cat_formula <- formula(paste0("~", paste(names_twoway_new, collapse = "+")))

  ## Subset pop_totals ONCE to the selected model's margins
  pop_totals_sel <- create_v3(pop_totals, colnames(sparse.model.matrix(cat_formula, dt_agg_aou)))

  ## LIFO stepwise raking — ascending, stop at first failure
  n_base <- length(names_univar)
  n_max  <- length(names_twoway_new)

  weights_rails <- rep(NA_real_, nrow(dt_agg_aou))
  n_used        <- 0
  theta_prev    <- theta1

  if (n_max == n_base) warning("fun.rails.twoway: no two-way term selected — RAILS weights are NA.")

  n_uni <- n_base
  while (n_uni < n_max) {
    n_uni <- n_uni + 1
    message(format(Sys.time(), "%H:%M:%S"), "  LIFO step ", n_uni - n_base, "/",
            n_max - n_base, ": adding '", names_twoway_new[n_uni], "'")
    cat_temp <- formula(paste0("~", paste(names_twoway_new[seq_len(n_uni)], collapse = "+")))
    temp_d   <- fun.nps(cat_temp, dt_agg_aou, dt_agg_pums, nsiz = nsiz, theta_init = theta_prev)
    theta_prev <- attr(temp_d, "theta")

    temp <- dt_agg_aou %>%
      mutate(count = temp_d * weight, count = count / sum(count) * nsiz)

    clus     <- svydesign(id = ~1, weights = ~count, data = temp)
    names_ps <- cal_names(cat_temp, clus)
    vars_cal <- create_v3(pop_totals_sel, names_ps)

    step_ok <- tryCatch({
      clus_cal <- survey::calibrate(
        clus, cat_temp, vars_cal,
        calfun = "raking", maxit = 1e3,
        epsilon = rep(1e-7, length(vars_cal))
      )
      weights_rails <- weights(clus_cal)
      n_used        <- n_uni
      TRUE
    }, error = function(e) FALSE)

    if (!step_ok) {
      warning(
        "fun.rails.twoway: covariate combination '", names_twoway_new[n_uni],
        "' fails to pile up to higher order — keeping the previous model."
      )
      break
    }
  }

  ## -------------------------------------------------------------------
  ## Benchmark methods on the same cells (same columns as the three-way
  ## version). The FULL two-way NPS model may be singular on a small
  ## sample: then d_nps2 / d_nps2_rake are NA with a warning.
  ## -------------------------------------------------------------------
  twovars    <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  cat_twoway <- formula(paste0("~", paste(c(names_univar, twovars), collapse = "+")))

  d_nps1_cell <- {
    ct <- d_init * dt_agg_aou$weight        # one-way NPS already fitted above
    ct / sum(ct) * nsiz
  }
  d_nps2_cell <- tryCatch({
    d2 <- fun.nps(cat_twoway, dt_agg_aou, dt_agg_pums, nsiz = nsiz)
    ct <- d2 * dt_agg_aou$weight
    ct / sum(ct) * nsiz
  }, error = function(e) {
    warning("fun.rails.twoway: full two-way NPS benchmark not estimable (", conditionMessage(e), ")")
    rep(NA_real_, nrow(dt_agg_aou))
  })
  d_eq_cell <- dt_agg_aou$weight / sum(dt_agg_aou$weight) * nsiz

  ## Rake a cell design with the given starting cell totals to the margins
  ## of cat_f; NA (with warning) on non-convergence or missing start weights.
  fun.cal <- function(cat_f, start_cell) {
    if (any(is.na(start_cell))) return(rep(NA_real_, nrow(dt_agg_aou)))
    temp_c <- dt_agg_aou %>% mutate(count = start_cell)
    clus_c <- svydesign(id = ~1, weights = ~count, data = temp_c)
    vars_c <- create_v3(pop_totals, cal_names(cat_f, clus_c))
    tryCatch(
      weights(survey::calibrate(
        clus_c, cat_f, vars_c,
        calfun = "raking", maxit = 1e3,
        epsilon = rep(1e-7, length(vars_c))
      )),
      error = function(e) {
        warning("fun.rails.twoway: benchmark raking failed to converge (",
                conditionMessage(e), ")")
        rep(NA_real_, nrow(dt_agg_aou))
      }
    )
  }
  d_cal1_cell      <- fun.cal(cat_oneway, d_eq_cell)
  d_cal2_cell      <- fun.cal(cat_twoway, d_eq_cell)
  d_nps1_rake_cell <- fun.cal(cat_oneway, d_nps1_cell)
  d_nps2_rake_cell <- fun.cal(cat_twoway, d_nps2_cell)

  ## Cell totals -> per-individual weights (divide by cell count)
  dt_agg_aou %>%
    mutate(
      d_unweighted     = nsiz / sum(dt_agg_aou$weight),
      d_cal1           = d_cal1_cell / weight,
      d_cal2           = d_cal2_cell / weight,
      d_nps1           = d_nps1_cell / weight,
      d_nps2           = d_nps2_cell / weight,
      d_nps1_rake      = d_nps1_rake_cell / weight,
      d_nps2_rake      = d_nps2_rake_cell / weight,
      d_rails          = weights_rails / weight,
      selected_terms   = paste(names_twoway_new, collapse = " + "),
      calibrated_terms = if (n_used > 0) {
        paste(names_twoway_new[seq_len(n_used)], collapse = " + ")
      } else {
        NA_character_
      }
    )
}

##### Subgroup RAILS estimator, two-way variant
### Same loop as fun.sub.rails.threeway, but each stratum's population
### margins are one-way + two-way only and the estimator is fun.rails.twoway.
fun.sub.rails.twoway <- function(
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
    stop("fun.sub.rails.twoway: subgroup_var must not appear in names_univar")
  }
  if (!subgroup_var %in% names(dt_agg_aou) || !subgroup_var %in% names(dt_agg_pums)) {
    stop("fun.sub.rails.twoway: subgroup_var must be a column of both cell tables")
  }

  twovars     <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  max_formula <- formula(paste0("~", paste(c(names_univar, twovars), collapse = "+")))

  subgroup_levels <- sort(unique(dt_agg_aou[[subgroup_var]]))

  out <- lapply(subgroup_levels, function(lev) {

    message("=== Subgroup ", subgroup_var, " = ", lev, " (two-way) ===")

    agg_aou_sub  <- dt_agg_aou  %>% filter(.data[[subgroup_var]] == lev)
    agg_pums_sub <- dt_agg_pums %>% filter(.data[[subgroup_var]] == lev)

    if (nrow(agg_pums_sub) == 0) {
      warning("fun.sub.rails.twoway: no PUMS cells for ", subgroup_var,
              " = ", lev, " — skipped.")
      return(NULL)
    }

    mat_max    <- sparse.model.matrix(max_formula, data = agg_pums_sub, keep.order = TRUE)
    pop_totals <- setNames(
      as.numeric(Matrix::crossprod(mat_max, agg_pums_sub$weight)),
      colnames(mat_max)
    )

    fun.rails.twoway(
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
