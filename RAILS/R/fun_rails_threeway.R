#' Global RAILS estimator (two-way base + selected three-way interactions)
#'
#' Applies the full RAILS procedure to the entire population:
#' \enumerate{
#'   \item Fits a two-way non-probability propensity model to obtain a
#'     baseline log-likelihood.
#'   \item Forward likelihood-ratio selection of three-way interactions
#'     (via [fun.lkd()]) at significance level `alpha`.
#'   \item LIFO stepwise raking: the selected three-way terms are added one at
#'     a time in selection order, re-fitting the propensity model and raking to
#'     the population margins at each step; the walk stops at the first
#'     non-convergent step and keeps the last successful model's weights.
#'   \item Computes benchmark weighting methods on the same cells.
#' }
#' All computation is on aggregated covariate cells. Population margins are
#' passed in via `pop_totals` (precomputed once by the caller) rather than
#' recomputed internally.
#'
#' @param dt_agg_aou Aggregated non-probability (AoU) cell counts, with a
#'   column `weight` giving the cell size.
#' @param dt_agg_pums Aggregated reference-sample (PUMS) cell counts, with a
#'   column `weight`.
#' @param pop_totals Named numeric vector of population margins (one-way,
#'   two-way, and three-way of `names_univar`), e.g. from
#'   `Matrix::crossprod()` of a sparse model matrix on `dt_agg_pums`.
#' @param names_univar Character vector of main-effect variables. All two-way
#'   interactions form the base model and all three-way interactions are
#'   selection candidates.
#' @param alpha Significance threshold for forward LRT selection.
#' @param nsiz Target population total the weights are scaled to. Defaults to
#'   the total reference-sample weight; override when the reference sample
#'   differs from the source of `pop_totals` or when targeting a subgroup.
#'
#' @return `dt_agg_aou` with appended per-individual weight columns (cell
#'   totals divided by the cell count):
#'   `d_unweighted`, `d_cal1`, `d_cal2`, `d_nps1`, `d_nps2`, `d_nps1_rake`,
#'   `d_nps2_rake`, `d_rails`, plus `selected_terms` (all LRT-selected terms)
#'   and `calibrated_terms` (terms of the final converged model; `NA` if no
#'   step converged).
#'
#' @details Progress is reported via [message()]. A [warning()] names any
#'   covariate combination that fails to converge during raking; `d_cal2` may
#'   be `NA` (with a warning) when full two-way raking does not converge.
#'
#' @seealso [fun.sub.rails.threeway()] for the stratified wrapper.
#'
#' @export
fun.rails.threeway <- function(
    dt_agg_aou,
    dt_agg_pums,
    pop_totals,
    names_univar = c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
    alpha        = 0.05,
    nsiz         = sum(dt_agg_pums$weight)
) {
  ## Two-way base model and its likelihood
  twovars      <- utils::combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  names_twoway <- c(names_univar, twovars)
  cat_twoway   <- stats::formula(paste0("~", paste(names_twoway, collapse = "+")))

  d_init <- fun.nps(cat_twoway, dt_agg_aou, dt_agg_pums, nsiz = nsiz)
  m_pums <- Matrix::sparse.model.matrix(cat_twoway, dt_agg_pums)
  m_aou  <- Matrix::sparse.model.matrix(cat_twoway, dt_agg_aou)
  theta2 <- attr(d_init, "theta")
  LKHD2  <- sum(as.numeric(Matrix::crossprod(m_aou, dt_agg_aou$weight)) * theta2) -
    sum(dt_agg_pums$weight * log1p(exp(as.numeric(m_pums %*% theta2))))

  ## Forward three-way selection
  vars_new           <- utils::combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))
  names_threeway_new <- names_twoway
  LKHD_new           <- LKHD2
  m_s_new            <- m_pums

  index_add <- 1
  while (index_add != 0) {
    if (length(vars_new) == 0) break
    temp_p <- apply(
      as.data.frame(vars_new), 1, fun.lkd,
      names_univar = names_threeway_new, LKHD = LKHD_new,
      dt_s = dt_agg_pums, dt_aou = dt_agg_aou, m_s = m_s_new
    )
    colnames(temp_p) <- vars_new
    temp_sign <- temp_p[, temp_p[3, ] < alpha, drop = FALSE]
    if (ncol(temp_sign) == 0 || all(is.na(temp_sign))) {
      index_add <- 0
    } else {
      selected_var       <- colnames(temp_sign)[which.max(temp_sign[4, ])]
      message(format(Sys.time(), "%H:%M:%S"), "  GVS: selected '", selected_var, "'")
      names_threeway_new <- c(names_threeway_new, selected_var)
      LKHD_new           <- temp_sign[1, selected_var]
      m_s_new <- Matrix::sparse.model.matrix(
        stats::as.formula(paste0("~", paste(names_threeway_new, collapse = "+"))),
        dt_agg_pums
      )
      vars_new <- vars_new[vars_new != selected_var]
    }
  }

  cat_formula <- stats::formula(paste0("~", paste(names_threeway_new, collapse = "+")))

  ## Subset pop_totals ONCE to the full selected model's margins;
  ## each LIFO step then subsets this small vector.
  pop_totals <- create_v3(pop_totals, colnames(Matrix::sparse.model.matrix(cat_formula, dt_agg_aou)))

  ## LIFO stepwise raking — ascending, stop at first failure
  n_base <- length(names_twoway)
  n_max  <- length(names_threeway_new)

  weights_rails <- rep(NA_real_, nrow(dt_agg_aou))
  n_used        <- 0
  theta_prev    <- theta2

  if (n_max == n_base) warning("fun.rails.threeway: no three-way term selected — RAILS weights are NA.")

  n_uni <- n_base
  while (n_uni < n_max) {
    n_uni    <- n_uni + 1
    message(format(Sys.time(), "%H:%M:%S"), "  LIFO step ", n_uni - n_base, "/",
            n_max - n_base, ": adding '", names_threeway_new[n_uni], "'")
    cat_temp <- stats::formula(paste0("~", paste(names_threeway_new[seq_len(n_uni)], collapse = "+")))
    temp_d   <- fun.nps(cat_temp, dt_agg_aou, dt_agg_pums, nsiz = nsiz, theta_init = theta_prev)
    theta_prev <- attr(temp_d, "theta")

    temp <- dt_agg_aou %>%
      mutate(count = temp_d * .data$weight, count = .data$count / sum(.data$count) * nsiz)

    clus     <- survey::svydesign(id = ~1, weights = ~count, data = temp)
    names_ps <- survey::cal_names(cat_temp, clus)
    vars_cal <- create_v3(pop_totals, names_ps)

    step_ok <- tryCatch({
      clus_cal <- survey::calibrate(
        clus, cat_temp, vars_cal,
        calfun = "raking", maxit = 1e3,
        epsilon = rep(1e-7, length(vars_cal))
      )
      weights_rails <- stats::weights(clus_cal)
      n_used        <- n_uni
      TRUE
    }, error = function(e) FALSE)

    if (!step_ok) {
      warning(
        "fun.rails.threeway: covariate combination '", names_threeway_new[n_uni],
        "' fails to pile up to higher order — keeping the previous model."
      )
      break
    }
  }

  ## Benchmark methods (all on the same aggregated cells)
  cat_oneway <- stats::formula(paste0("~", paste(names_univar, collapse = "+")))

  d_nps1_cell <- {
    d1 <- fun.nps(cat_oneway, dt_agg_aou, dt_agg_pums, nsiz = nsiz)
    ct <- d1 * dt_agg_aou$weight
    ct / sum(ct) * nsiz
  }
  d_nps2_cell <- {
    ct <- d_init * dt_agg_aou$weight    # two-way NPS already fitted above
    ct / sum(ct) * nsiz
  }
  d_eq_cell <- dt_agg_aou$weight / sum(dt_agg_aou$weight) * nsiz

  ## Rake a cell design with the given starting cell totals to the margins of
  ## cat_f; NA (with warning) on non-convergence.
  fun.cal <- function(cat_f, start_cell) {
    temp_c <- dt_agg_aou %>% mutate(count = start_cell)
    clus_c <- survey::svydesign(id = ~1, weights = ~count, data = temp_c)
    vars_c <- create_v3(pop_totals, survey::cal_names(cat_f, clus_c))
    tryCatch(
      stats::weights(survey::calibrate(
        clus_c, cat_f, vars_c,
        calfun = "raking", maxit = 1e3,
        epsilon = rep(1e-7, length(vars_c))
      )),
      error = function(e) {
        warning("fun.rails.threeway: benchmark raking failed to converge (",
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
      d_cal1           = d_cal1_cell / .data$weight,
      d_cal2           = d_cal2_cell / .data$weight,
      d_nps1           = d_nps1_cell / .data$weight,
      d_nps2           = d_nps2_cell / .data$weight,
      d_nps1_rake      = d_nps1_rake_cell / .data$weight,
      d_nps2_rake      = d_nps2_rake_cell / .data$weight,
      d_rails          = weights_rails / .data$weight,
      selected_terms   = paste(names_threeway_new, collapse = " + "),
      calibrated_terms = if (n_used > 0) {
        paste(names_threeway_new[seq_len(n_used)], collapse = " + ")
      } else {
        NA_character_
      }
    )
}
