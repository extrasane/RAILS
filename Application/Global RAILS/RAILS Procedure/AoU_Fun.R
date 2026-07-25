#### AoU_Fun.R — Core functions for Global RAILS
#### Run inside the AoU Researcher Workbench

library(dplyr)
library(Matrix)
library(survey)

##### Weight diagnostics
fun.out <- function(w, ...) {
  temp <- numeric(6)
  temp[1] <- sum(w, na.rm = TRUE)
  temp[2] <- var(w, na.rm = TRUE)
  temp[3] <- mean(w <= 0, na.rm = TRUE)
  temp[4] <- sum(w < 1)
  temp[5] <- min(w, na.rm = TRUE)
  temp[6] <- max(w, na.rm = TRUE)
  names(temp) <- c("sum", "var", "non-positive ratio", "Less than 1", "min", "max")
  return(temp)
}

##### Forward likelihood-ratio test for one candidate interaction term
### Returns: c(log-likelihood, LRT statistic, p-value, LRT/df, df)
fun.lkd <- function(x, names_univar, LKHD, dt_s, dt_aou, m_s) {
  cat_formula <- formula(paste0("~", paste0(paste0(names_univar, collapse = "+"), "+", x)))
  t_s   <- sparse.model.matrix(cat_formula, dt_s)
  t_aou <- sparse.model.matrix(cat_formula, dt_aou)
  theta <- rep(0, ncol(t_aou))
  pia   <- rep(1 / 2, nrow(t_s))
  W_s   <- dt_s$weight
  W_aou <- dt_aou$weight
  U1    <- as.numeric(crossprod(t_aou, W_aou))

  if (any(U1 == 0)) return(rep(NA, 5))

  res       <- 1
  temp_LKHD <- NA
  iter      <- 0

  tryCatch({
    while (res >= 1e-3) {
      iter <- iter + 1
      if (iter > 1000) stop("non-convergence")
      W1      <- W_s * pia
      Uscore  <- U1 - as.numeric(crossprod(t_s, W1))
      Hsov    <- solve(as.matrix(crossprod(t_s, t_s * (W1 * (1 - pia)))))
      theta0  <- as.numeric(theta + Hsov %*% Uscore)
      pia     <- plogis(as.numeric(t_s %*% theta0))
      res     <- sum((Hsov %*% Uscore)^2)
      theta   <- theta0
      temp_LKHD <- sum(U1 * theta) - sum(W_s * log1p(exp(as.numeric(t_s %*% theta))))
    }
    df   <- ncol(t_s) - ncol(m_s)
    LRT  <- 2 * (temp_LKHD - LKHD)
    pval <- 1 - pchisq(LRT, df = df)
    return(c(temp_LKHD, LRT, pval, LRT / df, df))
  }, error = function(e) rep(NA, 5))
}

##### Pseudo-likelihood propensity score estimation (Newton-Raphson)
### Returns NPS weights scaled to nsiz; theta_init is an optional warm start
### (same unique maximum, fewer iterations).
fun.nps <- function(cat_temp, dt_aou, dt_s, nsiz = sum(dt_s$weight), theta_init = NULL) {
  m_aou  <- sparse.model.matrix(cat_temp, dt_aou)
  m_nhis <- sparse.model.matrix(cat_temp, dt_s)
  theta  <- if (is.null(theta_init)) rep(0, ncol(m_aou)) else
    c(theta_init, rep(0, ncol(m_aou) - length(theta_init)))
  pia    <- plogis(as.numeric(m_nhis %*% theta))
  W_s    <- dt_s$weight
  W_aou  <- dt_aou$weight
  U1     <- as.numeric(crossprod(m_aou, W_aou))

  res  <- 1
  iter <- 0
  while (res >= 1e-10) {
    iter <- iter + 1
    if (iter > 1000) stop("fun.nps: no convergence within 1000 iterations")
    W1     <- W_s * pia
    Uscore <- U1 - as.numeric(crossprod(m_nhis, W1))
    Hsov   <- solve(as.matrix(crossprod(m_nhis, m_nhis * (W1 * (1 - pia)))))
    theta0 <- as.numeric(theta + Hsov %*% Uscore)
    pia    <- plogis(as.numeric(m_nhis %*% theta0))
    res    <- sum((Hsov %*% Uscore)^2)
    theta  <- theta0
  }

  d <- 1 / plogis(as.numeric(m_aou %*% theta))
  d <- d / sum(d) * nsiz
  attr(d, "theta") <- theta
  return(d)
}

##### Helper: standardize interaction term name (sort components alphabetically)
standardize_string <- function(x) {
  paste(sort(unlist(strsplit(x, ":"))), collapse = ":")
}

##### Helper: match calibration targets to design columns (returned in v2 order)
create_v3 <- function(v1, v2) {
  v1_std <- vapply(names(v1), standardize_string, character(1), USE.NAMES = FALSE)
  v2_std <- vapply(v2,        standardize_string, character(1), USE.NAMES = FALSE)
  idx  <- match(v2_std, v1_std)
  keep <- !is.na(idx)
  setNames(as.numeric(v1[idx[keep]]), v2[keep])
}

##### Global RAILS estimator
### 1. Two-way NPS model for the baseline likelihood.
### 2. Forward three-way LRT selection (fun.lkd).
### 3. LIFO stepwise raking: add selected terms one at a time, stop at the
###    first non-convergent step, keep the last successful model's weights.
### The model space is specifically the two-way base + selected THREE-WAY
### interactions — hence the name. All computation is on aggregated cells;
### output weight columns are PER-INDIVIDUAL (cell totals / cell count).
### Besides RAILS, returns the benchmark methods:
###   d_unweighted              — equal weights nsiz / n (naive benchmark)
###   d_cal1 / d_cal2           — raking from equal weights, one-way / two-way margins
###   d_nps1 / d_nps2           — NPS weights, one-way / two-way model (no raking)
###   d_nps1_rake / d_nps2_rake — NPS starting weights raked to the same margins
###   d_rails                   — RAILS weights (stepwise raked, final model)
fun.rails.threeway <- function(
    dt_agg_aou,
    dt_agg_pums,
    pop_totals,
    names_univar = c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
    alpha        = 0.05,
    nsiz         = sum(dt_agg_pums$weight)
) {
  ## Two-way base model and its likelihood
  twovars     <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  names_twoway <- c(names_univar, twovars)
  cat_twoway  <- formula(paste0("~", paste(names_twoway, collapse = "+")))

  d_init <- fun.nps(cat_twoway, dt_agg_aou, dt_agg_pums, nsiz = nsiz)
  m_pums <- sparse.model.matrix(cat_twoway, dt_agg_pums)
  m_aou  <- sparse.model.matrix(cat_twoway, dt_agg_aou)
  theta2 <- attr(d_init, "theta")
  LKHD2  <- sum(as.numeric(crossprod(m_aou, dt_agg_aou$weight)) * theta2) -
    sum(dt_agg_pums$weight * log1p(exp(as.numeric(m_pums %*% theta2))))

  ## Forward three-way selection
  vars_new           <- combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))
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
      m_s_new <- sparse.model.matrix(
        as.formula(paste0("~", paste(names_threeway_new, collapse = "+"))),
        dt_agg_pums
      )
      vars_new <- vars_new[vars_new != selected_var]
    }
  }

  cat_formula <- formula(paste0("~", paste(names_threeway_new, collapse = "+")))

  ## Subset pop_totals ONCE to the full selected model's margins;
  ## each LIFO step then subsets this small vector.
  pop_totals <- create_v3(pop_totals, colnames(sparse.model.matrix(cat_formula, dt_agg_aou)))

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
    cat_temp <- formula(paste0("~", paste(names_threeway_new[seq_len(n_uni)], collapse = "+")))
    temp_d   <- fun.nps(cat_temp, dt_agg_aou, dt_agg_pums, nsiz = nsiz, theta_init = theta_prev)
    theta_prev <- attr(temp_d, "theta")

    temp <- dt_agg_aou %>%
      mutate(count = temp_d * weight, count = count / sum(count) * nsiz)

    clus     <- svydesign(id = ~1, weights = ~count, data = temp)
    names_ps <- cal_names(cat_temp, clus)
    vars_cal <- create_v3(pop_totals, names_ps)

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
        "fun.rails.threeway: covariate combination '", names_threeway_new[n_uni],
        "' fails to pile up to higher order — keeping the previous model."
      )
      break
    }
  }

  ## -------------------------------------------------------------------
  ## Benchmark methods (all on the same aggregated cells):
  ##   cal-1 / cal-2:           raking from equal starting weights to the
  ##                            one-way / one-way + two-way PUMS margins
  ##   nps-1 / nps-2:           pseudo-likelihood NPS weights, one-way /
  ##                            two-way model (no raking)
  ##   nps-1-rake / nps-2-rake: the corresponding NPS weights as starting
  ##                            point, raked to the same margins
  ## -------------------------------------------------------------------
  cat_oneway <- formula(paste0("~", paste(names_univar, collapse = "+")))

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

  ## Rake a cell design with the given starting cell totals to the margins
  ## of cat_f; NA (with warning) on non-convergence.
  fun.cal <- function(cat_f, start_cell) {
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
      d_cal1           = d_cal1_cell / weight,
      d_cal2           = d_cal2_cell / weight,
      d_nps1           = d_nps1_cell / weight,
      d_nps2           = d_nps2_cell / weight,
      d_nps1_rake      = d_nps1_rake_cell / weight,
      d_nps2_rake      = d_nps2_rake_cell / weight,
      d_rails          = weights_rails / weight,
      selected_terms   = paste(names_threeway_new, collapse = " + "),
      calibrated_terms = if (n_used > 0) {
        paste(names_threeway_new[seq_len(n_used)], collapse = " + ")
      } else {
        NA_character_
      }
    )
}
