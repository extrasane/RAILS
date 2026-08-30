# Master simulation driver for the subgroup calibration study.
# Author: Huiding Chen
#
# sim.fun() runs ONE replicate: it generates a population, draws the
# probability (s) and non-probability (aou) samples, and evaluates all eight
# competing estimators, returning a named list consumed by fun.rep()
# (Results/Res_subgroup/subgroup_summary.Rmd).
#
# Requires the shared functions in subgroup_functions.R (sourced by the caller):
#   sample.function, outcome.function, fun.nps, fun.gvs, fun.lifo,
#   fun.oneway.rake, fun.subgroup.rake, fun.out, fun.pre, transform_term*.
#
# The out list keys, in order, map to the eight reported methods:
#   unweighted (Naive), trueweight (Oracle), svy_oneway_rake (G-Raking),
#   nps (NPS), GVS, rails (G-RAILS), rake_subgroup (S-Raking),
#   rails_subgroup (S-RAILS); rails_subgroup_recalibrated + avg_size +
#   zero_cells are the trailing diagnostics fun.rep drops.

sim.fun <- function(seed,nsiz, alpha, beta, gamma, pars ,names_var,
                    pre_pop, true_y, tolerance = 1e-7) {
  #### 1) Generate the Data ####
  dt <- sample.function(seed,nsiz,z_prob = pars$z_prob,p1_list = pars$p1,p2_list = pars$p2,
                        p3_list = pars$p3, p4_list = pars$p4)
  dt <- outcome.function(dt,alpha,beta,gamma)
  temp_y <- mean(dt$y)            # Sample (biased) mean
  ntaou  <- sum(dt$aou)           # AoU size
  nts    <- sum(dt$s)             # Probability-sample size

  dt_aou <- dt[dt$aou == 1, ]
  dt_s   <- dt[dt$s   == 1, ]
  y_aou  <- as.matrix(dt_aou$y)

  #### 2) Build two-way indicators and remove zero cells ####
  num_univar   <- sum(stringr::str_detect(names(dt_aou), "^catx\\d+$"))
  names_univar <- paste0("catx", 1:num_univar)

  # all 2-way combos for AoU
  twovars <- apply(t(combn(names_univar, 2)), 1, paste, collapse = ":")
  #all_aou <- model.matrix(as.formula(paste0("~0+", paste0(twovars, collapse = "+"))), dt_aou)

  # zero-cell detection via survey margins
  NA_pop <- rep(NA_real_, 16)

  clus_zerodetect <- survey::svydesign(id = ~1, weights = ~1, fpc = ~fpc, data = dt_aou)
  margin_names <- cal_names(as.formula(paste0("~", paste0(twovars, collapse = "+"))), clus_zerodetect)
  transformed_names <- sapply(margin_names, transform_term_aou)
  poptotals_zerodetect <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  index_zerocell <- names(poptotals_zerodetect)[poptotals_zerodetect == 0]

  # remove interactions involving empty levels
  modify_string <- function(s) sub("x(\\d)\\d", "x\\1", s)
  if (length(index_zerocell) > 0) {
    bad <- unique(sapply(index_zerocell, modify_string))
    twovars <- twovars[!stringr::str_detect(twovars, paste0(bad, collapse = "|"))]
  }

  # sub-matrix of selected variable(s)
  if (length(names_var) == 1) {
    names_select <- names_var
  } else {
    names_select <- apply(t(combn(names_var, 2)), 1, paste, collapse = ":")
  }
  sub_aou <- model.matrix(as.formula(paste0("~0+", paste0(names_select, collapse = "+"))), dt_aou)
  n_totallevels <- ncol(sub_aou)
  NA_sub <- matrix(NA_real_, nrow = 7, ncol = n_totallevels)

  #### Output container ####
  out <- list()
  names_pop <- c("mean","variance","CI_lower","CI_upper",
                 "bias","bias2","Coverage","TrueCoverage","size","VarW","NegW",
                 "SvyVar","SvyLB","SvyUB","SvyCov","Memory")

  #### M1: Unweighted ####
  result_unweighted <- list()
  t0 <- proc.time()
  m <- mean(dt_aou$y); v <- var(dt_aou$y)
  ci <- c(m - qnorm(0.975) * sqrt(v/ntaou), m + qnorm(0.975) * sqrt(v/ntaou))
  pop_unweighted <- c(m, v, ci,
                      m - temp_y, (m - temp_y)^2,
                      as.numeric(temp_y >= ci[1] & temp_y <= ci[2]),
                      as.numeric(true_y >= ci[1] & true_y <= ci[2]),
                      ntaou, 0, 0, 0, 0, 0, 0, 0)
  names(pop_unweighted) <- names_pop
  result_unweighted$pop  <- pop_unweighted
  result_unweighted$sub  <- fun.pre(sub_aou, rep(1, nrow(sub_aou)), dt_aou$y, pre_pop)
  result_unweighted$time <- (proc.time() - t0)["elapsed"]
  out$unweighted <- result_unweighted

  #### M2: True weights (1/paou) ####
  result_trueweight <- list(); t0 <- proc.time()
  pop_trueweight <- fun.out(y_aou, 1/dt_aou$paou, true_y, temp_y, names_pop)
  result_trueweight$pop  <- pop_trueweight
  result_trueweight$sub  <- fun.pre(sub_aou, 1/dt_aou$paou, dt_aou$y, pre_pop)
  result_trueweight$time <- (proc.time() - t0)["elapsed"]
  out$trueweight <- result_trueweight

  #### M3: One-way survey calibration (raking) ####
  result_svy_oneway_rake <- list()
  t0 <- proc.time()
  clus_oneway_cal <- fun.oneway.rake(dt,dt_aou,nsiz,ntaou,names_univar,d = nsiz/ntaou)
  if(!is.null(clus_oneway_cal)){
    pop_svy_oneway_rake <- fun.out(y_aou, weights(clus_oneway_cal), true_y, temp_y, names_pop, clus_oneway_cal)
    sub_svy_oneway_rake <- fun.pre(sub_aou, weights(clus_oneway_cal), dt_aou$y, pre_pop)
    time_svy_oneway_rake <- (proc.time() - t0)["elapsed"]
    result_svy_oneway_rake$pop  <- pop_svy_oneway_rake
    result_svy_oneway_rake$sub  <- sub_svy_oneway_rake
    result_svy_oneway_rake$time <- time_svy_oneway_rake
  }else{
    result_svy_oneway_rake$pop  <- NA_pop
    result_svy_oneway_rake$sub  <- NA_sub
    result_svy_oneway_rake$time <- (proc.time() - t0)["elapsed"]
  }
  out$svy_oneway_rake <- result_svy_oneway_rake

  #### M4: NPS (intercept-only pseudo-likelihood) ####
  result_nps <- list()
  t0 <- proc.time()
  dt_s$w <- 1/dt_s$ps
  fit0 <- try(fun.nps(as.formula(paste0("~",paste0(names_univar,collapse = "+"))),dt, dt_aou, dt_s), silent = TRUE)
  if (inherits(fit0, "try-error")) {
    out$error <- "NPS Failed"; return(out)
  }
  d <- fit0$d; LKHD <- fit0$LKHD
  result_nps$pop  <- fun.out(y_aou, d, true_y, temp_y, names_pop)
  result_nps$sub  <- fun.pre(sub_aou, d, dt_aou$y, pre_pop)
  result_nps$time <- (proc.time() - t0)["elapsed"]
  out$nps <- result_nps

  #### M5: Greedy VS via cached RAILS (two-way) ####
  t0 <- proc.time()
  vs_fit <- suppressMessages(fun.gvs(names_univar,dt ,dt_aou, dt_s, 0.05,FALSE))
  if (is.na(vs_fit[[1]][1])) {
    out$error <- "Variable Selection Error"; return(out)
  }
  d_vs <- vs_fit$d_vs
  names_univar_new <- vs_fit$names_univar_new
  result_vs <- list()
  result_vs$pop  <- fun.out(y_aou, d_vs, true_y, temp_y, names_pop)
  result_vs$sub  <- fun.pre(sub_aou, d_vs, dt_aou$y, pre_pop)
  result_vs$time <- (proc.time() - t0)["elapsed"]
  out$GVS <- result_vs

  #### M6: Stepwise RAILS + raking over selected terms ####
  result_rails <- list()
  t0 <- proc.time()
  temp_result <- suppressMessages(fun.lifo(dt,dt_aou, dt_s,names_univar_new))
  if (is.null(temp_result)) {
    out$error <- "LIFO Failed"; return(out)
  }
  clus_cal <- temp_result$clus_cal
  w_rails <- weights(clus_cal)
  n_stepwise <- temp_result$n_stepwise
  model_rails <- temp_result$model

  pop_rails <- fun.out(y_aou, w_rails, true_y, temp_y, names_pop, clus_cal)
  sub_rails <- fun.pre(sub_aou, w_rails, dt_aou$y, pre_pop)
  time_rails <- (proc.time() - t0)["elapsed"]

  result_rails$pop <- pop_rails
  result_rails$sub <- sub_rails
  result_rails$time <- time_rails
  result_rails$n_stepwise <- n_stepwise
  result_rails$model <- model_rails

  out$rails <- result_rails

  #### M7: NPS national, subgroup raking (NR) ####
  result_mix_NR <- list()
  t0 <- proc.time()
  dt_NR <- fun.subgroup.rake(names_univar,names_var,dt,dt_aou,1)
  result_mix_NR$n_nonconvergence <- sum(dt_NR$index_nonconv)
  dt_NR <- dt_NR$out
  time_NR <- (proc.time() - t0)["elapsed"]
  result_mix_NR$pop <- fun.out(dt_NR$y, dt_NR$ori, true_y, temp_y)
  result_mix_NR$sub <- fun.pre(sub_aou, dt_NR$ori, dt_NR$y, pre_pop)
  result_mix_NR$time <- time_NR
  out$rake_subgroup <- result_mix_NR

  # National re-calibration after subgroup raking
  # clus <- survey::svydesign(id = ~1, weights = ~ori, data = dt_mix_NR)
  # margin_names <- cal_names(as.formula(paste0("~", paste0(names_univar, collapse = "+"))), clus)
  # transformed_names <- sapply(margin_names, transform_term)
  # poptotals <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  # result_mix_NRC <- list()
  # try({
  #   clus_cal <- survey::calibrate(clus, as.formula(paste0("~", paste0(names_univar, collapse = "+"))),
  #                                 poptotals, calfun = "raking", maxit = 1e2,
  #                                 epsilon = rep(1e-7, length(poptotals)))
  #   w_NRC <- weights(clus_cal)
  #   result_mix_NRC$pop <- fun.out(dt_mix_NR$y, w_NRC, true_y, temp_y)
  #   result_mix_NRC$sub <- fun.pre(sub_aou, w_NRC, dt_mix_NR$y, pre_pop)
  #   result_mix_NRC$time <- time_mix_NR
  #   result_mix_NRC$n_nonconvergence <- sum(index_nonconv)
  # }, silent = TRUE)
  # out$rake_subgroup_recalibrated <- result_mix_NRC

  #### M8: NPS national, subgroup RAILS (RR) with LIFO ####
  result_mix_RR <- list(); result_mix_RRC <- list(); t0 <- proc.time()
  dt_mix_RR <- cbind(dt_aou, weights = w_rails, ori = 1)
  dt_mix_RR_sub <- eval(parse(text = paste0("split(dt_mix_RR, dt_mix_RR$", names_var, ")")))
  names_sub <- names(dt_mix_RR_sub)
  names_univar_mix <- names_univar[!stringr::str_detect(names_univar, names_var)]
  index_nonconv <- numeric(length(names_sub))
  model_RR_ori  <- vector("list", length(names_sub))
  model_RR_lifo <- vector("list", length(names_sub))

  for (j in seq_along(names_sub)) {
    temp <- dt_mix_RR_sub[[j]]
    # Select within subgroup (uses dt_s for probability sample as in global)
    temp_s <- dt_s[dt_s[[names_var]] == names_sub[j], ]
    sel <- suppressMessages(fun.gvs(names_univar_mix,dt ,temp, temp_s, 0.05,FALSE))
    if (is.na(sel[[1]][1])) { index_nonconv[j] <- 1; temp$ori <- temp$weights; dt_mix_RR_sub[[j]] <- temp; next }
    model_RR_ori[[j]] <- sel$names_univar_new

    # LIFO back-off until raking converges at subgroup level
    names_sub_sel <- sel$names_univar_new
    n_uni <- length(names_sub_sel) + 1
    clus <- survey::svydesign(id = ~1, weights = ~weights, data = temp)
    index_con <- NA
    while (is.na(index_con) && n_uni != 1) {
      n_uni <- n_uni - 1
      cat_temp <- as.formula(paste0("~", paste0(names_sub_sel[seq_len(n_uni)], collapse = "+")))
      margin_names <- cal_names(cat_temp, clus)
      transformed_names <- sapply(margin_names, transform_term_sub, univar = names_var, j = names_sub[j])
      poptotals <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
      index_con <- 1
      tryCatch({
        clus_cal <- survey::calibrate(clus, cat_temp, poptotals, calfun = "raking", maxit = 1e3,
                                      epsilon = rep(tolerance, length(poptotals)))
        temp$ori <- weights(clus_cal)
        dt_mix_RR_sub[[j]] <- temp
        model_RR_lifo[[j]] <- names_sub_sel[seq_len(n_uni)]
      }, error = function(e) index_con <<- NA)
    }
    if (n_uni == 1 && is.na(index_con)) { index_nonconv[j] <- 1; temp$ori <- temp$weights; dt_mix_RR_sub[[j]] <- temp }
  }
  time_mix_RR <- (proc.time() - t0)["elapsed"]
  dt_mix_RR <- do.call(rbind, dt_mix_RR_sub)
  result_mix_RR$pop <- fun.out(dt_mix_RR$y, dt_mix_RR$ori, true_y, temp_y)
  result_mix_RR$sub <- fun.pre(sub_aou, dt_mix_RR$ori, dt_mix_RR$y, pre_pop)
  result_mix_RR$time <- time_mix_RR
  result_mix_RR$n_nonconvergence <- sum(index_nonconv)
  result_mix_RR$model_ori <- model_RR_ori
  result_mix_RR$model_lifo <- model_RR_lifo
  result_mix_RR$model_index <- is.null(unlist(model_RR_lifo))
  out$rails_subgroup <- result_mix_RR

  # National recalibration
  clus <- survey::svydesign(id = ~1, weights = ~ori, data = dt_mix_RR)
  margin_names <- cal_names(as.formula(paste0("~", paste0(names_univar, collapse = "+"))), clus)
  transformed_names <- sapply(margin_names, transform_term)
  poptotals <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  result_mix_RRC <- list()
  try({
    clus_cal <- survey::calibrate(clus, as.formula(paste0("~", paste0(names_univar, collapse = "+"))),
                                  poptotals, calfun = "raking", maxit = 1e2,
                                  epsilon = rep(1e-7, length(poptotals)))
    w_RRC <- weights(clus_cal)
    result_mix_RRC$pop <- fun.out(dt_mix_RR$y, w_RRC, true_y, temp_y)
    result_mix_RRC$sub <- fun.pre(sub_aou, w_RRC, dt_mix_RR$y, pre_pop)
    result_mix_RRC$time <- time_mix_RR
    result_mix_RRC$n_nonconvergence <- sum(index_nonconv)
  }, silent = TRUE)
  out$rails_subgroup_recalibrated <- result_mix_RRC


  #### Diagnostics ####
  out$avg_size   <- colSums(sub_aou)
  out$zero_cells <- index_zerocell
  out
}
