##### Simulation Function v13

### Link function:
### Logit for initial PS and logit for forward selection
### Add QR for twoway
### Add Estimating Equation residuals for each weight
### Add Comparison wrt svymean variance

### Report Subgroup Prevalence & Variance
### Simulation Function

sim.fun <- function(seed,nsiz,alpha,beta,gamma,names_var,true_model,...){
  set.seed(seed)
  
  ##### Variable Specification
  dt <- sample.function(nsiz,alpha,beta,gamma)
  # Mean for current simulation, biased
  temp_y <- mean(dt$y)
  
  # Size of AoU
  ntaou <- sum(dt$aou)
  # Size of Probability Samples
  nts <- sum(dt$s)
  dt_aou <- dt[dt$aou==1, ] 
  dt_s <- dt[dt$s==1, ]     
  y_aou <- as.matrix(dt_aou$y)
  
  ##### All two-way combinations #####
  num_univar <- sum(str_detect(names(dt_aou),"catx"))
  names_univar <- paste0("catx",1:num_univar)
  temp <- t(combn(names_univar,2))
  twovars <- temp %>% apply(1,paste,collapse = ":")
  all_aou <- model.matrix(as.formula(paste0("~0+",paste0(twovars,collapse = "+"))), dt_aou)
  NA_pop <- rep(NA,15)
  
  ### Rule out if there's any combination has zero cell
  clus_zerodetect <- svydesign(id=~1, weights = ~1, fpc =~fpc, data= dt_aou)
  margin_names <- cal_names(as.formula(paste0("~",paste0(twovars,collapse = "+"))),clus_zerodetect)
  transformed_names <- sapply(margin_names, transform_term_aou)
  poptotals_zerodetect <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  # Detect Zero Cells
  index_zerocell <- poptotals_zerodetect == 0
  index_zerocell <- names(poptotals_zerodetect)[index_zerocell]
  
  modify_string <- function(s) {
    gsub("x(\\d)\\d", "x\\1", s)
  }
  
  if(length(index_zerocell)>0){
    index_zerocell <- unique(sapply(index_zerocell, modify_string))
    twovars <- twovars[!str_detect(twovars,paste0(index_zerocell,collapse = "|"))]
  }
  
  ##### All two-way combinations for selected variables #####
  temp <- t(combn(names_var,2))
  names_select <- temp %>% apply(1,paste,collapse = ":")
  sub_aou <- model.matrix(as.formula(paste0("~0+",paste0(names_select,collapse = "+"))), dt_aou)
  n_totallevels <- dim(sub_aou)[2]
  
  
  ##### Set Output as the return list
  out <- list()
  names_pop <- c("mean","variance","CI_lower","CI_upper",
                 "bias","bias2","NomCoverage","size","VarW","NegW",
                 "SvyVar","SvyLB","SvyUB","SvyCov","Memory")
  
  ##### M1: Unweighted Result #####
  result_unweighted <- list()
  start_time <- proc.time()
  ### Population Estimation
  pop_unweighted <- c(mean(dt_aou$y),var(dt_aou$y))
  pop_unweighted <- c(pop_unweighted,
                      pop_unweighted[1] - qnorm(0.975) * sqrt(pop_unweighted[2]/ntaou),
                      pop_unweighted[1] + qnorm(0.975) * sqrt(pop_unweighted[2]/ntaou)
  )
  pop_unweighted <- c(pop_unweighted,
                      pop_unweighted[1] - temp_y,
                      (pop_unweighted[1] - temp_y)^2,
                      temp_y >= pop_unweighted[3] & temp_y <= pop_unweighted[4],
                      ntaou,0,0,0,0,0,0,0
  )
  names(pop_unweighted) <- names_pop
  result_unweighted$pop <- pop_unweighted                
  
  end_time <- proc.time()
  result_unweighted$time <- (end_time - start_time)["elapsed"]
  out[["unweighted"]] <- result_unweighted
  
  ##### M2: True Weighted Result #####
  result_trueweight <- list()
  start_time <- proc.time()
  ### Population Estimation
  pop_trueweight <- fun.out(y_aou,1/dt_aou$paou,temp_y,names_pop)
  result_trueweight$pop <- pop_trueweight
  
  
  end_time <- proc.time()
  result_trueweight$time <- (end_time - start_time)["elapsed"]
  out[["trueweight"]] <- result_trueweight
  
  ##### M3: One-way, Survey Calibration, linear & raking #####
  result_svy_oneway_linear <- list()
  result_svy_oneway_rake <- list()
  
  start_time <- proc.time()
  ### Weight Calculation
  dt_oneway <- cbind(dt_aou, weights = nsiz/ntaou)
  clus_oneway <- svydesign(id=~1, weights = ~weights, fpc =~fpc, data= dt_oneway)
  margin_names <- cal_names(as.formula(paste0("~",paste0(names_univar,collapse = "+"))),clus_oneway)
  transformed_names <- sapply(margin_names, transform_term)
  poptotals_oneway <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  
  ### Linear
  pop_svy_oneway_linear <- NA_pop
  time_svy_oneway_linear <- NA
  tryCatch(
    {
      clus_oneway_calirated <- survey::calibrate(clus_oneway, 
                                                 as.formula(paste0("~",paste0(names_univar,collapse = "+"))), 
                                                 poptotals_oneway,
                                                 calfun = "linear",maxit = 1e2) 
      pop_svy_oneway_linear <- fun.out(y_aou,weights(clus_oneway_calirated),
                                       temp_y,names_pop,
                                       clus_oneway_calirated)

      end_time <- proc.time()
      time_svy_oneway_linear <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_svy_oneway_linear$pop <- pop_svy_oneway_linear
  result_svy_oneway_linear$time <- time_svy_oneway_linear
  out[["svy_oneway_linear"]] <- result_svy_oneway_linear
  
  ### Raking
  start_time <- proc.time()
  pop_svy_oneway_rake <- NA_pop
  time_svy_oneway_rake <- NA
  tryCatch(
    {
      clus_oneway_calirated_rake <- survey::calibrate(clus_oneway, 
                                                      as.formula(paste0("~",paste0(names_univar,collapse = "+"))),
                                                      poptotals_oneway,
                                                      calfun = "raking",maxit = 1e2,
                                                      epsilon = rep(1e-7,length(poptotals_oneway))) 
      pop_svy_oneway_rake <- fun.out(y_aou,
                                     weights(clus_oneway_calirated_rake),
                                     temp_y,names_pop,
                                     clus_oneway_calirated_rake)
      end_time <- proc.time()
      time_svy_oneway_rake <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_svy_oneway_rake$pop <- pop_svy_oneway_rake
  result_svy_oneway_rake$time <- time_svy_oneway_rake
  out[["svy_oneway_rake"]] <- result_svy_oneway_rake
  
  ##### M4: Two-way (Univaraite + Two-way), Survey Calibration, linear & raking #####
  result_svy_twoway_linear <- list()
  result_svy_twoway_rake <- list()
  
  start_time <- proc.time()
  ### Weight Calculation
  ##### Calibrate on oneway + twoway
  names_twoway <- c(names_univar,twovars)
  dt_svy_twoway <- cbind(dt_aou, weights = nsiz/ntaou)
  clus_svy_twoway <- svydesign(id=~1, weights = ~weights, fpc =~fpc, data= dt_svy_twoway)
  margin_names <- cal_names(as.formula(paste0("~",paste0(names_twoway,collapse = "+"))),clus_svy_twoway)
  transformed_names <- sapply(margin_names, transform_term)
  poptotals_twoway <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  
  ### Linear
  pop_svy_twoway_linear <- NA_pop
  time_svy_twoway_linear <- NA
  
  # Remove Coliniearity
  formula_twoway <- as.formula(paste0("~",paste0(names_twoway,collapse = "+")))
  temp <- model.matrix(formula_twoway, dt_svy_twoway)
  qr_temp <- qr(temp)
  if(qr_temp$rank < ncol(temp)){
    tol <- 1e-7
    diag_R <- diag(qr.R(qr_temp))
    rank_temp <- sum(abs(diag_R) > tol)
    temp <- temp[,1:rank_temp]
    temp_names <- names(poptotals_twoway[(rank_temp+1):length(poptotals_twoway)])
    temp_names <- sapply(temp_names,  function(term) gsub("(catx)([0-9])[0-9]", "\\1\\2", term))
    temp_twowaynames <- names_twoway[!names_twoway%in%temp_names]
    temp_poptotals_twoway <- sapply(names(poptotals_twoway),
                                    function(term) gsub("(catx)([0-9])[0-9]", "\\1\\2", term))
    temp_poptotals_twoway <- poptotals_twoway[!temp_poptotals_twoway%in%temp_names]
  }
  
  tryCatch(
    {
      clus_svy_twoway_calirated <- survey::calibrate(clus_svy_twoway, 
                                                     as.formula(paste0("~",paste0(temp_twowaynames,collapse = "+"))), 
                                                     temp_poptotals_twoway,
                                                     calfun = "linear",maxit = 1e2) 
      pop_svy_twoway_linear <- fun.out(y_aou,weights(clus_svy_twoway_calirated),
                                       temp_y,
                                       names_pop,
                                       clus_svy_twoway_calirated)
      end_time <- proc.time()
      time_svy_twoway_linear <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_svy_twoway_linear$pop <- pop_svy_twoway_linear
  result_svy_twoway_linear$time <- time_svy_twoway_linear
  out[["svy_twoway_linear"]] <- result_svy_twoway_linear
  
  ### Raking
  start_time <- proc.time()
  pop_svy_twoway_rake <- NA_pop
  time_svy_twoway_rake <- NA
  tryCatch(
    {
      clus_svy_twoway_calirated_rake <- survey::calibrate(clus_svy_twoway, 
                                                          as.formula(paste0("~",paste0(temp_twowaynames,collapse = "+"))),
                                                          temp_poptotals_twoway,
                                                          calfun = "raking",maxit = 1e2,
                                                          epsilon = rep(1e-7,length(temp_poptotals_twoway))) 
      pop_svy_twoway_rake <- fun.out(y_aou,
                                     weights(clus_svy_twoway_calirated_rake),
                                     temp_y,names_pop,
                                     clus_svy_twoway_calirated_rake)
      end_time <- proc.time()
      time_svy_twoway_rake <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_svy_twoway_rake$pop <- pop_svy_twoway_rake
  result_svy_twoway_rake$time <- time_svy_twoway_rake
  out[["svy_twoway_rake"]] <- result_svy_twoway_rake
  
  ##### M5: Two-way (Univaraite + Two-way), Survey Calibration, use one-way raking if not convergent #####
  result_svy_twoway_rake_hybrid <- list()
  if(any(is.na(result_svy_twoway_rake$pop))){
    result_svy_twoway_rake_hybrid <- result_svy_oneway_rake
  }else{
    result_svy_twoway_rake_hybrid <- result_svy_twoway_rake
  }
  out[["svy_twoway_rake_hybrid"]] <- result_svy_twoway_rake_hybrid
  
  ##### M6: One-way, Pseudo-likelihood-based Propensity Score Weighting; raw, linear and raking #####
  result_ps_oneway_raw <- list()
  start_time <- proc.time()
  ##### Initial PS without variable selection
  dt_s[,"w"] <- 1/dt_s$ps
  m_s <- model.matrix(as.formula(paste0("~",paste0(names_univar,collapse = "+"))), dt_s)
  temp <- fun.vs(as.formula(paste0("~",paste0(names_univar,collapse = "+"))),dt_aou,dt_s)
  d <- temp$d
  LKHD <- temp$LKHD
  pop_ps_oneway_raw <- NA_pop
  time_ps_oneway_raw <- NA
  tryCatch(
    {
      pop_ps_oneway_raw <- fun.out(y_aou,d,temp_y,names_pop)
      end_time <- proc.time()
      time_ps_oneway_raw <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_ps_oneway_raw$pop <- pop_ps_oneway_raw
  result_ps_oneway_raw$time <- time_ps_oneway_raw
  out[["ps_oneway_raw"]] <- result_ps_oneway_raw
  
  ### Linear
  result_ps_oneway_linear <- list()
  start_time <- proc.time()
  pop_ps_oneway_linear <- NA_pop
  time_ps_oneway_linear <- NA
  dt_oneway_ps <- cbind(dt_aou, weights = d)
  clus_oneway_ps <- svydesign(id=~1, weights = ~weights, fpc =~fpc, data= dt_oneway_ps)
  
  tryCatch(
    {
      clus_ps_oneway_linear <- survey::calibrate(clus_oneway_ps, 
                                                 as.formula(paste0("~",paste0(names_univar,collapse = "+"))), 
                                                 poptotals_oneway,
                                                 calfun = "linear",maxit = 1e2) 
      pop_ps_oneway_linear <- fun.out(y_aou,weights(clus_ps_oneway_linear),
                                      
                                      temp_y,
                                      names_pop,
                                      clus_ps_oneway_linear)
      end_time <- proc.time()
      time_ps_oneway_linear <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_ps_oneway_linear$pop <- pop_ps_oneway_linear
  result_ps_oneway_linear$time <- time_ps_oneway_raw + time_ps_oneway_linear
  out[["ps_oneway_linear"]] <- result_ps_oneway_linear
  
  ### Raking
  result_ps_oneway_rake <- list()
  start_time <- proc.time()
  pop_ps_oneway_rake <- NA_pop
  time_ps_oneway_rake <- NA
  tryCatch(
    {
      clus_ps_oneway_rake <- survey::calibrate(clus_oneway_ps, 
                                               as.formula(paste0("~",paste0(names_univar,collapse = "+"))),
                                               poptotals_oneway,
                                               calfun = "raking",maxit = 1e2,
                                               epsilon = rep(1e-7,length(poptotals_oneway))) 
      pop_ps_oneway_rake <- fun.out(y_aou,weights(clus_ps_oneway_rake),temp_y,names_pop,clus_ps_oneway_rake)
      end_time <- proc.time()
      time_ps_oneway_rake <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_ps_oneway_rake$pop <- pop_ps_oneway_rake
  result_ps_oneway_rake$time <- time_ps_oneway_raw + time_ps_oneway_rake
  out[["ps_oneway_rake"]] <- result_ps_oneway_rake
  
  ##### M7: Two-way, Pseudo-likelihood-based Propensity Score Weighting; raw, linear and raking #####
  pop_ps_twoway_raw <- NA_pop
  pop_ps_twoway_linear <- NA_pop
  pop_ps_twoway_rake <- NA_pop
  time_ps_twoway_raw <- NA
  time_ps_twoway_linear <- NA
  time_ps_twoway_rake <- NA
  result_ps_twoway_raw <- list()
  result_ps_twoway_linear <- list()
  result_ps_twoway_rake <- list()
  
  start_time <- proc.time()
  m_s <- model.matrix(as.formula(paste0("~",paste0(names_twoway,collapse = "+"))), dt_s)
  d2 <- rep(NA,nrow(dt_aou))
  tryCatch({
    temp <- fun.vs(as.formula(paste0("~",paste0(names_twoway,collapse = "+"))),dt_aou,dt_s)
    d2 <- temp$d # PS two-way 
  },
  error = function(e) e
  )
  
  
  ### Return NA if non-convergent
  index_ps_twoway <- all(is.na(d2)) # If 1, indicates non-convergent, VS needed
  if(index_ps_twoway){
    result_ps_twoway_raw$pop <- pop_ps_twoway_raw
    result_ps_twoway_raw$time <- time_ps_twoway_raw
    out[["ps_twoway_raw"]] <- result_ps_twoway_raw
    result_ps_twoway_linear$pop <- pop_ps_twoway_linear
    result_ps_twoway_linear$time <- time_ps_twoway_linear
    out[["ps_twoway_linear"]] <- result_ps_twoway_linear
    result_ps_twoway_rake$pop <- pop_ps_twoway_rake
    result_ps_twoway_rake$time <- time_ps_twoway_rake
    out[["ps_twoway_rake"]] <- result_ps_twoway_rake
  } else{
    ### Raw
    pop_ps_twoway_raw <- fun.out(y_aou,d2,temp_y,names_pop)
    end_time <- proc.time()
    time_ps_twoway_raw <- (end_time - start_time)["elapsed"]
    result_ps_twoway_raw$pop <- pop_ps_twoway_raw
    result_ps_twoway_raw$time <- time_ps_twoway_raw
    out[["ps_twoway_raw"]] <- result_ps_twoway_raw
    
    ### Linear
    start_time <- proc.time()
    dt_twoway_ps <- cbind(dt_aou, weights = d2)
    clus_twoway_ps <- svydesign(id=~1, weights = ~weights, fpc =~fpc, data= dt_twoway_ps)
    tryCatch(
      {
        clus_ps_twoway_linear <- survey::calibrate(clus_twoway_ps, 
                                                   as.formula(paste0("~",paste0(names_twoway,collapse = "+"))), 
                                                   poptotals_twoway,
                                                   calfun = "linear",maxit = 1e2) 
        pop_ps_twoway_linear <- fun.out(y_aou,weights(clus_ps_twoway_linear),temp_y,names_pop,clus_ps_twoway_linear)
        end_time <- proc.time()
        time_ps_twoway_linear <- (end_time - start_time)["elapsed"]
      },
      error = function(e) e
    )
    result_ps_twoway_linear$pop <- pop_ps_twoway_linear
    result_ps_twoway_linear$time <- time_ps_twoway_raw + time_ps_twoway_linear 
    out[["ps_twoway_linear"]] <- result_ps_twoway_linear
    
    ### Raking
    start_time <- proc.time()
    tryCatch(
      {
        clus_ps_twoway_rake <- survey::calibrate(clus_twoway_ps, 
                                                 as.formula(paste0("~",paste0(names_twoway,collapse = "+"))),
                                                 poptotals_twoway,
                                                 calfun = "raking",maxit = 1e2,
                                                 epsilon = rep(1e-7,length(poptotals_twoway))) 
        pop_ps_twoway_rake <- fun.out(y_aou,weights(clus_ps_twoway_rake),temp_y,names_pop,clus_ps_twoway_rake)
        end_time <- proc.time()
        time_ps_twoway_rake <- (end_time - start_time)["elapsed"]
      },
      error = function(e) e
    )
    result_ps_twoway_rake$pop <- pop_ps_twoway_rake
    result_ps_twoway_rake$time <- time_ps_twoway_raw + time_ps_twoway_rake
    out[["ps_twoway_rake"]] <- result_ps_twoway_rake
  }
  
  
  ##### M8: Two-way, Pseudo-likelihood-based Propensity Score Weighting; Variable Selection, raw, linear and raking #####
  ### Variable Selection
  start_time <- proc.time()
  twovars_new <- twovars
  names_univar_new <- names_univar
  LKHD_new <- LKHD
  m_aou <- model.matrix(as.formula(paste0("~",paste0(names_univar,collapse = "+"))) , dt_aou)
  m_s <- model.matrix(as.formula(paste0("~",paste0(names_univar,collapse = "+"))), dt_s)
  m_s_new <- m_s
  index_add <- 1
  
  while(index_add !=0){
    temp_p <- apply(as.data.frame(twovars_new),1,fun.lkd,
                    names_univar = names_univar_new,LKHD = LKHD_new,
                    dt_s = dt_s,dt_aou = dt_aou, m_s = m_s_new)
    colnames(temp_p) <- twovars_new
    temp <- temp_p 
    
    temp_sign <- temp[,temp[3,] <0.05,drop = F]
    if(length(temp_sign)==0|all(is.na(temp_sign))){
      index_add <- 0
    }else{
      index_add <- which.max(temp_sign[4,])
      names_univar_new <- c(names_univar_new,names(index_add))
      LKHD_new <- temp_sign[1,index_add]
      m_s_new <- model.matrix(as.formula(paste0("~",paste0(names_univar_new,collapse = "+"))), dt_s)
      twovars_new <- twovars_new[!names(index_add)==twovars_new]
    }
  }
  temp <- fun.vs(as.formula(paste0("~",paste0(names_univar_new,collapse = "+"))),dt_aou,dt_s)
  d_vs <- temp$d # PS two-way
  
  ### Raw
  result_ps_vs_twoway_raw <- list()
  pop_ps_vs_twoway_raw <- NA_pop
  time_ps_vs_twoway_raw <- NA
  n_vs <- rep(NA, 2)
  tryCatch(
    {
      pop_ps_vs_twoway_raw <- fun.out(y_aou,d_vs,temp_y,names_pop)
      end_time <- proc.time()
      time_ps_vs_twoway_raw <- (end_time - start_time)["elapsed"]
      n_vs <- c(unselect = length(names_twoway) - length(names_univar_new),
                select = length(names_univar_new) - length(names_univar)
      )
    },
    error = function(e) e
  )
  result_ps_vs_twoway_raw$pop <- pop_ps_vs_twoway_raw
  result_ps_vs_twoway_raw$time <- time_ps_vs_twoway_raw
  result_ps_vs_twoway_raw$n_vs <- n_vs
  out[["ps_vs_twoway_raw"]] <- result_ps_vs_twoway_raw
  
  ### Linear
  result_ps_vs_twoway_linear <- list()
  pop_ps_vs_twoway_linear <- NA_pop
  time_ps_vs_twoway_linear <- NA
  start_time <- proc.time()
  if(!all(is.na(d_vs))){
    dt_twoway_ps_vs <- cbind(dt_aou, weights = d_vs)
    clus_ps_vs_twoway <- svydesign(id=~1, weights = ~weights, fpc =~fpc, data= dt_twoway_ps_vs)
    margin_names <- cal_names(as.formula(paste0("~",paste0(names_univar_new,collapse = "+"))),clus_ps_vs_twoway)
    transformed_names <- sapply(margin_names, transform_term)
    poptotals_ps_vs <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  }
  
  tryCatch(
    {
      clus_ps_twoway_linear_vs_calirated <- survey::calibrate(clus_ps_vs_twoway, 
                                                              as.formula(paste0("~",paste0(names_univar_new,collapse = "+"))), 
                                                              poptotals_ps_vs,
                                                              calfun = "linear",maxit = 1e2) 
      pop_ps_vs_twoway_linear <- fun.out(y_aou,weights(clus_ps_twoway_linear_vs_calirated),
                                         temp_y,
                                         names_pop,
                                         clus_ps_twoway_linear_vs_calirated)
      end_time <- proc.time()
      time_ps_vs_twoway_linear <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_ps_vs_twoway_linear$pop <- pop_ps_vs_twoway_linear
  result_ps_vs_twoway_linear$time <- time_ps_vs_twoway_raw + time_ps_vs_twoway_linear
  out[["ps_vs_twoway_linear"]] <- result_ps_vs_twoway_linear
  
  ### Raking
  result_ps_vs_twoway_rake <- list()
  pop_ps_vs_twoway_rake <- NA_pop
  time_ps_vs_twoway_rake <- NA
  start_time <- proc.time()
  tryCatch(
    {
      clus_ps_vs_twoway_rake <- survey::calibrate(clus_ps_vs_twoway, 
                                                  as.formula(paste0("~",paste0(names_univar_new,collapse = "+"))),
                                                  poptotals_ps_vs,
                                                  calfun = "raking",maxit = 1e2,
                                                  epsilon = rep(1e-7,length(poptotals_ps_vs))) 
      pop_ps_vs_twoway_rake <- fun.out(y_aou,
                                       weights(clus_ps_vs_twoway_rake),
                                       temp_y,names_pop,
                                       clus_ps_vs_twoway_rake)
      end_time <- proc.time()
      time_ps_vs_twoway_rake <- (end_time - start_time)["elapsed"]
    },
    error = function(e) e
  )
  result_ps_vs_twoway_rake$pop <- pop_ps_vs_twoway_rake
  result_ps_vs_twoway_rake$time <- time_ps_vs_twoway_raw + time_ps_vs_twoway_rake
  out[["ps_vs_twoway_rake"]] <- result_ps_vs_twoway_rake
  
  ##### M9: Two-way, Pseudo-likelihood-based Propensity Score Weighting; Variable Selection + Selection For Convergence, raw, linear and raking #####
  result_ps_vs_stepwise_twoway_rake <- list()
  start_time <- proc.time()
  n_uni <- length(names_univar)
  n_max <- length(names_univar_new)
  index_con <- 1
  pop_ps_vs_stepwise_twoway_rake <- NA_pop
  time_ps_vs_stepwise_twoway_rake <- NA
  n_stepwise <- rep(NA, 3)
  while(!is.na(index_con) & n_uni < n_max){
    n_uni <- n_uni + 1
    cat_temp <- formula(paste0("~",paste0(names_univar_new[seq(n_uni)],collapse = "+")))
    temp <- dt_aou
    temp_out <- fun.vs(cat_temp,dt_aou,dt_s)
    temp_d <- temp_out$d
    theta <- temp_out$theta
    dt_stepwise <- cbind(temp, weights = temp_d)
    
    clus_stepwise <- svydesign(id=~1, weights = ~weights, data= dt_stepwise)
    
    margin_names <- cal_names(as.formula(paste0("~",paste0(names_univar_new[seq(n_uni)],collapse = "+"))),clus_stepwise)
    transformed_names <- sapply(margin_names, transform_term)
    poptotals_stepwise <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
    tryCatch(
      {
        clus_stepwise_calibrated <- survey::calibrate(clus_stepwise, 
                                                      cat_temp,
                                                      poptotals_stepwise,
                                                      calfun = "raking",maxit = 1e2,
                                                      epsilon = rep(1e-7,length(poptotals_stepwise))) 
        pop_ps_vs_stepwise_twoway_rake <- fun.out(y_aou,
                                                  weights(clus_stepwise_calibrated),
                                                  
                                                  temp_y,names_pop,
                                                  clus_stepwise_calibrated)
        end_time <- proc.time()
        time_ps_vs_stepwise_twoway_rake <- (end_time - start_time)["elapsed"]
        n_stepwise <- c(total = n_max,
                        unselect = n_max - n_uni,
                        select = n_uni - length(names_univar)
        )
      },
      error = function(e) {
        index_con <<- NA
      }
    )
  }
  result_ps_vs_stepwise_twoway_rake$pop <- pop_ps_vs_stepwise_twoway_rake
  result_ps_vs_stepwise_twoway_rake$time <- time_ps_vs_stepwise_twoway_rake
  result_ps_vs_stepwise_twoway_rake$n_stepwise <- n_stepwise
  out[["ps_vs_stepwise_twoway_rake"]] <- result_ps_vs_stepwise_twoway_rake
  
  ### Variance Correction
  result_varcor <- NA
  tryCatch(
    {
      result_varcor <- fun.rails.var(dt_aou = dt_aou,
                                     dt_s = dt_s,
                                     y_aou = y_aou,
                                     cat_temp = cat_temp,
                                     w_NP =weights(clus_stepwise_calibrated),
                                     d_P = dt_s$w,
                                     pi_NP = temp_out$pi_aou,
                                     pi_P = temp_out$pia,
                                     T_pop = poptotals_stepwise,
                                     tempy = temp_y
      )
    },error = function(e) e
  )
  out[["var_correction"]] <- result_varcor
  
  
  # M10. No variable selection, use oracle model
  result_oracle_rails <- list()
  start_time <- proc.time()
  pop_oracle_rails <- NA_pop
  time_oracle_rails <- NA
  tryCatch(
    {
      cat_temp <- true_model
      temp <- dt_aou
      temp_out <- fun.vs(cat_temp,dt_aou,dt_s)
      temp_d <- temp_out$d
      theta <- temp_out$theta
      dt_oracle_rails <- cbind(temp, weights = temp_d)
      
      clus_oracle_rails <- svydesign(id=~1, weights = ~weights, data= dt_oracle_rails)
      
      margin_names <- cal_names(true_model,clus_stepwise)
      transformed_names <- sapply(margin_names, transform_term)
      poptotals_oracle_rails <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
      clus_oracle_rails <- survey::calibrate(clus_oracle_rails, 
                                             cat_temp,
                                             poptotals_oracle_rails,
                                             calfun = "raking",maxit = 1e2,
                                             epsilon = rep(1e-7,length(poptotals_oracle_rails))) 
      pop_oracle_rails <- fun.out(y_aou,
                                  weights(clus_oracle_rails),
                                  
                                  temp_y,names_pop,
                                  clus_stepwise_calibrated)
      end_time <- proc.time()
      time_oracle_rails <- (end_time - start_time)["elapsed"]
      
    },
    error = function(e) {
      index_con <<- NA
    }
  )
  result_oracle_rails$pop <- pop_oracle_rails
  result_oracle_rails$time <- time_oracle_rails
  out[["oracle_rails"]] <- result_oracle_rails
  
  
  
  ### Variance without VS
  result_varcor_oracle_rails <- NA
  tryCatch(
    {
      cat_temp <- true_model
      temp <- dt_aou
      temp_out <- fun.vs(cat_temp,dt_aou,dt_s)
      temp_d <- temp_out$d
      theta <- temp_out$theta
      dt_oracle_rails <- cbind(temp, weights = temp_d)
      result_varcor_oracle_rails <- fun.rails.var(dt_aou = dt_aou,
                                     dt_s = dt_s,
                                     y_aou = y_aou,
                                     cat_temp = cat_temp,
                                     w_NP =weights(clus_oracle_rails),
                                     d_P = dt_s$w,
                                     pi_NP = temp_out$pi_aou,
                                     pi_P = temp_out$pia,
                                     T_pop = poptotals_oracle_rails,
                                     tempy = temp_y
      )
    },error = function(e) e
  )
  out[["var_correction_oracle_rails"]] <- result_varcor_oracle_rails
  
  #### Raking with stacked variance
  #### Univariate only
  result_varcor_svy_oneway_rake <- NA
  cat_temp <- as.formula(paste0("~",paste0(names_univar,collapse = "+")))
  tryCatch(
    {
      result_varcor_svy_oneway_rake <- fun.svy.var(dt_aou = dt_aou,
                                                   y_aou = y_aou,
                                                   cat_temp = cat_temp,
                                                   w_NP = weights(clus_oneway_calirated_rake),
                                                   T_pop = poptotals_oneway,
                                                   tempy = temp_y
        
      )
    },error = function(e) e
  )
  out[["var_correction_svy_oneway_rake"]] <- result_varcor_svy_oneway_rake
  
  out$out_size <- ntaou
  
  ##### Number of Zero Cells
  out$zero_cells <- index_zerocell
  
  return(out)
}
