fun.rep <- function(x = sim_result,n){
  # x as the input list
  # v as the variable names to be reported
  # n as the number of simulations
  x <- sim_result

  
  # Names of the variables
  names_var <- colnames(x[[1]]$unweighted$sub)
  n_var <- length(names_var)
  
  # Names of the Methods
  names_method <- names(x[[1]])[1:(length(x[[1]]))]
  names_method <- subset(names_method, names_method != "var_correction" 
                         & names_method != "avg_size" 
                         & names_method != "error" 
                         & names_method != "zero_cells"
                         & names_method !="var_correction_oracle_rails"
                         & names_method !="var_correction_svy_oneway_rake")
  names_method <- c(
    names_method[names_method != "ps_vs_stepwise_twoway_rake"],
    "ps_vs_stepwise_twoway_rake"
  )
  n_method <- length(names_method)
  names_method_fix <- c("naive","oracle","G-Raking","NPS","GVS",
                        "G-RAILS","S-Raking","S-RAILS")
  # Index for population estimate
  names_popest <- names(x[[1]]$unweighted$pop)
  n_popest <- length(names_popest)
  
  # Index for subgroup estimate
  names_subest <- rownames(x[[1]]$unweighted$sub)
  n_subest <- length(names_subest)
  
  # Collect the number of error messages
  error_count <- c()
  error_index <- c()
  for (i in 1:n){
    if(!is.null(x[[i]]$error)){
      error_count <- c(error_count,x[[i]]$error)
      error_index <- c(error_index,i)
    }
  }
  error_count <- as.vector(table(error_count))
  if(length(error_count)==0){
    error_count <- 0
  }
  
  
  index <- c()
  s = 0
  # Correct VS_STEPWISE if non-convergent
  for(i in 1:n){
    if(is.na(x[[i]]$ps_vs_stepwise_twoway_rake[1])){
      index[s] <- i
      s = s + 1
      x[[i]]$ps_vs_stepwise_twoway_rake <- x[[i]]$ps_vs_twoway_rake
    }
  }
  
  
  # Population Prevalence
  count_non_na <- numeric(length = n_method)
  sum_pop <- matrix(0,ncol = n_popest,nrow = n_method)
  est_pop <- matrix(0,nrow = n,ncol = n_method)
  sd_pop <- matrix(0,nrow = n,ncol = n_method)
  colnames(sum_pop) <- names_popest
  colnames(est_pop) <- names_method
  colnames(sd_pop) <- names_method
  rownames(sum_pop) <- names_method
  for(i in setdiff(1:n,error_index)){
    for(j in 1:n_method){
      temp <- x[[i]][[names_method[j]]]$pop
      count_non_na[j] <- count_non_na[j] + all(!is.na(temp[1]))
      sum_pop[j,] <- sum_pop[j,] + ifelse((is.na(temp)),0,temp)
      est_pop[i,j] <-  ifelse(is.na(temp[1]),NA,temp[1])
      sd_pop[i,j] <-  ifelse(is.na(temp[2]),NA,sqrt(temp[2]))
    }
  }
  out_mean <-  sum_pop/count_non_na
  out_mean <- cbind(out_mean, "Ratio of NAs" = 1 - count_non_na/n)
  out_mean <- as.data.frame(out_mean)
  empsd <-      apply(est_pop,2,sd,na.rm = TRUE)
  est_median <- apply(est_pop,2,median,na.rm = TRUE)
  exp_sd <-     apply(sd_pop,2,mean,na.rm = TRUE)
  sd_median <-  apply(sd_pop,2,median,na.rm = TRUE)
  out_mad <-    apply(est_pop,2,mad,na.rm = TRUE)
  
  truey <- out_mean["trueweight",1]
  
  # Bias v.s. true y
  out_mean <- out_mean %>%
    mutate(median = est_median,.after = "mean") %>%
    mutate("Bias(Mean - True)" = mean - truey,.after = bias2) %>%
    mutate("Relative Bias(Mean)" = (mean - truey)/truey,.after = "median") %>%
    mutate("Relative Bias(Median)" = (median - truey)/truey,.after = "Relative Bias(Mean)") %>%
    mutate("Avg SD" = exp_sd,.after = "Relative Bias(Median)") %>%
    mutate("Emp SD" = empsd,.after = "Avg SD") %>%
    mutate("Median SD" = sd_median,.after = "Emp SD") %>%
    mutate("MAD" = out_mad,.after = "Emp SD")
  
  ### Summarize modified variance 
  out_varcor <- vector("list", length = n)
  out_var_mod <- vector("list", length = n)
  out_varcor_oracle <- vector("list", length = n)
  out_varcor_svy <- vector("list", length = n)
  
  for(i in 1:n){
    # var_correction
    temp <- x[[i]]$var_correction
    out_varcor[[i]] <- ifelse(is.na(temp), rep(NA, 24), temp)
    out_var_mod[[i]] <- ifelse(is.na(temp[1:4]), NA,
                               c(sum(temp[c(16,17,20)]),
                                 sqrt(sum(temp[c(16,17,20)])),
                                 sum(temp[c(16,17,20)]) / temp[2],
                                 sqrt(sum(temp[c(16,17,20)])) / temp[3]))
    
    # var_correction_oracle_rails
    temp_oracle <- x[[i]]$var_correction_oracle_rails
    out_varcor_oracle[[i]] <- ifelse(is.na(temp_oracle), rep(NA, 24), temp_oracle)
    
    # var_correction_svy_oneway_rake
    temp_svy <- x[[i]]$var_correction_svy_oneway_rake
    out_varcor_svy[[i]] <- ifelse(is.na(temp_svy), rep(NA, 24), temp_svy)
  }
  
  # var_correction summaries
  index_varcor <- !sapply(out_varcor, function(x) length(x) == 1 && is.na(x))
  res_varcor   <- Reduce(`+`, out_varcor[index_varcor]) / sum(index_varcor)
  pos_count    <- Reduce(`+`, lapply(out_varcor[index_varcor], function(x) x > 0))
  mod_count    <- Reduce(`+`, out_var_mod[index_varcor]) / sum(index_varcor)
  names(mod_count) <- c("var_mod", "se_mod", "ratio_var_mod", "ratio_se_mod")
  
  # var_correction_oracle_rails summaries
  index_varcor_oracle <- !sapply(out_varcor_oracle, function(x) length(x) == 1 && is.na(x))
  res_varcor_oracle   <- Reduce(`+`, out_varcor_oracle[index_varcor_oracle]) / sum(index_varcor_oracle)
  pos_count_oracle    <- Reduce(`+`, lapply(out_varcor_oracle[index_varcor_oracle], function(x) x > 0))
  
  # var_correction_svy_oneway_rake summaries
  index_varcor_svy <- !sapply(out_varcor_svy, function(x) length(x) == 1 && is.na(x))
  res_varcor_svy   <- Reduce(`+`, out_varcor_svy[index_varcor_svy]) / sum(index_varcor_svy)
  pos_count_svy    <- Reduce(`+`, lapply(out_varcor_svy[index_varcor_svy], function(x) x > 0))
  first_elements <- sapply(out_varcor, `[[`, 1)
  emp_sd_varcor <- sd(first_elements, na.rm = TRUE)
  
  # Oracle Coverage
  index_svy_oneway_rake <- which(names_method == "svy_oneway_rake")
  index_oracle_rails <- which(names_method == "oracle_rails")
  safe_get <- function(obj, idx = 6) {
    if (is.null(obj) || length(obj) < idx || all(is.na(obj))) NA else obj[[idx]]
  }
  out_oracov <- numeric(length = n_method + 4)
  for(i in 1:n){
    for(j in 1:(n_method + 4)){
      if(j < (n_method + 1)){
        temp <- x[[i]][[names_method[j]]]$pop
        debias_est <- temp[1] - out_mean$`Bias(Mean - True)`[j]
        temp_ci <- debias_est + c(-1, 1) * 1.96 * sqrt(temp[2])
        temp_oracov <- truey >= temp_ci[1] & truey <= temp_ci[2]
        out_oracov[j] <- out_oracov[j] + ifelse(is.na(temp_oracov), 0, temp_oracov)
      } else {
        temp <- if (j == n_method + 1) {
          x[[i]][[names_method[index_svy_oneway_rake]]]$pop
        } else if (j == n_method + 2) {
          x[[i]][[names_method[n_method]]]$pop        # Original RAILS
        } else if (j == n_method + 3) {
          x[[i]][[names_method[index_oracle_rails]]]$pop
        } else {
          # j == n_method + 4: Raking_1 Survey variance
          x[[i]][[names_method[index_svy_oneway_rake]]]$pop
        }
        
        temp_se <- if (j == n_method + 1) {
          safe_get(x[[i]][["var_correction_svy_oneway_rake"]])
        } else if (j == n_method + 2) {
          safe_get(x[[i]][["var_correction"]])
        } else if (j == n_method + 3) {
          safe_get(x[[i]][["var_correction_oracle_rails"]])
        } else {
          # j == n_method + 4: use SvyVar directly (temp[2])
          temp[2]
        }
        
        debias_est <- if (j == n_method + 1) {
          temp[1] - out_mean$`Bias(Mean - True)`[index_svy_oneway_rake]
        } else if (j == n_method + 2) {
          temp[1] - out_mean$`Bias(Mean - True)`[n_method]
        } else if (j == n_method + 3) {
          temp[1] - out_mean$`Bias(Mean - True)`[index_oracle_rails]
        } else {
          # j == n_method + 4
          temp[1] - out_mean$`Bias(Mean - True)`[index_svy_oneway_rake]
        }
        
        temp_ci <- debias_est + c(-1, 1) * 1.96 * sqrt(temp_se)
        temp_oracov <- truey >= temp_ci[1] & truey <= temp_ci[2]
        out_oracov[j] <- out_oracov[j] + ifelse(is.na(temp_oracov), 0, temp_oracov)
      }
    }
  }
  
  count_varcor <- sum(index_varcor)
  
  count_non_na_2 <- c(count_non_na,
                      count_non_na[index_svy_oneway_rake],   # n_method + 1: svy stacked
                      count_varcor,                           # n_method + 2: orig stacked
                      count_non_na[index_oracle_rails],       # n_method + 3: oracle stacked
                      count_non_na[index_svy_oneway_rake]     # n_method + 4: svy survey
  )
  out_oracov <- out_oracov / count_non_na_2
  out_mean <- out_mean %>%
    mutate("Oracle Coverage" = out_oracov[1:n_method],.after = NomCoverage)
  
  
  # Elapsed Time
  out_time <- numeric(length = n_method)
  names(out_time) <- names_method
  for(i in 1:n){
    for(j in 1:n_method){
      temp <- x[[i]][[names_method[j]]]$time
      out_time[j] <- out_time[j] + ifelse((is.na(temp)),0,temp)
    }
  }
  out_time <- out_time/count_non_na
  
  # Average Size
  out_size <- numeric(length = n_var)
  names(out_size) <- names_var
  for(i in 1:n){
    temp <- x[[i]]$avg_size
    out_size <- out_size + temp
  }
  out_size <- out_size/n
  
  # Average Variable Selected
  out_vs_n <- numeric(2)
  for(i in 1:n){
    temp <- x[[i]]$ps_vs_twoway_raw$n_vs
    out_vs_n <- out_vs_n + ifelse((is.na(temp)),0,temp)
  }
  out_vs_n <- out_vs_n/count_non_na[names_method=="ps_vs_twoway_raw"]
  out_vs_n <- c(total= sum(out_vs_n),out_vs_n)
  # Average Variable for Stepwise Selection
  out_vs_n_step <- numeric(3)
  for(i in 1:n){
    temp <- x[[i]]$ps_vs_stepwise_twoway_rake$n_stepwise
    out_vs_n_step <- out_vs_n_step + ifelse((is.na(temp)),0,temp)
  }
  out_vs_n_step <- out_vs_n_step/count_non_na[names_method=="ps_vs_stepwise_twoway_rake"]
  
  
  # Plot
  bias <- out_mean$`Bias(Mean - True)`
  tt <- data.frame(matrix(0,nrow = n,ncol = n_method))
  colnames(tt) <- names_method
  
  for(i in 1:n){
    for(j in 1:n_method){
      temp <- x[[i]][[names_method[j]]]$pop[1]
      tt[i,j] <-  ifelse((is.na(temp)),NA,temp) - truey
    }
  }
  
  
  ### tt wide to long format
  tmat <- data.frame(Method = names_method,MeanBias = bias) %>% na.omit()
  level_order <- tmat$Method
  tmp <- tt %>%
    as.data.frame() %>%
    select(where(~ !all(is.na(.)))) %>%
    tidyr::pivot_longer(everything(), names_to = "Method",values_to = "value") %>%
    right_join(tmat) %>%
    mutate(Method = factor(Method, levels = level_order)) %>%
    mutate(Method = droplevels(Method))
  
  out_varcor_cov <- out_oracov[n_method : (n_method + 4)]
  names(out_varcor_cov) <- c("original_rails",
                             "svy_one_rake_stacked",
                             "original_rails_stacked",
                             "oracle_rails_stacked",
                             "svy_one_rake_survey")
  
  out <- list(out_mean = out_mean,
              out_time = out_time,
              out_size = out_size,
              out_vs_n = out_vs_n,
              out_vs_n_step = out_vs_n_step,
              out_p = tmp,
              out_level_order = level_order,
              out_varcor = res_varcor,
              out_varcor_oracle = res_varcor_oracle,
              out_varcor_svy = res_varcor_svy,
              out_pos_count = pos_count,
              out_varcor_cov = out_varcor_cov,
              out_mod_count = mod_count,
              out_varcor_empsd = emp_sd_varcor)
  return(out)
}
