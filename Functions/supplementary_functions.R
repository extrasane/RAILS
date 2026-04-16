# Supplementary Functions for the main script

#### expit function
expit <- function(x){
  exp(x)/(1+exp(x))
}

#### String Transformation
transform_term <- function(term) {
  if (term == "(Intercept)") {
    return("`(Intercept)`= nsiz")
  } else if (str_detect(term, ":")) {
    # Handle interaction terms
    parts <- str_split(term, ":")[[1]]
    transformed_parts <- lapply(parts, function(part) {
      match <- str_match(part, "catx(\\d)(\\d)")
      var <- match[2]
      val <- match[3]
      paste0("dt$catx", var, " == ", val)
    })
    return(paste0("`", term, "` = sum(", paste(transformed_parts, collapse = " & "), ")"))
  } else {
    # Handle single catx terms
    match <- str_match(term, "catx(\\d)(\\d)")
    var <- match[2]
    val <- match[3]
    return(paste0("`", term, "` = sum(dt$catx", var, " == ", val, ")"))
  }
}

transform_term_aou <- function(term) {
  if (term == "(Intercept)") {
    return("`(Intercept)`= ntaou")
  } else if (str_detect(term, ":")) {
    # Handle interaction terms
    parts <- str_split(term, ":")[[1]]
    transformed_parts <- lapply(parts, function(part) {
      match <- str_match(part, "catx(\\d)(\\d)")
      var <- match[2]
      val <- match[3]
      paste0("dt_aou$catx", var, " == ", val)
    })
    return(paste0("`", term, "` = sum(", paste(transformed_parts, collapse = " & "), ")"))
  } else {
    # Handle single catx terms
    match <- str_match(term, "catx(\\d)(\\d)")
    var <- match[2]
    val <- match[3]
    return(paste0("`", term, "` = sum(dt_aou$catx", var, " == ", val, ")"))
  }
}


##### Variable Selection Function
fun.vs <- function(cat_temp,dt_aou,dt_s){
  m_aou <- model.matrix(cat_temp, dt_aou)
  m_s <- model.matrix(cat_temp, dt_s)
  theta <- rep(0,dim(m_aou)[2])
  pia <- rep(1/2,dim(dt_s)[1])
  U1 <- colSums(m_aou) 
  res <- 1
  while(res>=1e-10){
    Uscore <- U1 - colSums(dt_s$w * pia * m_s)
    Htheta <- t(m_s) %*% (dt_s$w * pia * (1 - pia)*m_s)
    theta0 <- theta + solve(Htheta) %*% Uscore
    pia <- as.numeric(m_s %*% theta0)
    pia <- exp(pia)/(1+exp(pia))
    #pia <- pt(pia,df = length(Uscore) - 1,lower.tail = T)
    res <- sum((solve(Htheta) %*% Uscore)^2)
    theta <- theta0
    LKHD <- sum(m_aou %*% theta) - sum(dt_s$w*log(1+exp(m_s %*% theta)))
  }
  pi_aou <- as.numeric(m_aou %*% theta)
  pi_aou <- exp(pi_aou)/(1+exp(pi_aou))
  EENP_1 <- colSums(m_aou) 
  EEP_1 <- - (t(dt_s$w*m_s) %*% pia)
  d <-  1/pi_aou
  d <- d/sum(d)*nsiz
  return(list(d = d, LKHD = LKHD,EEP_1 = EEP_1,EENP_1 = EENP_1,theta = theta,
              pi_aou = pi_aou,pia = pia))
}

##### Forward selection based likelihood ratio test
### For each two-way interaction, return the likelihood ratio test result
fun.lkd <- function(x,names_univar,LKHD,dt_s,dt_aou, m_s){
  #### x as a newly added variable to the model
  cat_formula <- formula(paste0("~",paste0(paste0(names_univar,collapse = "+"),"+",x)))
  t_s <- model.matrix(cat_formula, dt_s)
  t_aou <- model.matrix(cat_formula, dt_aou)
  theta <- rep(0,dim(t_aou)[2])
  pia <- rep(1/2,dim(dt_s)[1])
  U1 <- colSums(t_aou) 
  if(any(U1==0)){
    return(rep(NA,5))
  }
  res <- 1
  res0 <- 0
  #while(abs(res-res0)>=1e-8){
  tryCatch(
    {
      while(res>=1e-10){
        res0 <- res
        Uscore <- U1 - colSums(dt_s$w * pia * t_s)
        Htheta <- t(t_s) %*% (dt_s$w * pia * (1 - pia)*t_s)
        theta0 <- theta + solve(Htheta) %*% Uscore
        pia <- as.numeric(t_s %*% theta0)
        pia <- exp(pia)/(1+exp(pia))
        #pia <- pt(pia,df = length(Uscore) - 1,lower.tail = T)
        res <- sum((solve(Htheta) %*% Uscore)^2)
        theta <- theta0
        temp_LKHD <- sum(t_aou %*% theta) - sum(dt_s$w*log(1+exp(t_s %*% theta)))
      }
      # Likelihood ratio test based on difference in variables numbers
      LRT <- 2 * (temp_LKHD - LKHD)
      pval <- 1 - pchisq(LRT,df = dim(t_s)[2]- dim(m_s)[2] )
      return(c(temp_LKHD,
               LRT,
               pval,
               LRT/(dim(t_s)[2]- dim(m_s)[2]),
               dim(t_s)[2]- dim(m_s)[2]
      ))
    }, error = function(e) {
      return(rep(NA,5))
    }
  )
}

#==============================================================================
# RAILS Variance Estimation Function
#==============================================================================

fun.rails.var <- function(
    # Data
  dt_aou,           # Biobank data frame
  dt_s,             # Reference data frame
  y_aou,            # Outcome vector
  
  # Formula
  cat_temp,         # Formula for PS and calibration (~x1 + x2 + ...)
  
  # Weights and PS
  w_NP,             # Final calibrated weights
  d_P,              # Survey weights for reference
  pi_NP,            # Propensity scores for biobank
  pi_P,             # Propensity scores for reference
  
  # Population totals
  T_pop,            # Calibration totals
  
  # True values for coverage
  tempy = NULL      # Current simulation mean
) {
  
  #============================================================================
  # SETUP
  #============================================================================
  
  # Ensure vectors
  y_aou <- as.vector(y_aou)
  w_NP  <- as.vector(w_NP)
  pi_NP <- as.vector(pi_NP)
  pi_P  <- as.vector(pi_P)
  d_P   <- as.vector(d_P)
  
  n_NP   <- length(w_NP)
  n_P    <- length(d_P)
  N_hat  <- sum(w_NP)
  mu_hat <- sum(w_NP * y_aou) / N_hat
  
  # Design matrices
  X_ps_NP <- model.matrix(cat_temp, dt_aou)
  X_ps_P  <- model.matrix(cat_temp, dt_s)
  X_cal   <- model.matrix(cat_temp, dt_aou)
  
  p <- ncol(X_ps_NP)
  q <- ncol(X_cal)
  
  #============================================================================
  # INDIVIDUAL PSI VECTORS
  #============================================================================
  
  # Reference (P): only Stage 1 is non-zero
  psi1_P <- -X_ps_P * d_P * pi_P
  
  # Biobank (NP)
  psi1_NP <- X_ps_NP
  
  psi2_NP <- X_cal * w_NP -
    matrix(T_pop / n_NP, nrow = n_NP, ncol = q, byrow = TRUE)
  
  # CORRECTED: psi3 = w_NP * (y_i - mu_hat) / N_hat
  # Previously: w_NP * (y_i - mu_hat)  [missing 1/N_hat scaling]
  psi3_NP <- w_NP * (y_aou - mu_hat) / N_hat
  
  #============================================================================
  # MEAT MATRIX M
  #============================================================================
  
  M_11 <- crossprod(psi1_P) + crossprod(psi1_NP)
  M_12 <- crossprod(psi1_NP, psi2_NP)
  M_13 <- crossprod(psi1_NP, psi3_NP)
  M_22 <- crossprod(psi2_NP)
  M_23 <- crossprod(psi2_NP, psi3_NP)
  M_33 <- sum(psi3_NP^2)
  # Note: M_33 = sum( w_NP^2 * (y_aou - mu_hat)^2 ) / N_hat^2
  
  #============================================================================
  # BREAD MATRIX B
  #============================================================================
  
  # B_11: positive definite (unchanged)
  B_11 <- crossprod(X_ps_P * sqrt(d_P * pi_P * (1 - pi_P)))
  
  # B_22: negative definite (unchanged)
  B_22 <- -crossprod(X_cal * sqrt(w_NP))
  
  # CORRECTED: B_33 = 1  (was N_hat; dividing Psi_3 by N_hat makes the
  # derivative w.r.t. mu equal to -sum(w_NP)/N_hat = -1, so B_33 = 1)
  B_33 <- 1
  
  # B_21: unchanged
  B_21 <- t(X_cal) %*% (X_ps_NP * w_NP * (1 - pi_NP))
  
  # CORRECTED: B_31 gains a 1/N_hat factor
  # Previously: multiplier_31 <- w_NP * (1 - pi_NP) * (y_aou - mu_hat)
  multiplier_31 <- w_NP * (1 - pi_NP) * (y_aou - mu_hat) / N_hat
  B_31 <- matrix(colSums(sweep(X_ps_NP, 1, multiplier_31, "*")), nrow = 1)
  
  # CORRECTED: B_32 gains a 1/N_hat factor
  # Previously: multiplier_32 <- w_NP * (y_aou - mu_hat)
  multiplier_32 <- w_NP * (y_aou - mu_hat) / N_hat
  B_32 <- matrix(-colSums(sweep(X_cal, 1, multiplier_32, "*")), nrow = 1)
  
  #============================================================================
  # BREAD INVERSE
  #============================================================================
  
  B_11_inv <- solve(B_11)
  B_22_inv <- solve(B_22)
  
  # CORRECTED: C_33 = 1/B_33 = 1
  C_33 <- 1 / B_33          # = 1
  C_32 <- -C_33 * B_32 %*% B_22_inv
  C_31 <- -C_33 * (B_31 - B_32 %*% B_22_inv %*% B_21) %*% B_11_inv
  
  #============================================================================
  # VARIANCE COMPONENTS
  #============================================================================
  
  term_33 <- C_33^2 * M_33                               # Naive = M_33
  term_11 <- as.numeric(C_31 %*% M_11 %*% t(C_31))      # Stage 1
  term_22 <- as.numeric(C_32 %*% M_22 %*% t(C_32))      # Stage 2
  term_12 <- as.numeric(2 * C_31 %*% M_12 %*% t(C_32))  # Cross 1-2
  term_13 <- as.numeric(2 * C_33 * C_31 %*% M_13)       # Cross 1-3
  term_23 <- as.numeric(2 * C_33 * C_32 %*% M_23)       # Cross 2-3
  
  var_naive <- term_33
  var_full  <- term_33 + term_11 + term_22 + term_12 + term_13 + term_23
  
  # Handle potential negative variance (numerical issues)
  if (var_full < 0) {
    warning("Full variance is negative! Using absolute value.")
    var_full <- abs(var_full)
  }
  
  se_naive <- sqrt(var_naive)
  se_full  <- sqrt(var_full)
  
  #============================================================================
  # CONFIDENCE INTERVALS
  #============================================================================
  
  z_alpha <- qnorm(0.975)
  
  ci_naive_lower <- mu_hat - z_alpha * se_naive
  ci_naive_upper <- mu_hat + z_alpha * se_naive
  
  ci_full_lower  <- mu_hat - z_alpha * se_full
  ci_full_upper  <- mu_hat + z_alpha * se_full
  
  #============================================================================
  # COVERAGE (if true values provided)
  #============================================================================
  

  
  if (!is.null(tempy)) {
    cover_naive_temp <- (tempy >= ci_naive_lower) & (tempy <= ci_naive_upper)
    cover_full_temp  <- (tempy >= ci_full_lower)  & (tempy <= ci_full_upper)
  } else {
    cover_naive_temp <- NA
    cover_full_temp  <- NA
  }
  
  #============================================================================
  # OUTPUT
  #============================================================================
  
  result <- c(
    # Point estimate
    mu_hat = mu_hat,
    
    # Naive variance
    var_naive      = var_naive,
    se_naive       = se_naive,
    ci_naive_lower = ci_naive_lower,
    ci_naive_upper = ci_naive_upper,
    
    # Full variance
    var_full      = var_full,
    se_full       = se_full,
    ci_full_lower = ci_full_lower,
    ci_full_upper = ci_full_upper,
    
    # Variance ratio
    ratio_full_naive = se_full / se_naive,
    
    
    # Coverage (temp/simulation mean)
    cover_naive_temp = as.numeric(cover_naive_temp),
    cover_full_temp  = as.numeric(cover_full_temp),
    
    # Variance components
    term_33 = term_33,
    term_11 = term_11,
    term_22 = term_22,
    term_12 = term_12,
    term_13 = term_13,
    term_23 = term_23,
    
    # Other info
    N_hat = N_hat,
    n_NP  = n_NP,
    n_P   = n_P
  )
  
  return(result)
}

fun.svy.var <- function(
    # Data
  dt_aou,           # Biobank data frame
  y_aou,            # Outcome vector
  
  # Formula
  cat_temp,         # Formula for calibration (~x1 + x2 + ...)
  
  # Weights
  w_NP,             # Final calibrated weights
  
  # Population totals
  T_pop,            # Calibration totals
  
  # True values for coverage
  tempy = NULL      # Current simulation mean
) {
  
  #============================================================================
  # SETUP
  #============================================================================
  
  # Ensure vectors
  y_aou <- as.vector(y_aou)
  w_NP  <- as.vector(w_NP)
  
  n_NP   <- length(w_NP)
  N_hat  <- sum(w_NP)
  mu_hat <- sum(w_NP * y_aou) / N_hat
  
  # Design matrix
  X_cal <- model.matrix(cat_temp, dt_aou)
  q     <- ncol(X_cal)
  
  #============================================================================
  # INDIVIDUAL PSI VECTORS
  #============================================================================
  
  psi2_NP <- X_cal * w_NP -
    matrix(T_pop / n_NP, nrow = n_NP, ncol = q, byrow = TRUE)
  
  psi3_NP <- w_NP * (y_aou - mu_hat) / N_hat
  
  #============================================================================
  # MEAT MATRIX M
  #============================================================================
  
  M_22 <- crossprod(psi2_NP)
  M_23 <- crossprod(psi2_NP, psi3_NP)
  M_33 <- sum(psi3_NP^2)
  
  #============================================================================
  # BREAD MATRIX B
  #============================================================================
  
  # B_22: negative definite
  B_22 <- -crossprod(X_cal * sqrt(w_NP))
  
  # B_33 = 1
  B_33 <- 1
  
  # B_32
  multiplier_32 <- w_NP * (y_aou - mu_hat) / N_hat
  B_32 <- matrix(-colSums(sweep(X_cal, 1, multiplier_32, "*")), nrow = 1)
  
  #============================================================================
  # BREAD INVERSE
  #============================================================================
  
  B_22_inv <- solve(B_22)
  
  C_33 <- 1 / B_33   # = 1
  C_32 <- -C_33 * B_32 %*% B_22_inv
  
  #============================================================================
  # VARIANCE COMPONENTS
  #============================================================================
  
  term_33 <- C_33^2 * M_33
  term_22 <- as.numeric(C_32 %*% M_22 %*% t(C_32))
  term_23 <- as.numeric(2 * C_33 * C_32 %*% M_23)
  
  var_naive <- term_33
  var_full  <- term_33 + term_22 + term_23
  
  if (var_full < 0) {
    warning("Full variance is negative! Using absolute value.")
    var_full <- abs(var_full)
  }
  
  se_naive <- sqrt(var_naive)
  se_full  <- sqrt(var_full)
  
  #============================================================================
  # CONFIDENCE INTERVALS
  #============================================================================
  
  z_alpha <- qnorm(0.975)
  
  ci_naive_lower <- mu_hat - z_alpha * se_naive
  ci_naive_upper <- mu_hat + z_alpha * se_naive
  
  ci_full_lower  <- mu_hat - z_alpha * se_full
  ci_full_upper  <- mu_hat + z_alpha * se_full
  
  #============================================================================
  # COVERAGE (if true values provided)
  #============================================================================

  if (!is.null(tempy)) {
    cover_naive_temp <- (tempy >= ci_naive_lower) & (tempy <= ci_naive_upper)
    cover_full_temp  <- (tempy >= ci_full_lower)  & (tempy <= ci_full_upper)
  } else {
    cover_naive_temp <- NA
    cover_full_temp  <- NA
  }
  
  #============================================================================
  # OUTPUT
  #============================================================================
  
  result <- c(
    # Point estimate
    mu_hat = mu_hat,
    
    # Naive variance
    var_naive      = var_naive,
    se_naive       = se_naive,
    ci_naive_lower = ci_naive_lower,
    ci_naive_upper = ci_naive_upper,
    
    # Full variance
    var_full      = var_full,
    se_full       = se_full,
    ci_full_lower = ci_full_lower,
    ci_full_upper = ci_full_upper,
    
    # Variance ratio
    ratio_full_naive = se_full / se_naive,
    
    # Coverage (temp/simulation mean)
    cover_naive_temp = as.numeric(cover_naive_temp),
    cover_full_temp  = as.numeric(cover_full_temp),
    
    # Variance components
    term_33 = term_33,
    term_22 = term_22,
    term_23 = term_23,
    
    # Other info
    N_hat = N_hat,
    n_NP  = n_NP
  )
  
  return(result)
}