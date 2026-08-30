# Supplementary Functions for the main script — refined (cached design matrices)
# Author: Huiding Chen
# Notes: Caches model matrices and updates by column binding to avoid recomputation in
#        greedy selection. Also fixes DF accounting and improves numerical stability.
# - Requires: stringr, dplyr, purrr, stats, Matrix, survey
# - Style: base-R with light tidyverse; defensive checks; zero side effects.

###############################
# 0) Imports & utilities
###############################

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(purrr)
  library(survey)
})

# Stable logistic (expit)
expit <- function(x) 1/(1 + exp(-x))

# Log1pexp for numerical stability
.log1pexp <- function(z) {
  out <- numeric(length(z))
  big <- z > 0
  out[big]  <- z[big] + log1p(exp(-z[big]))
  out[!big] <- log1p(exp(z[!big]))
  out
}

# Solve (Z'Z + λI) v = u, with fallback
.solve_ridge_mat <- function(Z, u, ridge = 1e-8) {
  p <- ncol(Z)
  G <- crossprod(Z) + diag(ridge, p)
  tryCatch({
    drop(chol2inv(chol(G)) %*% u)
  }, error = function(e) {
    drop(solve(G, u, tol = 1e-12))
  })
}

# Truncated Pareto sampler (helper), shape>0, lower>0, optional upper
# rtruncpareto <- function(n, lower, shape, upper = Inf) {
#   stopifnot(n >= 1, lower > 0, shape > 0)
#   U <- runif(n)
#   if (is.finite(upper)) {
#     Fl <- 1 - (lower/upper)^shape
#     x <- lower / (1 - Fl * U)^(1/shape)
#   } else {
#     x <- lower / (U)^(1/shape)
#   }
#   x
# }

###############################
# 1) String transformation helpers
###############################
# These convert model terms to counting statements for margins

.transform_single <- function(term, prefix = "dt$") {
  m <- str_match(term, "catx(\\d+)(\\d+)")
  if (anyNA(m)) stop("Unrecognized term: ", term)
  var <- m[, 2]
  val <- m[, 3]
  paste0("`", term, "` = sum(", prefix, "catx", var, " == ", val, ")")
}

.transform_interaction <- function(term, prefix = "dt$") {
  parts <- str_split(term, ":", simplify = TRUE)
  conds <- apply(parts, 1, function(part) {
    m <- str_match(part, "catx(\\d+)(\\d+)")
    if (anyNA(m)) stop("Unrecognized interaction part: ", part)
    paste0(prefix, "catx", m[,2], " == ", m[,3])
  })
  paste0("`", term, "` = sum(", paste(conds, collapse = " & "), ")")
}

transform_term <- function(term) {
  if (term == "(Intercept)") {
    return("`(Intercept)` = nrow(dt)")
  } else if (str_detect(term, ":")) {
    return(.transform_interaction(term, prefix = "dt$"))
  } else {
    return(.transform_single(term, prefix = "dt$"))
  }
}

transform_term_aou <- function(term) {
  if (term == "(Intercept)") {
    return("`(Intercept)` = nrow(dt_aou)")
  } else if (str_detect(term, ":")) {
    return(.transform_interaction(term, prefix = "dt_aou$"))
  } else {
    return(.transform_single(term, prefix = "dt_aou$"))
  }
}

transform_term_sub <- function(term, univar, j) {
  stopifnot(is.character(univar), length(j) == 1)
  if (term == "(Intercept)") {
    return(paste0("`(Intercept)` = sum(dt$", univar, " == ", j, ")"))
  } else if (str_detect(term, ":")) {
    parts <- str_split(term, ":", simplify = TRUE)
    conds <- apply(parts, 1, function(part) {
      m <- str_match(part, "catx(\\d+)(\\d+)")
      if (anyNA(m)) stop("Unrecognized interaction part: ", part)
      paste0("dt$catx", m[,2], " == ", m[,3])
    })
    return(paste0("`", term, "` = sum(", paste(conds, collapse = " & "), " & dt$", univar, " == ", j, ")"))
  } else {
    m <- str_match(term, "catx(\\d+)(\\d+)")
    if (anyNA(m)) stop("Unrecognized term: ", term)
    var <- m[, 2]; val <- m[, 3]
    return(paste0("`", term, "` = sum(dt$catx", var, " == ", val, " & dt$", univar, " == ", j, ")"))
  }
}

###############################
# 2) Newton scoring for pseudo-likelihood (matrix-input version)
###############################
# Core fitter using cached design matrices: solves U(θ)=0 with NR.

# Formula-input version retained for external calls if needed
fun.nps <- function(cat_formula, dt, dt_aou, dt_s) {
  m_aou <- model.matrix(cat_formula , dt_aou)
  m_s <- model.matrix(cat_formula, dt_s)
  theta <- rep(0,dim(m_aou)[2])
  pia <- rep(1/2,dim(dt_s)[1])
  U1 <- colSums(m_aou)
  w <- as.numeric(dt_s$w)
  res <- 1
  while(res>=1e-10){
    Uscore <- U1 - colSums(w * pia * m_s)
    Htheta <- t(m_s) %*% (w * pia * (1 - pia)*m_s)
    theta0 <- theta + solve(Htheta) %*% Uscore
    pia <- as.numeric(m_s %*% theta0)
    pia <- exp(pia)/(1+exp(pia))
    #pia <- pt(pia,df = length(Uscore) - 1,lower.tail = T)
    res <- sum((solve(Htheta) %*% Uscore)^2)
    theta <- theta0
    LKHD <- sum(m_aou %*% theta) - sum(w*log(1+exp(m_s %*% theta)))
  }
  pi_aou <- as.numeric(m_aou %*% theta)
  pi_aou <- exp(pi_aou)/(1+exp(pi_aou))
  d <-  1/pi_aou
  d <- d/sum(d)*nsiz
  return(list(d = d, LKHD = LKHD))
}

###############################
# 3) Calibration Raking
###############################
fun.oneway.rake <- function(dt,dt_aou,nsiz,ntaou,names_univar,d = 1){
  dt_oneway <- cbind(dt_aou, weights = d)
  clus_oneway <- survey::svydesign(id = ~1, weights = ~weights, fpc = ~fpc, data = dt_oneway)
  margin_names <- cal_names(as.formula(paste0("~", paste0(names_univar, collapse = "+"))), clus_oneway)
  transformed_names <- sapply(margin_names, transform_term)
  poptotals_oneway <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
  clus_oneway_cal <- NULL
  clus_oneway_cal <- tryCatch({
    survey::calibrate(clus_oneway,
                      as.formula(paste0("~", paste0(names_univar, collapse = "+"))), poptotals_oneway,
                      calfun = "raking", maxit = 1e2)
  }, error = function(e) {
    NULL
  })
  return(clus_oneway_cal)
}


###############################
# 4) Greedy forward selection via LRT with cached matrices
###############################
# Caches full two-way design once and selects by adding term columns.
fun.lkd <- function(x,names_univar,LKHD,dt_s,dt_aou, m_s){
  #### x as a newly added variable to the model
  if(is.null(names_univar)){
    cat_formula <- formula(paste0("~",x))
  }else{
    cat_formula <- formula(paste0("~",paste0(paste0(names_univar,collapse = "+"),"+",x)))
  }
  
  t_s <- model.matrix(cat_formula, dt_s)
  t_aou <- model.matrix(cat_formula, dt_aou)
  theta <- rep(0,dim(t_aou)[2])
  pia <- rep(1/2,dim(dt_s)[1])
  U1 <- colSums(t_aou) 
  w <- as.numeric(dt_s$w)
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
        Uscore <- U1 - colSums(w * pia * t_s)
        Htheta <- t(t_s) %*% (w * pia * (1 - pia)*t_s)
        theta0 <- theta + solve(Htheta) %*% Uscore
        pia <- as.numeric(t_s %*% theta0)
        pia <- exp(pia)/(1+exp(pia))
        #pia <- pt(pia,df = length(Uscore) - 1,lower.tail = T)
        res <- sum((solve(Htheta) %*% Uscore)^2)
        theta <- theta0
        temp_LKHD <- sum(t_aou %*% theta) - sum(w*log(1+exp(t_s %*% theta)))
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


# GVS with main-effects start (interactions as candidates)
fun.gvs <- function(names_univar, dt,dt_aou, dt_s, alpha = 0.05, main_effect_start = TRUE) {
  if(main_effect_start == TRUE){
    # Base Model with main effects
    m_s <- model.matrix(as.formula(paste0("~",paste0(names_univar,collapse = "+"))), dt_s)
    m_aou <- model.matrix(as.formula(paste0("~",paste0(names_univar,collapse = "+"))) , dt_aou)
    tryCatch(
      {
        temp <- fun.nps(as.formula(paste0("~",paste0(names_univar,collapse = "+"))),dt,dt_aou,dt_s)
        d <- temp$d
        LKHD <- temp$LKHD
      },
      error = function(e) {
        LKHD <<- NA
      }
    )
    if(is.na(LKHD)){
      return(NA)
    }
  }else{
    # Base Model with main effects
    m_s <- model.matrix(~1, dt_s)
    m_aou <- model.matrix(~1 , dt_aou)
    tryCatch(
      {
        temp <- fun.nps(~1,dt,dt_aou,dt_s)
        d <- temp$d
        LKHD <- temp$LKHD
      },
      error = function(e) {
        LKHD <<- NA
      }
    )
    if(is.na(LKHD)){
      return(NA)
    }
  }
  
  # Greedy Variable Selection
  temp <- t(combn(names_univar,2))
  twovars <- temp %>% apply(1,paste,collapse = ":")
  # twovars_new <- twovars
  # names_univar_new <- names_univar
  
  
  if(main_effect_start == TRUE){
    pool_vars <- twovars
    cur_vars <- names_univar
    
  }else{
    pool_vars <- c(names_univar,twovars)
    cur_vars <- NULL
  }
  LKHD_new <- LKHD
  m_s_new <- m_s
  index_add <- 1
  
  while(index_add !=0){
    if(length(pool_vars)==0){
      index_add <- 0
      break
    }
    temp_p <- apply(as.data.frame(pool_vars),1,fun.lkd,
                    names_univar = cur_vars,LKHD = LKHD_new,
                    dt_s = dt_s,dt_aou = dt_aou, m_s = m_s_new)
    colnames(temp_p) <- pool_vars
    #temp <- temp_p[, colSums(is.na(temp_p)) == 0] 
    
    temp_sign <- temp_p[,temp_p[3,] < alpha,drop = F]
    if(length(temp_sign)==0|all(is.na(temp_sign))){
      index_add <- 0
    }else{
      index_add <- which.max(temp_sign[4,])
      cur_vars <- c(cur_vars,names(index_add))
      LKHD_new <- temp_sign[1,index_add]
      m_s_new <- model.matrix(as.formula(paste0("~",paste0(cur_vars,collapse = "+"))), dt_s)
      pool_vars <- pool_vars[!names(index_add)==pool_vars]
    }
  }
  tryCatch(
    {
      temp <- fun.nps(as.formula(paste0("~",paste0(cur_vars,collapse = "+"))),dt,dt_aou,dt_s)
      d_vs <- temp$d # PS two-way
    },
    error = function(e)
      d_vs <<- NA
  )
  
  if(any(is.na(d_vs))){
    return(NA)
  }
  
  return(list(
    d = d,
    LKHD = LKHD,
    d_vs = d_vs,
    names_univar_new = cur_vars))
}


###############################
# 5) LIFO on GVS selected model
###############################


fun.lifo <- function(dt,dt_aou, dt_s,names_univar_new){
  out <- list()
  n_uni <- 0; n_max <- length(names_univar_new)
  clus_cal <- NULL
  n_stepwise <- rep(NA_real_, 3)
  while (is.null(clus_cal) && n_uni < n_max) {
    cat_temp <- as.formula(paste0("~", paste0(names_univar_new[seq_len(n_max-n_uni)], collapse = "+")))
    n_uni <- n_uni + 1
    temp_d <- fun.nps(cat_temp, dt,dt_aou, dt_s)$d
    dt_step <- cbind(dt_aou, weights = temp_d)
    clus <- survey::svydesign(id = ~1, weights = ~weights, data = dt_step)
    margin_names <- cal_names(cat_temp, clus)
    transformed_names <- sapply(margin_names, transform_term)
    poptotals <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
    clus_cal <- tryCatch(
      survey::calibrate(clus, cat_temp, poptotals, calfun = "raking", maxit = 1e2,
                        epsilon = rep(1e-3, length(poptotals)))
      ,error = function(e){
        return(NULL)
      }
    )
  }
  n_stepwise <- c(total = n_max, unselect = n_uni, select = n_max - n_uni)
  model_rails <- names_univar_new[seq_len(n_max-n_uni)]
  return(list(clus_cal = clus_cal,
              n_stepwise = n_stepwise,
              model_rails = model_rails))
}

###############################
# 6) Subgroup Raking
###############################
fun.subgroup.rake <- function(names_univar,names_var,dt,dt_aou,d){
  result <- list()
  out <- cbind(dt_aou, weights = d, ori = 1)
  out_sub <- eval(parse(text = paste0("split(out, out$", names_var, ")")))
  names_sub <- names(out_sub)
  names_mix_NR <- names_univar[!stringr::str_detect(names_univar, names_var)]
  index_nonconv <- numeric(length(names_sub))
  for (j in seq_along(names_sub)) {
    temp <- out_sub[[j]]
    clus <- survey::svydesign(id = ~1, weights = ~weights, data = temp)
    margin_names <- cal_names(as.formula(paste0("~", paste0(names_mix_NR, collapse = "+"))), clus)
    transformed_names <- sapply(margin_names, transform_term_sub, univar = names_var, j = names_sub[j])
    poptotals <- eval(parse(text = paste0("c(", paste(transformed_names, collapse = ", "), ")")))
    
    clus_cal <- tryCatch({
      survey::calibrate(clus, as.formula(paste0("~", paste0(names_mix_NR, collapse = "+"))),
                        poptotals, calfun = "raking", maxit = 1e2,
                        epsilon = rep(1e-7, length(poptotals)))
    },error = function(e){
      NULL
    })
    if(is.null(clus_cal)){
      index_nonconv[j] <- 1
    } else
      temp$ori <- weights(clus_cal)
    out_sub[[j]] <- temp
  } 
  
  out <- do.call(rbind, out_sub)
  return(list(out = out, index_nonconv = index_nonconv))
}

###############################
# 7) Data generation
###############################
sample.function <- function(seed = NULL,nsiz = 1e4,
                            z_prob = c(0.15, 0.35, 0.495, 0.005),
                            # Named parameter lists by x3 level ("1","2","3","4")
                            p1_list = c("1"=0.65, "2"=0.62, "3"=0.58, "4"=0.55),
                            p2_list = c("1"=0.65, "2"=0.62, "3"=0.58, "4"=0.55),
                            p3_list = list(
                              "1" = c(0.15, 0.35, 0.3, 0.2),
                              "2" = c(0.22, 0.37, 0.245, 0.165),
                              "3" = c(0.1, 0.24, 0.295, 0.365),
                              "4" = c(0.27, 0.18, 0.295, 0.255)
                            ),
                            p4_list = list(
                              "1" = c(0.00, 0.10),  # base, slope
                              "2" = c(0.02, 0.08),
                              "3" = c(0.04, 0.07),
                              "4" = c(0.05, 0.06)
                            )
) {
  if (!is.null(seed)) set.seed(seed)
  # alpha, beta, gamma are lists (or named vectors) controlling signals by z stratum
  # Returns a data.frame with factors catx1..catx4 and y, s, aou, ps, paou
  
  # x3c <- rtruncpareto(nsiz, lower = x3_lower, shape = x3_shape, upper = x3_upper)
  # qx3 <- quantile(x3c, probs = x3_cut_probs)
  # brks <- c(-Inf, qx3, Inf)
  # x3 <- as.integer(cut(x3c, breaks = brks, labels = FALSE))
  # catx3 <- factor(x3)
  # dt <- data.frame(id = seq_len(nsiz), x3, catx3, stringsAsFactors = FALSE)
  
  # Z as stratification variable
  z <- sample(1:4, nsiz, replace = TRUE, prob = z_prob)
  catz <- factor(z, levels = 1:4)
  dt <- data.frame(id = seq_len(nsiz), z, catz, stringsAsFactors = FALSE)
  
  # Pre-allocate
  x1 <- numeric(nsiz); catx1 <- integer(nsiz)
  x2 <- integer(nsiz)
  x3 <- integer(nsiz); catx3 <- integer(nsiz)
  x4 <- integer(nsiz); catx4 <- integer(nsiz)
  
  # 2) Generate by stratum (levels 1..4)
  for (k in 1:4) {
    idx <- which(z == k)
    n_k <- length(idx)
    if (n_k == 0) next
    
    key <- as.character(k)
    
    # x1 parameters for this stratum
    p1_k <- as.numeric(p1_list[[key]])
    if (length(p1_k) != 1 || is.na(p1_k)) stop("p1_list must give a single prob per level")
    x1_k <- rbinom(n_k, size = 1L, prob = p1_k) 
    
    # x2 ~ Bernoulli with level-specific prob
    p2_k <- as.numeric(p2_list[[key]])
    if (length(p2_k) != 1 || is.na(p2_k)) stop("p2_list must give a single prob per level")
    x2_k <- rbinom(n_k, size = 1L, prob = p2_k)
    
    # x3 ~ {0,1,2,3} with level-specific probs
    probs3 <- p3_list[[key]]
    if (length(probs3) != 4 || abs(sum(probs3) - 1) > 1e-8)
      stop("Each p3_list[[level]] must be length-4 and sum to 1.")
    x3_k <- sample(0:3, n_k, replace = TRUE, prob = probs3)  # 0,..,3
    
    # x4 ~ Bernoulli with prob depending on x5 and level
    p4_kv <- p4_list[[key]]
    if (is.list(p4_kv)) {
      base <- as.numeric(p4_kv$base); slope <- as.numeric(p4_kv$slope)
    } else {
      base <- as.numeric(p4_kv[1]); slope <- as.numeric(p4_kv[2])
    }
    if (anyNA(c(base, slope))) stop("p4_list[[level]] needs base and slope")
    prob_x4 <- pmin(0.99, pmax(0, base + slope * k))
    x4_k <- rbinom(n_k, size = 1L, prob = prob_x4)
    
    # Assign back
    x1[idx]     <- x1_k
    x2[idx]     <- x2_k
    x3[idx]     <- x3_k
    x4[idx]     <- x4_k
  }
  # Factors & standardization
  catx1 <- factor(x1, levels = 0:1)
  catx2 <- factor(x2, levels = 0:1)
  catx3 <- factor(x3, levels = 0:3)
  catx4 <- factor(x4, levels = 0:1)
  catx5 <- catz
  x1s <- as.numeric(scale(x1, center = FALSE, scale = TRUE))
  x2s <- as.numeric(scale(x2, center = FALSE, scale = TRUE))
  x3s <- as.numeric(scale(x3, center = TRUE, scale = TRUE))
  x4s <- as.numeric(scale(x4, center = FALSE, scale = TRUE))
  x5s <- as.numeric(scale(z, center = FALSE, scale = TRUE))
  # Output
  data.frame(
    x1, x2,x3,x4,x5 = z, 
    x1s, x2s,x3s,x4s,x5s,
    catx1,catx2,catx3,catx4,catx5,
    stringsAsFactors = FALSE
  )
}

outcome.function <- function(dt,alpha,beta,gamma){
  dt <- dt %>%
    group_split(catx5, .keep = TRUE) %>%
    map_dfr(function(subdf) {
      k <- as.character(unique(subdf$catx5))
      g <- alpha[[k]]
      pa = g[1] + g[2] * subdf$x1 + g[3]*subdf$x2 + g[4]*subdf$x3 + g[5]*subdf$x4+
                 g[6] * subdf$x1*subdf$x2 + g[7]*subdf$x2*subdf$x4
      mutate(subdf, py = expit(pa))
    })%>%
    mutate(y = rbinom(n(), size = 1, prob = py)) %>%
    mutate(ps = expit(beta[1] + beta[2]* x1 + beta[3]*x2  + beta[4]*x3 + beta[5]*x4+
                        beta[6] * x1*x2 + beta[7]*x2*x4
    )) %>%
    mutate(s = rbinom(n(), size = 1, prob = ps)) %>%
    mutate(fpc = n()) %>%
    group_split(catx5, .keep = TRUE) %>%
    map_dfr(function(subdf) {
      k <- as.character(unique(subdf$catx5))
      g <- gamma[[k]]
      pa = g[1] + g[2] * subdf$x1 + g[3]*subdf$x2 + g[4]*subdf$x3 + g[5]*subdf$x4+
        g[6] * subdf$x1*subdf$x2 + g[7]*subdf$x2*subdf$x4 
      return(mutate(subdf, paou = expit(pa)))
    }) %>%
    mutate(aou = rbinom(n(), size = 1, prob = paou))
  return(dt)
}

###############################
# 8) Output/diagnostic helpers
###############################

# Simple weighted mean & WR variance CI; "truth" is population mean truey.
fun.out <- function(mm, w, truey, targety, names_pop = NULL, svydesign = NULL) {
  m <- sum(w * mm) / sum(w)
  v <- sum(w^2 * (mm - m)^2) / (sum(w)^2)
  ci <- c(m - qnorm(0.975) * sqrt(v), m + qnorm(0.975) * sqrt(v))
  out <- numeric(16)
  out[1:4] <- c(m, v, ci)
  out[5] <- m - targety
  out[6] <- out[5]^2
  out[7] <- as.numeric(truey >= ci[1] & truey <= ci[2])
  out[8] <- as.numeric(targety >= ci[1] & targety <= ci[2])
  out[9] <- sum(w)
  out[10] <- var(w)
  out[11] <- mean(w <= 0)
  if (!is.null(svydesign)) {
    est <- survey::svymean(~ y, svydesign)
    out[12] <- as.numeric(vcov(est))
    ci2 <- confint(est)
    out[13:14] <- ci2
    out[15] <- as.numeric(truey >= ci2[1] & truey <= ci2[2])
    out[16] <- as.numeric(object.size(svydesign))
  } else out[12:16] <- NA_real_
  if (!is.null(names_pop)) names(out) <- names_pop
  out
}

# Survey-only convenience
svymean.out <- function(svydesign, truey) {
  est <- survey::svymean(~ y, svydesign)
  m <- as.numeric(coef(est))
  v <- as.numeric(vcov(est))
  ci <- confint(est)
  w  <- survey::weights(svydesign)
  c(m, v, ci[1], ci[2], m - truey, (m - truey)^2, as.numeric(truey >= ci[1] & truey <= ci[2]),
    sum(w), var(w), mean(w <= 0))
}

# Subgroup evaluation: columns of G are subgroup indicators (0/1)
# Returns for each subgroup: size, mean, variance, LB, UB, Coverage (vs pre_pop)
fun.pre0 <- function(G, w, y, pre_pop) {
  stopifnot(nrow(G) == length(w), length(w) == length(y))
  K <- ncol(G)
  res <- matrix(NA_real_, nrow = 6, ncol = K)
  for (k in seq_len(K)) {
    idx <- as.logical(G[,k])
    nk <- sum(idx)
    if (nk == 0) next
    wk <- w[idx]
    yk <- y[idx]
    m  <- sum(wk * yk) / sum(wk)
    v  <- sum(wk^2 * (yk - m)^2) / (sum(wk)^2)
    ci <- c(m - qnorm(0.975) * sqrt(v), m + qnorm(0.975) * sqrt(v))
    res[,k] <- c(nk, m, v, ci, as.numeric(pre_pop[k] >= ci[1] & pre_pop[k] <= ci[2]))
  }
  rownames(res) <- c("size","mean","variance","LB","UB","Coverage")
  as.data.frame(res)
}

# Same, but also returns overall pre_all-matching check vector for EE
fun.pre <- function(G, w, y, pre_pop) {
  fun.pre0(G, w, y, pre_pop)
}

# Estimating equation residual vector: one per margin column
# If pre_all are means, multiply by sum(w) to convert to totals
fun.ee <- function(G, w, y, pre_all) {
  y <- as.numeric(y)
  # moment residuals per column
  left  <- as.numeric(t(G) %*% (w * y))
  right <- as.numeric(sum(w) * pre_all)
  left - right
}



