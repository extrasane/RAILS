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
  m_aou <- model.matrix(cat_temp , dt_aou)
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
  d <-  1/pi_aou
  d <- d/sum(d)*nsiz
  return(list(d = d, LKHD = LKHD))
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
    return(NA)
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
