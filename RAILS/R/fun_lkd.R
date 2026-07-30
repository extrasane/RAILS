#' Forward likelihood-ratio test for one candidate interaction
#'
#' Fits the pseudo-likelihood propensity model with one candidate interaction
#' term added and returns a likelihood-ratio test of that term against the
#' current model. Used by [fun.rails.threeway()] for forward selection.
#'
#' @param x Candidate term, e.g. `"agegroup:sex:income"`.
#' @param names_univar Character vector of terms already in the current model.
#' @param LKHD Log-likelihood of the current model.
#' @param dt_s Aggregated reference-sample cell counts (column `weight`).
#' @param dt_aou Aggregated non-probability cell counts (column `weight`).
#' @param m_s Model matrix for the current model on `dt_s` (used only for its
#'   column count, to compute degrees of freedom).
#'
#' @return A length-5 numeric vector
#'   `c(log-likelihood, LRT, p-value, LRT/df, df)`, or `rep(NA, 5)` if any
#'   cell count is zero or the solver fails.
#'
#' @export
fun.lkd <- function(x, names_univar, LKHD, dt_s, dt_aou, m_s) {
  cat_formula <- stats::formula(paste0("~", paste0(paste0(names_univar, collapse = "+"), "+", x)))
  t_s   <- Matrix::sparse.model.matrix(cat_formula, dt_s)
  t_aou <- Matrix::sparse.model.matrix(cat_formula, dt_aou)
  theta <- rep(0, ncol(t_aou))
  pia   <- rep(1 / 2, nrow(t_s))
  W_s   <- dt_s$weight
  W_aou <- dt_aou$weight
  U1    <- as.numeric(Matrix::crossprod(t_aou, W_aou))

  if (any(U1 == 0)) return(rep(NA, 5))

  res       <- 1
  temp_LKHD <- NA
  iter      <- 0

  tryCatch({
    while (res >= 1e-3) {
      iter <- iter + 1
      if (iter > 1000) stop("non-convergence")
      W1      <- W_s * pia
      Uscore  <- U1 - as.numeric(Matrix::crossprod(t_s, W1))
      Hsov    <- solve(as.matrix(Matrix::crossprod(t_s, t_s * (W1 * (1 - pia)))))
      theta0  <- as.numeric(theta + Hsov %*% Uscore)
      pia     <- stats::plogis(as.numeric(t_s %*% theta0))
      res     <- sum((Hsov %*% Uscore)^2)
      theta   <- theta0
      temp_LKHD <- sum(U1 * theta) - sum(W_s * log1p(exp(as.numeric(t_s %*% theta))))
    }
    df   <- ncol(t_s) - ncol(m_s)
    LRT  <- 2 * (temp_LKHD - LKHD)
    pval <- 1 - stats::pchisq(LRT, df = df)
    c(temp_LKHD, LRT, pval, LRT / df, df)
  }, error = function(e) rep(NA, 5))
}
