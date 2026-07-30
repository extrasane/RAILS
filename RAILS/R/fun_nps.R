#' Pseudo-likelihood propensity score weights
#'
#' Estimates non-probability-sample (NPS) weights by maximizing the
#' pseudo-likelihood propensity score model via Newton-Raphson on aggregated
#' covariate cells. The returned weights are scaled to sum to `nsiz`.
#'
#' @param cat_temp A model formula, e.g. `~ agegroup + sex + agegroup:sex`.
#' @param dt_aou Aggregated non-probability (AoU) cell counts, with a column
#'   `weight` giving the cell size.
#' @param dt_s Aggregated reference-sample cell counts, with a column
#'   `weight`.
#' @param nsiz Target population total the weights are scaled to. Defaults to
#'   the total reference-sample weight, `sum(dt_s$weight)`.
#' @param theta_init Optional warm-start coefficient vector (e.g. the previous
#'   step's solution, zero-padded for new columns). The pseudo-likelihood is
#'   globally concave, so the warm start only reduces the iteration count.
#'
#' @return A numeric vector of NPS weights (one per cell), scaled to sum to
#'   `nsiz`, with the converged coefficient vector attached as attribute
#'   `"theta"`.
#'
#' @details The Newton-Raphson loop is capped at 1000 iterations and errors if
#'   it fails to converge.
#'
#' @export
fun.nps <- function(cat_temp, dt_aou, dt_s, nsiz = sum(dt_s$weight), theta_init = NULL) {
  m_aou  <- Matrix::sparse.model.matrix(cat_temp, dt_aou)
  m_nhis <- Matrix::sparse.model.matrix(cat_temp, dt_s)
  theta  <- if (is.null(theta_init)) rep(0, ncol(m_aou)) else
    c(theta_init, rep(0, ncol(m_aou) - length(theta_init)))
  pia    <- stats::plogis(as.numeric(m_nhis %*% theta))
  W_s    <- dt_s$weight
  W_aou  <- dt_aou$weight
  U1     <- as.numeric(Matrix::crossprod(m_aou, W_aou))

  res  <- 1
  iter <- 0
  while (res >= 1e-10) {
    iter <- iter + 1
    if (iter > 1000) stop("fun.nps: no convergence within 1000 iterations")
    W1     <- W_s * pia
    Uscore <- U1 - as.numeric(Matrix::crossprod(m_nhis, W1))
    Hsov   <- solve(as.matrix(Matrix::crossprod(m_nhis, m_nhis * (W1 * (1 - pia)))))
    theta0 <- as.numeric(theta + Hsov %*% Uscore)
    pia    <- stats::plogis(as.numeric(m_nhis %*% theta0))
    res    <- sum((Hsov %*% Uscore)^2)
    theta  <- theta0
  }

  d <- 1 / stats::plogis(as.numeric(m_aou %*% theta))
  d <- d / sum(d) * nsiz
  attr(d, "theta") <- theta
  d
}
