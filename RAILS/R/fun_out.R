#' Weight diagnostics
#'
#' Summarizes a vector of weights for quality checking after estimation.
#'
#' @param w Numeric vector of weights.
#' @param ... Ignored; present for call compatibility.
#'
#' @return A named numeric vector with elements `sum`, `var`,
#'   `non-positive ratio`, `Less than 1`, `min`, and `max`.
#'
#' @examples
#' fun.out(c(1.2, 0.8, 3.4, 0.5))
#'
#' @export
fun.out <- function(w, ...) {
  temp <- numeric(6)
  temp[1] <- sum(w, na.rm = TRUE)
  temp[2] <- stats::var(w, na.rm = TRUE)
  temp[3] <- mean(w <= 0, na.rm = TRUE)
  temp[4] <- sum(w < 1)
  temp[5] <- min(w, na.rm = TRUE)
  temp[6] <- max(w, na.rm = TRUE)
  names(temp) <- c("sum", "var", "non-positive ratio", "Less than 1", "min", "max")
  temp
}
