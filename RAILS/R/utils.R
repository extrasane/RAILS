## Internal helpers — not exported.

#' Standardize an interaction term name
#'
#' Sorts the colon-separated components of a term alphabetically so that,
#' e.g., `"sex:region"` and `"region:sex"` compare equal.
#'
#' @param x A single character string, possibly containing `":"`.
#' @return The term with components sorted, rejoined by `":"`.
#' @keywords internal
#' @noRd
standardize_string <- function(x) {
  paste(sort(unlist(strsplit(x, ":"))), collapse = ":")
}

#' Match calibration targets to design columns
#'
#' Subsets a named vector of population margins (`v1`) to the calibration
#' column names of a survey design (`v2`), matching by standardized term name
#' and returning the values in `v2` order (which is the order
#' [survey::calibrate()] expects).
#'
#' @param v1 Named numeric vector of population margins.
#' @param v2 Character vector of calibration column names.
#' @return Numeric vector named by `v2`, containing the matched `v1` values.
#' @keywords internal
#' @noRd
create_v3 <- function(v1, v2) {
  v1_std <- vapply(names(v1), standardize_string, character(1), USE.NAMES = FALSE)
  v2_std <- vapply(v2,        standardize_string, character(1), USE.NAMES = FALSE)
  idx  <- match(v2_std, v1_std)
  keep <- !is.na(idx)
  stats::setNames(as.numeric(v1[idx[keep]]), v2[keep])
}
