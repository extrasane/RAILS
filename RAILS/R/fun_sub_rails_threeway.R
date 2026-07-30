#' Subgroup RAILS estimator
#'
#' Runs [fun.rails.threeway()] separately within each level of a subgroup
#' variable. For each level it subsets both aggregated cell tables, builds the
#' subgroup's own population margins from the reference cell subset, scales
#' weights to that subgroup's total, and tags the result with `subgroup_run`.
#'
#' @param dt_agg_aou Aggregated non-probability cell counts. `subgroup_var`
#'   must be one of its aggregation dimensions (a column).
#' @param dt_agg_pums Aggregated reference-sample cell counts, likewise
#'   containing `subgroup_var`.
#' @param subgroup_var Name of the stratifying variable. Automatically excluded
#'   from `names_univar`.
#' @param names_univar Main-effect variables within each subgroup. Defaults to
#'   the seven shared covariates minus `subgroup_var`.
#' @param alpha Significance threshold for forward LRT selection, passed to
#'   [fun.rails.threeway()].
#'
#' @return A data frame of the row-bound per-subgroup results: every column of
#'   [fun.rails.threeway()]'s output plus `subgroup_run`. Selection and the
#'   LIFO walk run independently within each level, so `selected_terms` and
#'   `calibrated_terms` may differ across subgroups.
#'
#' @details Factor levels are preserved (no `droplevels()`) so that model
#'   matrix columns match between the AoU and reference subsets. Levels with no
#'   reference cells are skipped with a warning.
#'
#' @seealso [fun.rails.threeway()].
#'
#' @export
fun.sub.rails.threeway <- function(
    dt_agg_aou,
    dt_agg_pums,
    subgroup_var = "region",
    names_univar = setdiff(
      c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
      subgroup_var
    ),
    alpha        = 0.05
) {
  if (subgroup_var %in% names_univar) {
    stop("fun.sub.rails.threeway: subgroup_var must not appear in names_univar")
  }
  if (!subgroup_var %in% names(dt_agg_aou) || !subgroup_var %in% names(dt_agg_pums)) {
    stop("fun.sub.rails.threeway: subgroup_var must be a column of both cell tables")
  }

  ## Max cross-tabulation formula for the subgroup-level population margins
  twovars     <- utils::combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
  threevars   <- utils::combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))
  max_formula <- stats::formula(paste0("~", paste(c(names_univar, twovars, threevars), collapse = "+")))

  subgroup_levels <- sort(unique(dt_agg_aou[[subgroup_var]]))

  out <- lapply(subgroup_levels, function(lev) {

    message("=== Subgroup ", subgroup_var, " = ", lev, " ===")

    agg_aou_sub  <- dt_agg_aou  %>% filter(.data[[subgroup_var]] == lev)
    agg_pums_sub <- dt_agg_pums %>% filter(.data[[subgroup_var]] == lev)

    if (nrow(agg_pums_sub) == 0) {
      warning("fun.sub.rails.threeway: no reference cells for ", subgroup_var,
              " = ", lev, " — skipped.")
      return(NULL)
    }

    ## Subgroup population margins, computed once from the reference cell subset
    mat_max    <- Matrix::sparse.model.matrix(max_formula, data = agg_pums_sub, keep.order = TRUE)
    pop_totals <- stats::setNames(
      as.numeric(Matrix::crossprod(mat_max, agg_pums_sub$weight)),
      colnames(mat_max)
    )

    fun.rails.threeway(
      dt_agg_aou   = agg_aou_sub,
      dt_agg_pums  = agg_pums_sub,
      pop_totals   = pop_totals,
      names_univar = names_univar,
      alpha        = alpha,
      nsiz         = sum(agg_pums_sub$weight)
    ) %>%
      mutate(subgroup_run = lev)
  })

  bind_rows(out)
}
