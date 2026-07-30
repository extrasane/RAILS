#' @keywords internal
"_PACKAGE"

## Package-level imports. External calls are namespace-qualified in each
## function (Matrix::, survey::, stats::, utils::); only the pipe, the .data
## pronoun, and the dplyr verbs used inside pipes are imported here.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter bind_rows
NULL
