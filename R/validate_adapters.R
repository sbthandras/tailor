#' Validate adapter data frame
#' 
#' This function checks if the input data frame has the required structure and 
#' variable types to be used as an adapter data frame. It ensures that the 
#' necessary columns are present and of the correct type, and it assigns the 
#' "adapter" class to the data frame if it passes validation.
#' @param adapters a data frame that is expected to contain adapter information.
#' @return A data frame of class "adapter" if the input is valid.,If the input 
#' is not valid, the function will stop and return an error message indicating 
#' the issue.
#' @noRd
validate_adapters <- function(adapters) {
  if (inherits(adapters, "adapter")) {
    return(adapters)
  }
  if (!inherits(adapters, "data.frame")) {
    stop("'adapters' must be a data frame.")
  }
  required_cols <- c(
    "pattern_id", "subject_id", "start", "end", "mean_score", "pident"
  )
  missing_cols <- setdiff(required_cols, colnames(adapters))
  if (length(missing_cols) > 0) {
    stop(
      "The following required variables are missing from 'adapters': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.character(adapters$pattern_id)) {
    stop("'pattern_id' variable must be of type character.")
  }
  if (!is.character(adapters$subject_id)) {
    stop("'subject_id' variable must be of type character.")
  }
  if (!is.numeric(adapters$start)) {
    stop("'start' variable must be of type numeric.")
  }
  if (!is.numeric(adapters$end)) {
    stop("'end' variable must be of type numeric.")
  }
  if (!is.numeric(adapters$pident)) {
    stop("'pident' variable must be of type numeric.")
  }
  if (!is.numeric(adapters$mean_score)) {
    stop("'mean_score' variable must be of type numeric.")
  }
  # remove duplicate rows
  adapters <- adapters |> dplyr::distinct()
  # order pairs alphabetically and drop reverse duplicates with identical values
  adapters <- adapters |>
    dplyr::mutate(
      .id1 = pmin(pattern_id, subject_id),
      .id2 = pmax(pattern_id, subject_id),
      .rev = pattern_id > subject_id
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      pattern_id = ifelse(.rev, .id2, .id1),
      subject_id = ifelse(.rev, .id1, .id2)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.id1, -.id2, -.rev) |>
    dplyr::distinct(
      pattern_id, subject_id, start, end, mean_score, pident, .keep_all = TRUE
    )
  # check that pair appears more than once
  pairs <- adapters |>
    dplyr::select(pattern_id, subject_id) |>
    dplyr::distinct()
  if (nrow(pairs) != nrow(adapters)) {
    stop("Pairs must be unique.")
  }
  class(adapters) <- c("adapter", class(adapters))
  return(adapters)
}