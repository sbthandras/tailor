#' Validate RBP data set
#'
#' This function can be used to check the RBP data set for any consistency
#' issues that would produce errors in downstream analysis.
#' @param rbps data.frame
#' @param id_var name of the variable in \code{rbps} that contains unique
#' identifiers for each RBP. Must be unique and of class "character".
#' @param seq_var name of the variable in \code{rbps} that contains the amino
#' acid sequences. Must be of class "character" and cannot contain empty strings
#' or unknown characters.
#' @param verbose should the function return verbose messages to console?
#' @return The function returns \code{TRUE} if the tails data set is consistent
#' and returns \code{FALSE} if it is not.
#' @details The function checks whether id_var is present in the data set and is
#' unique, and whether seq_var is present in the data set, is of class
#' "character", and does not contain empty strings. The function also checks
#' whether seq_var contains unknown characters, i.e. characters other than the
#' 20 standard amino acids and the ambiguous characters B, Z, and X. If any of
#' these checks fail, the function returns \code{FALSE} and, if
#' \code{verbose = TRUE}, also returns a warning message to console.
#' @examples
#' data(rbps)
#' validate_rbps(rbps, verbose = TRUE)
#' @export
validate_rbps <- function(
    rbps,
    id_var = "Core_ORF",
    seq_var = "translation",
    verbose = getOption("verbose")
  ) {
  if (!inherits(rbps, "data.frame")) {
    stop("Argument 'rbps' must be a data.frame.")
  }
  if (id_var %in% names(rbps) == FALSE) {
    if (verbose) warning(paste0("Variable '", id_var, "' not found."))
    return(FALSE)
  }
  if (!inherits(rbps[[id_var]], "character")) {
    if (verbose) warning(paste0("Variable '", id_var, "' is not character."))
    return(FALSE)
  }
  if (any(is.na(rbps[[id_var]]))) {
    if (verbose) warning(paste0("Variable '", id_var, "' contains NA values."))
    return(FALSE)
  }
  if (length(unique(rbps[[id_var]])) != nrow(rbps)) {
    if (verbose) warning(paste0("Variable '", id_var, "' not unique."))
    return(FALSE)
  }
  if (seq_var %in% names(rbps) == FALSE) {
    if (verbose) warning(paste0("Variable '", seq_var, "' not found."))
    return(FALSE)
  }
  if (!inherits(rbps[[seq_var]], "character")) {
    if (verbose) warning(paste0("Variable '", seq_var, "' is not character."))
    return(FALSE)
  }
  if (any(is.na(rbps[[seq_var]]))) {
    if (verbose) warning(paste0("Variable '", seq_var, "' contains NA values."))
    return(FALSE)
  }
  if (0 %in% nchar(rbps[[seq_var]])) {
    index <- which(nchar(rbps[[seq_var]]) == 0)
    if (verbose) {
      index <- paste(index, collapse = ", ")
      warning(paste0("Variable '", seq_var, "' contains an empty string in the ",
                     "following rows: ", index))
    }
    return(FALSE)
  }
  AA <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V", "B", "Z", "X")
  out <- sapply(seq_along(rbps[[seq_var]]), function(x) {
    mean(strsplit(rbps[[seq_var]][x], "")[[1]] %in% AA)
  })
  index <- which(out < 1)
  if (length(index) > 0) {
    if (verbose) {
      index <- paste(index, collapse = ", ")
      warning(paste0("Variable '", seq_var, "' contains unknown character(s) in ",
                     "the following rows: ", index))
    }
    return(FALSE)
  }
  return(TRUE)
}
