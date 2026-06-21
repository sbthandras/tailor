#' Create adapter matrix
#'
#' Convert identified adapters into matrix form.
#' @importFrom dplyr distinct
#' @param adapters adapter; a data frame of class "adapter" containing
#' identified adapter sequences.
#' @param ids character; ids within the 'adapters' data frame, in the same
#' order as they should appear in the matrix. If NULL (default), order will be
#' determined by the order of unique ids in the 'adapters' data frame. If
#' specified, all ids must be present in the'adapters' data frame and all ids
#' in the 'adapters' data frame must be present in the 'ids' vector.
#' @param value character; a numeric variable within the adapter data frame to 
#' use as matrix values.
#' @param verbose logical; should verbose messages be printed to the console?
#' If TRUE, a progress bar is also displayed.
#' @return A symmetric similarity matrix. If value = "pident", the diagonal is 
#' set to 1. The matrix is of class "adapter_matrix" and can be used with the 
#' plot() function to visualize the similarity between sequences based on shared 
#' adapters.
#' @note The function handles incomplete adapter data frames where some sequence
#' pairs are missing. Missing pairs are assumed to lack shared adapters and are
#' filled with 0 values in the matrix. However, sequences with no shared
#' adapters across any pair are not recoverable if they were dropped from the
#' data frame.
#' @examples
#' # import example data
#' data(rbps)
#' # import adapters identified from example data
#' data(adapters)
#' # convert to matrix and plot
#' adapters |> adapter_matrix() |> plot()
#' # use same order as in the original data frame
#' adapters |> adapter_matrix(ids = rbps$Core_ORF) |> plot()
#' @importFrom foreach %dopar%
#' @export
adapter_matrix <- function(
    adapters,
    ids = NULL,
    value = "pident",
    verbose = getOption("verbose")
  ) {
  adapters <- validate_adapters(adapters)
  if (is.null(ids)) {
    ids <- unique(c(adapters$pattern_id, adapters$subject_id))
  } else {
    A <- setdiff(ids, unique(c(adapters$pattern_id, adapters$subject_id)))
    B <- setdiff(unique(c(adapters$pattern_id, adapters$subject_id)), ids)
    if( length(A) > 0) {
      msg <- paste0(
        "The following ids were not found in the 'adapter' data frame: ",
        paste(A, collapse = ", ")
      )
      stop(msg)
    }
    if (length(B) > 0) {
      msg <- paste0(
        "The following ids were not found in the 'ids' vector: ",
        paste(B, collapse = ", ")
      )
      stop(msg)
    }
  }
  valid_values <- names(adapters)[
    !names(adapters) %in% c("pattern_id", "subject_id") &
    sapply(adapters, is.numeric)
  ]
  value <- match.arg(value, choices = valid_values)
  dmat <- matrix(0, nrow = length(ids), ncol = length(ids))
  rownames(dmat) = ids
  colnames(dmat) = ids
  if (verbose) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = nrow(adapters),
      style = 3
    )
  }
  for (i in seq_len(nrow(adapters))) {
    index_x <- which(ids == adapters$pattern_id[i])
    index_y <- which(ids == adapters$subject_id[i])
    hit <- ifelse(is.na(adapters[[value]][i]), 0, adapters[[value]][i])
    if (length(hit) > 1) {
      stop(
        "Multiple values found for ids: ",
        adapters$pattern_id[i],
        ", ", adapters$subject_id[i], "."
      )
    }
    dmat[index_x, index_y] <- hit
    dmat[index_y, index_x] <- hit
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)
  if (value == "pident") diag(dmat) <- 1
  rownames(dmat) <- ids
  colnames(dmat) <- ids
  class(dmat) <- c("adapter_matrix", class(dmat))
  return(dmat)
}
