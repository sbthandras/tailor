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
#' @param value character; the variable within the adapter data frame to use as
#' values in the matrix, "pident", "mean_score", or "end".
#' @param cores integer; number of CPU cores to use.
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
    cores = 1
  ) {
  if (!inherits(adapters, "adapter")) {
    msg <- paste0(
      "'adapters' must be a data frame of class 'adapter'. ",
      "Use find_adapter() to create a compatible data frame."
    )
    stop(msg)
  }
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
  value <- match.arg(value, c("pident", "mean_score", "end"))
  cl <- parallel::makeCluster(cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  dmat <- matrix(NA, nrow = length(ids), ncol = length(ids))
  rownames(dmat) = ids
  colnames(dmat) = ids
  dmat <- foreach::foreach(
    i = seq_len(nrow(dmat)), .combine = 'rbind',.inorder=TRUE) %dopar% {
      i <- i
      for (j in seq_len(ncol(dmat))){
        if (i == j) next()
        adapter <- adapters |>
          filter_comparisons(ids[i]) |>
          filter_comparisons(ids[j])
        if (nrow(adapter) == 0) {
          stop("No comparisons found for ids: ", ids[i], ", ", ids[j], ".")
        } else if (nrow(adapter) == 1) {
          dmat[i, j] <- ifelse(is.na(adapter[[value]]), 0, adapter[[value]])
        } else {
          stop("Multiple comparisons found for ids: ", ids[i], ", ", ids[j], ".")
        }
      }
      return(dmat[i,])
    }
  if (value == "pident") diag(dmat) <- 1
  rownames(dmat) <- ids
  colnames(dmat) <- ids
  class(dmat) <- c("adapter_matrix", class(dmat))
  return(dmat)
}
