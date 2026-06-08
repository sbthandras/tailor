#' Find all adapters for a set of RBPs
#'
#' Look at all pairs between a set of RBPs, align, detect breakpoints, find
#' adapters and combine them into a single data frame.
#' @param ids character; the names of RBPs to compare.
#' @param data data.frame; a data frame that contains IDs and sequences
#' @param id_var character; variable within \code{data} that stores IDs.
#' @param seq_var character; variable within \code{data} that stores sequences.
#' @param submat character; the substitution matrix. See
#' \code{position_scores()} for more information.
#' @param max_end_vars character; optional vector of variable names in
#' \code{data} that contain numeric values representing positions which the
#' adapters should not overlap with and terminate earlier (e.g., domain
#' start positions). See Details for more information.
#' @param max_end integer; the maximum position to look for breakpoints. If
#' `NULL`, search the full alignment, if an integer, limit the search to that
#' position. If both this argument and `position_scores$max_end` are set, the
#' smaller integer is used.
#' @param method character; the method to use for finding breakpoints. Either
#' "cemean", "plateau", or "window". See \code{find_breakpoints()} for more
#' information
#' @param min_pident numeric; the lowest accepted percentage identity for a
#' region to be considered a "conserved" region.
#' @param max_start integer; the maximum accepted starting position for the
#' first conserved region.
#' @param merge_beginning logical; whether to merge the positions before the
#' first conserved region. See Details for more information.
#' @param merge_conserved logical; whether neighboring conserved regions should
#' be merged. See Details for more information.
#' @param cores integer; number of CPU cores to use.
#' @param verbose logical; should verbose messages be printed to the console?
#' @param ... additional arguments to pass to breakpoint detection. See
#' \code{find_breakpoints()} for more information.
#' @return A data.frame of class "adapter".
#' @details If you provide one or more `max_end_var` names, the function will
#' look for their values in both the pattern and subject rows, calculate their
#' minimum  and return it in the output object. The value will then be
#' respected by `find_breakpoints()` and the function will not look at
#' positions beyond this value.
#' @note When comparing all pairs of sequences, most pairs will not share an
#' adapter, resulting in a sparse data frame. It makes sense to remove NA rows
#' and export only those with shared adapters. The `adapter_matrix()` function
#' assumes missing pairs did not have a shared adaptor and restores these pairs
#' with 0 values in the matrix. However, sequences with no shared adapters to
#' any other sequence are dropped from the data frame and cannot be recovered
#' in the matrix.
#' @examples
#' data(rbps)
#' find_all_adapters(rbps$Core_ORF[1:3], data = rbps)
#' @export
find_all_adapters <- function(
    ids,
    data,
    id_var = "Core_ORF",
    seq_var = "translation",
    submat = "BLOSUM80",
    max_end_vars = NULL,
    max_end = NULL,
    method = "cemean",
    min_pident = 0.4,
    max_start = 10,
    merge_beginning = TRUE,
    merge_conserved = TRUE,
    cores = 1,
    verbose = getOption("verbose"),
    ...
  ) {

  pairs <- data |>
    dplyr::filter(!!rlang::sym(id_var) %in% ids) |>
    dplyr::pull(!!rlang::sym(id_var)) |>
    utils::combn(2) |>
    t() |>
    as.data.frame()
  names(pairs) <- c("pattern", "subject")
  adapters <- data.frame()

  if (cores == 1) {
    # Use regular for loop when cores = 1
    adapters_list <- list()
    for (i in seq_len(nrow(pairs))) {
      pattern <- pairs$pattern[i]
      subject <- pairs$subject[i]
      ps <- position_scores(
        pattern_id = pattern,
        subject_id = subject,
        id_var = id_var,
        seq_var = seq_var,
        submat = submat,
        max_end_vars = max_end_vars,
        data = data,
        verbose = verbose
      )
      adapter <- ps |>
        find_breakpoints(max_end = max_end, method = method, ...) |>
        find_adapter(
          min_pident = min_pident,
          max_start = max_start,
          merge_beginning = merge_beginning,
          merge_conserved = merge_conserved
        )
      adapters_list[[i]] <- adapter
    }
    adapters <- do.call(rbind, adapters_list)
  } else {
    # Use parallelization when cores > 1
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    adapters <- foreach::foreach(
      i = seq_len(nrow(pairs)), .combine = 'rbind'
      ) %dopar% {
        i <- i
        pattern <- pairs$pattern[i]
        subject <- pairs$subject[i]
        ps <- position_scores(
          pattern_id = pattern,
          subject_id = subject,
          id_var = id_var,
          seq_var = seq_var,
          submat = submat,
          max_end_vars = max_end_vars,
          data = data,
          verbose = verbose
        )
        adapter <- ps |>
          find_breakpoints(max_end = max_end, method = method, ...) |>
          find_adapter(
            min_pident = min_pident,
            max_start = max_start,
            merge_beginning = merge_beginning,
            merge_conserved = merge_conserved
          )
        return(adapter)
    }
  }
  return(adapters)
}
