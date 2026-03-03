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
#' @param method character; the method to use for finding breakpoints. Either
#' "cemean", "plateau", or "window". See \code{find_breakpoints()} for more
#' information
#' @param pident_threshold numeric; the lowest accepted percentage identity for
#' a region to be considered a "conserved" region.
#' @param start_threshold numeric; the highest accepted starting position for
#' conserved region to be considered an "N-terminal" region, expressed as the
#' ratio of the aligned sequence length.
#' @param merge logical; whether neighboring conserved regions should be merged.
#' @param cores integer; number of CPU cores to use.
#' @param verbose logical; should verbose messages be printed to the console?
#' @param ... additional arguments to pass to breakpoint detection. See
#' \code{find_breakpoints()} for more information.
#' @return A data.frame of class "adapter".
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
    method = "cemean",
    pident_threshold = 0.4,
    start_threshold = 0.25,
    merge = TRUE,
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
        data = data,
        verbose = verbose
      )
      adapter <- ps |>
        find_breakpoints(method = method, ...) |>
        find_adapter(
          pident_threshold = pident_threshold,
          start_threshold = start_threshold,
          merge = merge
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
          data = data,
          verbose = verbose
        )
        adapter <- ps |>
          find_breakpoints(method = method, ...) |>
          find_adapter(
            pident_threshold = pident_threshold,
            start_threshold = start_threshold,
            merge = merge
          )
        return(adapter)
    }
  }
  return(adapters)
}
