#' Group RBPs into adapter clusters
#'
#' @param mat adapter_matrix; a matrix of class "adapter_matrix". Use 
#' adapter_matrix() to create a compatible matrix.
#' @param k_min integer; minimum number of clusters to split the RBPs into.
#' @param k_max integer; maximum number of clusters to split the RBPs into.
#' @param h hclust object; a hierarchical clustering object created from the 
#' distance matrix. If NULL, the function will create a hierarchical clustering 
#' object from the input matrix using the average method.
#' @param homogeneity_thr numeric; minimum percentage identity of the identified
#' adapter sequence to be considered a shared adapter for homogeneity.
#' @param completeness_thr numeric; minimum percentage identity of the
#' identified adapter sequence to be considered a shared adapter for
#' completeness.
#' @param cores integer; number of CPU cores to use.
#' @param verbose logical; should verbose messages be printed to the console?
#' If TRUE, a progress bar is also displayed.
#' @return a data frame with two columns, id and cluster, containing cluster
#' assignments.
#' @examples
#' data(rbps)
#' data(adapters)
#' # convert to pident matrix and order by the original data frame
#' amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
#' # cluster RBPs into 2 to 5 clusters
#' cluster_adapters(amat, k_min = 2, k_max = 5)
#' @importFrom foreach %dopar%
#' @export
cluster_adapters <- function(
    mat,
    k_min = 1,
    k_max = NULL,
    h = NULL,
    homogeneity_thr = 0.75,
    completeness_thr = 0.75,
    cores = 1,
    verbose = getOption("verbose")
  ) {
  if (!inherits(mat, "adapter_matrix")) {
    msg <- paste0(
      "mat must be a matrix of class 'adapter_matrix'. ",
      "Use adapter_matrix() to create a compatible matrix."
    )
    stop(msg)
  }
  if (is.null(h)) {
    if (verbose) message("Running hierarchical clustering.")
    h <- mat |> stats::dist() |> stats::hclust(method = "average")
  }
  if (is.null(k_max)) k_max <- nrow(mat)
  if (verbose) message("Optimising number of clusters.")
  if (cores == 1) {
    if (verbose) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = length(k_min:k_max),
        style = 3
      )
    }
    schocos_list <- list()
    for (k in k_min:k_max) {
      schocos_list[[length(schocos_list) + 1]] <- schoco(
        mat, k, h, homogeneity_thr, completeness_thr, verbose = FALSE
      )
      if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)
    schocos <- do.call(rbind, schocos_list)
  } else {
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    schocos <- foreach::foreach(k = k_min:k_max, .combine = 'rbind') %dopar% {
      schoco(mat, k, h, homogeneity_thr, completeness_thr, verbose = FALSE)
    }
  }

  schocos <- split(schocos, schocos$nclust) |> unname()

  schoco_sums <- sapply(schocos, function(x) sum(x$schoco, na.rm = TRUE))
  index <- min(which(schoco_sums == max(schoco_sums)))

  k <- schocos[[index]]$nclust |> unique()

  cluster <- paste0("ACL ", stats::cutree(tree = h, k = k))

  out <- data.frame(
    id = row.names(mat),
    cluster = cluster
  )

  return(out)
}
