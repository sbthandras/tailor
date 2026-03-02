#' Group RBPs into adapter clusters
#'
#' @param mat adapter_matrix; a matrix of class "adapter_matrix". Use 
#' adapter_matrix() to create a compatible matrix.
#' @param k_min integer; minimum number of clusters to split the RBPs into.
#' @param k_max integer; maximum number of clusters to split the RBPs into.
#' @param homogeneity_thr numeric; minimum percentage identity of the identified
#' adapter sequence to be considered a shared adapter for homogeneity.
#' @param completeness_thr numeric; minimum percentage identity of the
#' identified adapter sequence to be considered a shared adapter for
#' completeness.
#' @param cores integer; number of CPU cores to use.
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
    homogeneity_thr = 0.75,
    completeness_thr = 0.75,
    cores = 1
  ) {
  if (!inherits(mat, "adapter_matrix")) {
    msg <- paste0(
      "mat must be a matrix of class 'adapter_matrix'. ",
      "Use adapter_matrix() to create a compatible matrix."
    )
    stop(msg)
  }
  h <- mat |> stats::dist() |> stats::hclust(method = "average")

  if (is.null(k_max)) k_max <- nrow(mat)

  cl <- parallel::makeCluster(cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  schocos <- foreach::foreach(k = k_min:k_max, .combine = 'rbind') %dopar% {
    schoco(mat, k, homogeneity_thr, completeness_thr)
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
