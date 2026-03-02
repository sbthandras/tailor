#' Completeness
#'
#' Completeness measures measures how much similar samples are put together by
#' the clustering algorithm. In the context of adapter sequences of RBPs,
#' completeness measures whether RBPs that share adapters with members of a
#' cluster are also included within that cluster. High completeness score for a
#' cluster indicates that when cluster members share adapters with other RBPs,
#' those other RBPs tend to be in the same cluster as well. On the other hand,
#' low scores suggests that cluster members often share adapters with RBPs that
#' do not belong to the cluster.
#' @param mat adapter_matrix; a matrix of class "adapter_matrix". Use 
#' adapter_matrix() to create a compatible matrix.
#' @param index numeric; indices of proteins that are included within a cluster.
#' @param threshold numeric; minimum percentage identity of the identified
#' adapter sequence to be considered a shared adapter.
#' @return A numeric value between 0 and 1. A value of 1 indicates that all
#' proteins sharing features with cluster members are contained within the
#' cluster, while a value of 0 indicates that shared features are found equally
#' inside and outside the cluster.
#' @note Note that \code{mat} is different for homogeneity and completeness.
#' When calculating homogeneity, we are only interested in comparisons within
#' a clusters and so \code{mat} represents a cluster. However, when calculating
#' completeness, we are also interested in comparisons outside a cluster, and
#' so \code{mat} must contain all comparisons, not just the cluster in question.
#' @examples
#' data(rbps)
#' data(adapters)
#' # convert to pident matrix and order by the original data frame
#' amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
#' # split to two clusters
#' clusters <- cluster_adapters(amat, k_min = 2, k_max = 2)
#' # calculate completeness for ACL 1
#' ids <- clusters |> dplyr::filter(cluster == "ACL 1") |> dplyr::pull(id)
#' index <- which(rownames(amat) %in% ids)
#' completeness(mat = amat, index = index)
#' @references https://medium.com/data-science/v-measure-an-homogeneous-and-complete-clustering-ab5b1823d0ad
#' @export
completeness <- function(mat, index, threshold = 0.75) {
  # if (!inherits(mat, "adapter_matrix")) {
  #   msg <- paste0(
  #     "mat must be a matrix of class 'adapter_matrix'. ",
  #     "Use adapter_matrix() to create a compatible matrix."
  #   )
  #   stop(msg)
  # }
  if (length(index) == 0){
    stop("cluster must contain at least one element. Please provide indices.")
  }
  imat <- mat[index, index]
  emat <- mat[index, -index]
  amat <- mat[index, ]
  hit <- sum(emat >= threshold)
  all <- sum(amat >= threshold)-length(index)-(sum(imat >= threshold)-length(index))/2
  co <- round(1-hit/all, 3)
  return(co)
}
