#' Homogeneity
#'
#' Homogeneity measures how much the samples in a group are similar. The
#' function calculates the percentage of pairs that are similar within the
#' group, relative to the total number of pairs in the group. In the context of
#' adapter sequences of RBPs, homogeneity measures the similarity of adapter
#' sequences within a cluster of RBPs. A high homogeneity score for the cluster
#' indicates that the RBPs in the cluster tend to have similar adapter
#' sequences, while a low score suggests that the adapter sequences are more
#' diverse within the cluster.
#' @param mat adapter_matrix; a matrix of class "adapter_matrix". Use 
#' adapter_matrix() to create a compatible matrix.
#' @param threshold numeric; minimum percentage identity of the identified
#' adapter sequence to be considered a shared adapter.
#' @return A numeric value between 0 and 1. A value of 1 indicates that all
#' samples in the group are similar, while a value of 0 indicates that no
#' samples in the group are similar.
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
#' # calculate homogeneity for ACL 1
#' ids <- clusters |> dplyr::filter(cluster == "ACL 1") |> dplyr::pull(id)
#' index <- which(rownames(amat) %in% ids)
#' homogeneity(mat = amat[index, index])
#' @references https://medium.com/data-science/v-measure-an-homogeneous-and-complete-clustering-ab5b1823d0ad
#' @export
homogeneity <- function(mat, threshold = 0.75) {
  # if (!inherits(mat, "adapter_matrix")) {
  #   msg <- paste0(
  #     "mat must be a matrix of class 'adapter_matrix'. ",
  #     "Use adapter_matrix() to create a compatible matrix."
  #   )
  #   stop(msg)
  # }
  mat <- as.matrix(mat)
  if (nrow(mat) == 1) {
    ho <- NA
  }
  if (nrow(mat) >1) {
    hit <- sum(mat >= threshold)-nrow(mat)
    all <- nrow(mat)^2-nrow(mat)
    ho <- round(hit/all, 3)
  }
  return(ho)
}
