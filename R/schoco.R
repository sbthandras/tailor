#' SCores of HOmogeneity and COmpleteness
#'
#' Calculate homogeneity and completeness scores for clusters of RBPs based on
#' the percentage identity of shared adapter sequences. The "schoco" score is a
#' weighted sum of the homogeneity and completeness scores, where homogeneity is
#' weighted by the size of the cluster and completeness is weighted by the log2
#' of the size of the cluster.
#' @param mat adapter_matrix; a matrix of class "adapter_matrix". Use
#' adapter_matrix() to create a compatible matrix.
#' @param k integer; the number of clusters to split the RBPs into.
#' @param homogeneity_thr numeric; minimum percentage identity of the identified
#' adapter sequence to be considered a shared adapter for homogeneity.
#' @param completeness_thr numeric; minimum percentage identity of the
#' identified adapter sequence to be considered a shared adapter for
#' completeness.
#' @return A data frame with 6 columns: nclust (number of clusters requested),
#' cluster (cluser identifier), homogeneity (homogeneity score for the cluster),
#' completeness (completeness score for the cluster), schoco (sum of weighted
#' homogeneity and completeness scores, where homogeneity is multiplied by the
#' size of the cluster and completeness is multiplied by log2 of the size of the
#' cluster. The data frame contains one row for each cluster.
#' @examples
#' data(rbps)
#' data(adapters)
#' # convert to pident matrix and order by the original data frame
#' amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
#' # split matrix into two clusters and calculate schoco values
#' schoco(amat, k = 2)
#' @export
schoco <- function(
    mat,
    k,
    homogeneity_thr = 0.75,
    completeness_thr = 0.75
  ) {
  if (!inherits(mat, "adapter_matrix")) {
    msg <- paste0(
      "mat must be a matrix of class 'adapter_matrix'. ",
      "Use adapter_matrix() to create a compatible matrix."
    )
    stop(msg)
  }
  h <- mat |> stats::dist() |> stats::hclust(method = "average")
  cluster <- paste0("ACL ", stats::cutree(tree = h, k = k))
  unique_clusters <- unique(cluster)
  df <- data.frame()
  for (j in unique_clusters) {
    index <- which(cluster == j)
    smat <- mat[index, index]
    class(smat) <- c("adapter_matrix", class(smat))
    ho <- homogeneity(smat, threshold = homogeneity_thr)
    ho_weighted <- round(length(index)*ho, 3)
    co <- completeness(mat, index, threshold = completeness_thr)
    co_weighted = round(log(length(index),2)*co, 3)
    df <- dplyr::bind_rows(
      df,
      data.frame(
        nclust = k,
        cluster = j,
        count = length(index),
        homogeneity = ho,
        completeness = co,
        schoco = sum(ho_weighted, co_weighted, na.rm = TRUE)
      )
    )
  }
  return(df)
}
