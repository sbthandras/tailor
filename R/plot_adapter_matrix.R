#' Plot adapter matrix
#'
#' Take an adapter matrix and optionally data frame of cluster designations and
#' plot them on a simple heatmap.
#' @import ggplot2 patchwork
#' @param x adapter_matrix; a matrix of class "adapter_matrix". Use
#' adapter_matrix() to create a compatible matrix.
#' @param clusters data.frame; a table of RBP IDs and cluster designations.
#' @param ... additional arguments (not used).
#' @return A heatmap. If clusters are provided, also show cluster assignments.
#' @examples
#' data(rbps)
#' data(adapters)
#' # convert to pident matrix and order by the original data frame
#' amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
#' # plot pident matrix
#' plot(amat)
#' # split RBPs into two clusters
#' clusters <- cluster_adapters(amat, k_min = 2, k_max = 2)
#' # plot pident matrix with cluster assignments
#' plot(amat, clusters = clusters)
#' @export
#' @method plot adapter_matrix
plot.adapter_matrix <- function(x, clusters = NULL, ...) {
  mat <- x
  mat_long <- reshape2::melt(mat, value.name = "pident")
  p_main <- ggplot2::ggplot(
      mat_long,
      ggplot2::aes(
        x = !!rlang::sym("Var1"), 
        y = !!rlang::sym("Var2"),
        fill = !!rlang::sym("pident"))
    ) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = round(!!rlang::sym("pident"), 1))) +
    ggplot2::scale_fill_gradient(low = "white", high = "salmon") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "", y = "", fill = "") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  if (is.null(clusters)) return(p_main)
  row_anno <- clusters |>
    dplyr::rename(Var2 = !!rlang::sym("id")) |>
    dplyr::mutate(
      Var2 = factor(!!rlang::sym("Var2"), levels = unique(mat_long$Var2)))
  p_row <- ggplot2::ggplot(
    row_anno, 
    ggplot2::aes(
      x = 1, y = !!rlang::sym("Var2"), fill = !!rlang::sym("cluster"))) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_discrete(
      limits = levels(mat_long$Var2), expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_discrete() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "none"
    )
  col_anno <- clusters |>
    dplyr::rename(Var1 = !!rlang::sym("id")) |>
    dplyr::mutate(Var1 = factor(!!rlang::sym("Var1"), levels = unique(mat_long$Var1)))
  p_col <- ggplot2::ggplot(
    col_anno, ggplot2::aes(x = !!rlang::sym("Var1"), y = 1, fill = !!rlang::sym("cluster"))) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(limits = levels(mat_long$Var1), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_discrete() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "none"
    )
  fig <- patchwork::wrap_plots(
    p_col, p_main, p_row,
    design = c(
      patchwork::area(1, 1, 1, 19),
      patchwork::area(2, 1, 20, 19),
      patchwork::area(2, 20, 20, 20)
    )
  )
  return(fig)
}

