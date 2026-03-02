#' Plot pairwise alignment scores with adapter regions
#'
#' This is a convenience function for looking at pairwise alignments. The
#' function will plot alignment scores along the length of the alignment. At the
#' same time, the function will also attempt to find an adapter region using 
#' various approaches and add the end of the conserved region to the plot.
#' @importFrom rlang .env
#' @param x an object of class `ps`, returned by `position_scores()`.
#' @param type character; strategy for plotting data points. Either "indiv",
#' "ma", or "cusum". See Details for more information.
#' @param method character; the method to use for finding the conserved N
#' terminal region. Can be one or more of "cemean", "cusum", "ewma", "xbarone",
#' and "window".
#' @param highlight character; the method to use for highlighting the conserved
#' N terminal region. Must be one of the methods used in `method`.
#' @param size integer; the number of successive data point to use for
#' calculating moving averages. Only used if `type == "ma"`.
#' @param ... additional arguments (not used).
#' @return a plot.
#' @details When `type = "indiv"` each alignment score will be plotted
#' separately. When `type = "ma"` the plot will show moving averages. When
#' `type = "cusum"` the plot will show the cumulative sum of scores by position.
#' @examples
#' data(rbps)
#' ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
#' plot(ps)
#' @import ggplot2
#' @export
#' @method plot ps
plot.ps <- function(
  x,
  type = "ma",
  method = NULL,
  highlight = NULL,
  size = 10,
  ...
  ) {
  ps <- x
  type <- match.arg(type, choices = c("indiv", "ma", "cusum"))
  if (is.null(method)) {
    method <- c("cemean", "cusum", "ewma", "xbarone", "window")
  } else {
    method <- match.arg(
      method,
      choices = c("cemean", "cusum", "ewma", "xbarone", "window"),
      several.ok = TRUE
    )
  }
  if (!is.null(highlight)) {
    highlight <- match.arg(highlight, choices = method)
  }
  min1row <- function(df) {
    if (nrow(df) == 0) {
      out <- data.frame(
        pattern_id = NA,
        subject_id = NA,
        start = NA,
        end = NA,
        mean_score = NA,
        pident = NA
      )
    } else {
      out <- df
    }
    return(out)
  }
  if ("cemean" %in% method) {
    cemean <- find_breakpoints(ps, method = "cemean") |>
      find_adapter() |>
      min1row()
  }
  if ("window" %in% method) {
    window <- find_breakpoints(ps, method = "window") |>
      find_adapter() |>
      min1row()
  }
  if ("cusum" %in% method) {
    cusum <- find_breakpoints(
      ps,
      method = "plateau",
      type = "cusum"
    ) |> find_adapter() |> min1row()
  }
  if ("ewma" %in% method) {
    ewma <- find_breakpoints(
      ps,
      method = "plateau",
      type = "ewma"
    ) |> find_adapter() |> min1row()
  }
  if ("xbarone" %in% method) {
    xbarone <- find_breakpoints(
      ps,
      method = "plateau",
      type = "xbar.one"
    ) |> find_adapter() |> min1row()
  }
  end_wide <- data.frame(
    cemean = ifelse(
      "cemean" %in% method,
      cemean$end[1], NA),
    cusum = ifelse(
      "cusum" %in% method,
      cusum$end[1], NA),
    ewma = ifelse(
      "ewma" %in% method,
      ewma$end[1], NA),
    xbarone = ifelse(
      "xbarone" %in% method,
      xbarone$end[1], NA),
    window = ifelse(
      "window" %in% method,
      window$end[1], NA)
  )
  end_long <- end_wide %>%
    tidyr::pivot_longer(
      cols = names(end_wide),
      names_to = "method",
      values_to = "position"
    ) %>%
    dplyr::arrange("position") %>%
    dplyr::filter(.data$method %in% .env$method)
  if (is.null(highlight)) {
    highlight_value <- end_long |>
      dplyr::pull(.data$position) |>
      stats::na.omit() |>
      max()
  } else {
    highlight_value <- end_long |>
      dplyr::filter(.data$method == highlight) |>
      dplyr::pull(.data$position)
  }
  # moving averages
  avg <- vector()
  for (i in size:length(ps$position_scores$score-size+1)) {
    avg <- c(avg, mean(ps$position_scores$score[(i-size+1):i]))
  }
  # first s values are identical
  avg <- c(rep(avg[1], times = size-1), avg)
  res <- data.frame(
    position = seq_along(ps$position_scores$score),
    score = ps$position_scores$score,
    ma = avg,
    cumsum = cumsum(ps$position_scores$score)
  )
  if(is.na(highlight_value)) {
    res$domain <- "C-terminal"
  } else {
    res$domain <- ifelse(
      res$position <= highlight_value,
      "N-terminal",
      "C-terminal"
    )
  }
  if (type == "indiv") {
    g <- ggplot(res, aes(.data$position, .data$score))
  }
  if (type == "ma") {
    g <- ggplot(res, aes(.data$position, .data$ma))
  }
  if (type == "cusum") {
    g <- ggplot(res, aes(.data$position, .data$cumsum))
  }
  g <- g +
    geom_point(aes(color = .data$domain), alpha = 0.5) +
    scale_color_manual(
      values = c("N-terminal" = "red", "C-terminal" = "blue"),
      na.translate = FALSE
    ) +
    geom_vline(
      data = end_long,
      aes(
        xintercept = .data$position,
        group = .data$method,
        linetype = .data$method
      ),
      key_glyph = "path"
    )
  if (length(method) == 1) {
    g <- g +
      guides(
        color = guide_legend(title = "Domain", order = 1),
        linetype = "none"
      )
  } else {
    g <- g +
      guides(
        color = guide_legend(title = "Domain", order = 1),
        linetype = guide_legend(title = "Method", order = 2)
      )
  }
  g <- g +
    xlab("Amino acid position from N terminal") +
    theme_bw()

  # Add dynamic y-axis label
  if (type == "indiv") {
    g <- g +
      ylab(paste0(
        "Substitution score\n(",
        ps$substitution_matrix,
        ")"
      ))
  }
  if (type == "ma") {
    g <- g +
      ylab(paste0(
        size,
        " point moving average of substitution scores\n(",
        ps$substitution_matrix,
        ")"
      ))
  }
  if (type == "cusum") {
    g <- g +
      ylab(paste0(
        "Cumulative sum of substitution scores\n(",
        ps$substitution_matrix,
        ")"
      ))
  }
  return(g)
}
