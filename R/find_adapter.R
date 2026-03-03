#' Find adapter region
#'
#' This functions takes a data frame of breakpoints from a pairwise global
#' alignment and detects conserved N-terminal regions, optinally merging
#' neighboring conserved regions.
#' @param bps data.frame; a data frame of breakpoints
#' @param pident_threshold numeric; the lowest accepted percentage identity for
#' a region to be considered a "conserved" region.
#' @param start_threshold numeric; the highest accepted starting position for
#' conserved region to be considered an "N-terminal" region, expressed as the
#' ratio of the aligned sequence length.
#' @param merge logical; whether neighboring conserved regions should be merged.
#' @return A data frame of class "adapter".
#' @details The function returns a single row for each pairwise comparison. If
#' multiple conserved N-terminal regions are identified, the function returns
#' the first one. If no conserved N-terminal regions are identified, the
#' function returns a row with NA values for the start, end, mean_score, and
#' pident columns.
#' @examples
#' data(rbps)
#' ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
#' bps <- find_breakpoints(ps, Nmax = 5)
#' find_adapter(bps)
#' @export
find_adapter <- function(
    bps,
    pident_threshold = 0.4,
    start_threshold = 0.25,
    merge = TRUE
  ) {
  if (!inherits(bps, "breakpoints")) {
    msg <- paste0(
      "bps must be a data frame of class 'breakpoints'. ",
      "Use find_breakpoints() to create a compatible data frame."
    )
    stop(msg)
  }
  class(bps) <- class(bps)[-which(class(bps) == "breakpoints")]
  if (start_threshold < 0 | start_threshold > 1) {
    stop("start_threshold must be between 0 and 1.")
  }
  cnts <- bps[bps$pident >= pident_threshold,]
  if (nrow(cnts) == 0) {
    cnts <- data.frame()
  } else {
    # optionally merge consecutive conserved regions
    if (nrow(cnts) > 1 & merge == TRUE) {
      for (i in 1:(nrow(cnts)-1)) {
        if (cnts$end[i]+1 == cnts$start[i+1]) {
          sum_score1 <- (cnts$end[i]-cnts$start[i]+1)*cnts$mean_score[i]
          sum_pident1 <- (cnts$end[i]-cnts$start[i]+1)*cnts$pident[i]
          sum_score2 <- (cnts$end[i+1]-cnts$start[i+1]+1)*cnts$mean_score[i+1]
          sum_pident2 <- (cnts$end[i+1]-cnts$start[i+1]+1)*cnts$pident[i+1]
          cnts$start[i+1] <- cnts$start[i]
          cnts[i,] <- NA
          newlength <- cnts$end[i+1]-cnts$start[i+1]+1
          cnts$mean_score[i+1] <- round((sum_score1+sum_score2)/newlength, 3)
          cnts$pident[i+1] <- round((sum_pident1+sum_pident2)/newlength, 3)
        }
      }
      cnts <- cnts[stats::complete.cases(cnts),]
    }
    cnt_start_max <- start_threshold*max(bps$end)
    cnts <- cnts[cnts$start <= cnt_start_max,]
  }
  if (nrow(cnts) == 0) {
    cnts <- data.frame(
      pattern_id = unique(bps$pattern_id),
      subject_id = unique(bps$subject_id),
      start = NA,
      end = NA,
      mean_score = NA,
      pident = NA
    )
  }
  cnts <- cnts |> dplyr::slice(1)
  class(cnts) <- c("adapter", class(cnts))
  return(cnts)
}
