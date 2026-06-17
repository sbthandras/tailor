#' Find adapter region
#'
#' This functions takes a data frame of breakpoints from a pairwise global
#' alignment and detects conserved N-terminal regions, optinally merging
#' neighboring conserved regions.
#' @param bps data.frame; a data frame of breakpoints
#' @param min_pident numeric; the lowest accepted percentage identity for a
#' region to be considered a "conserved" region.
#' @param max_start integer; the maximum accepted starting position for the
#' first conserved region.
#' @param min_width integer; the minimum accepted width for a conserved region
#' to be considered an adapter.
#' @param merge_beginning logical; whether to merge the positions before the
#' first conserved region. See Details for more information.
#' @param merge_conserved logical; whether neighboring conserved regions should
#' be merged. See Details for more information.
#' @return A data frame of class "adapter".
#' @details The function returns a single row for each pairwise comparison. If
#' no conserved N-terminal regions are identified, the function returns a row
#' with NA values for the start, end, mean_score, and pident columns. If a
#' single conserved N_terminal region is identified, the function returns the
#' start and end positions, mean score, and percentage identity of this region.
#' If multiple conserved N-terminal regions are identified, behaviour depends on
#' the values of `merge_conserved` and `merge_beginning`.
#' @details If `merge_conserved` is TRUE, neighboring conserved regions will be
#' merged and the mean score and percentage identity of the merged region will
#' be recalculated. If  `merge_conserved` is FALSE, the function will only keep
#' the first conserved region.
#' @details If `merge_beginning` is TRUE, and there is a conserved region
#' (merged or unmerged) which starts within `max_start`, but not at position 1,
#' then the function will also merge the positions before this conserved region
#' and recalculate the mean score and percentage identity of the merged region.
#' If the percentage identity of the merged region is below `min_pident`, the
#' function will drop the merged region and return NA values for the start, end,
#' mean_score, and pident columns.
#' @examples
#' data(rbps)
#' ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
#' bps <- find_breakpoints(ps, Nmax = 5)
#' find_adapter(bps)
#' @export
find_adapter <- function(
    bps,
    min_pident = 0.4,
    max_start = 10,
    min_width = 7,
    merge_beginning = TRUE,
    merge_conserved = TRUE
  ) {
  if (!inherits(bps, "breakpoints")) {
    msg <- paste0(
      "bps must be a data frame of class 'breakpoints'. ",
      "Use find_breakpoints() to create a compatible data frame."
    )
    stop(msg)
  }
  class(bps) <- class(bps)[-which(class(bps) == "breakpoints")]
  # process conserved regions
  consreg <- bps[bps$pident >= min_pident,]
  if (nrow(consreg > 1)) {
    # merge consecutive conserved regions
    if (merge_conserved) consreg <- merge_regions(consreg)
    # if there are still multiple, keep the first
    consreg <- consreg[1, ]
  }
  if (nrow(consreg) == 0 || consreg$start > max_start) {
    adapter <- data.frame(
      pattern_id = unique(bps$pattern_id),
      subject_id = unique(bps$subject_id),
      start = NA_integer_,
      end = NA_integer_,
      mean_score = NA_real_,
      pident = NA_real_
    )
  } else if (consreg$start != 1 && merge_beginning) {
    adapter <- merge_regions(
      dplyr::bind_rows(
        bps[which(bps$end < consreg$start), ],
        consreg
      )
    )
  } else {
    adapter <- consreg
  }
  if (!is.na(adapter$pident) && adapter$pident < min_pident) {
    adapter <- data.frame(
      pattern_id = unique(bps$pattern_id),
      subject_id = unique(bps$subject_id),
      start = NA_integer_,
      end = NA_integer_,
      mean_score = NA_real_,
      pident = NA_real_
    )
  }
  if (!is.na(adapter$start) && (adapter$end - adapter$start + 1) < min_width) {
    adapter <- data.frame(
      pattern_id = unique(bps$pattern_id),
      subject_id = unique(bps$subject_id),
      start = NA_integer_,
      end = NA_integer_,
      mean_score = NA_real_,
      pident = NA_real_
    )
  }
  class(adapter) <- c("adapter", class(adapter))
  return(adapter)
}

merge_regions <- function(df) {
  for (i in seq_len(nrow(df) - 1)) {
    if (
      !is.na(df$end[i]) &&
      !is.na(df$start[i + 1]) &&
      df$end[i] + 1 == df$start[i + 1]
    ) {
      sum_score1 <- (df$end[i] - df$start[i] + 1) * df$mean_score[i]
      sum_pident1 <- (df$end[i] - df$start[i] + 1) * df$pident[i]
      sum_score2 <- (df$end[i + 1] - df$start[i + 1] + 1) * df$mean_score[i + 1]
      sum_pident2 <- (df$end[i + 1] - df$start[i + 1] + 1) * df$pident[i + 1]
      df$start[i + 1] <- df$start[i]
      df[i, ] <- NA
      newlength <- df$end[i + 1] - df$start[i + 1] + 1
      df$mean_score[i + 1] <- round((sum_score1 + sum_score2) / newlength, 3)
      df$pident[i + 1] <- round((sum_pident1 + sum_pident2) / newlength, 3)
    }
  }
  df |> stats::na.omit()
}
