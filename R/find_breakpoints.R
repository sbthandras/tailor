#' Find breakpoints in global alignment scores
#'
#' This function looks at global alignment scores by position and attempts to
#' find locations along the length of the aligned sequences where these scores
#' tend to change substantially.
#' @name find_breakpoints
#' @param position_scores an object of class `ps`, returned by
#' `position_scores()`.
#' @param max_end integer; the maximum position to look for breakpoints. If 
#' `NULL`, search the full alignment, if an integer, limit the search to that 
#' position. If both this argument and `position_scores$max_end` are set, the 
#' smaller integer is used.
#' @param method character; the method to use for finding breakpoints. Either
#' "cemean", "plateau", or "window".
#' @param ... additional arguments depending on the selected method. When
#' `method = "plateau"`, additional arguments passed on to the selected control
#' chart. See `?qcc::cusum()`, `?qcc::ewma()`, or `?qcc::qcc()` for a full list
#' of control chart parameters.
#' @return A data frame with 6 columns:
#' \itemize{
#'   \item `pattern_id`: pattern ID
#'   \item `subject_id`: subject ID
#'   \item `start`: first position for a section of the alignment
#'   \item `end`: last position for a section of the alignment
#'   \item `mean_score`: mean score for amino acid pairs within the section
#'   \item `pident`: ratio of matching amino acid pairs within the section
#' }
#' @seealso [position_scores()]
#' @export
find_breakpoints <- function(
  position_scores, 
  max_end = NULL,
  method = "cemean", 
  ...
) {
  if (!inherits(position_scores, "ps")) {
    stop("position_scores must be an object created by position_scores().")
  }
  scores_df <- position_scores$position_scores
  n <- nrow(scores_df)
  if (!is.null(max_end)) {
    max_end <- as.integer(max_end)
    if (
      length(max_end) != 1L || 
      !is.integer(max_end) || 
      max_end < 1L || 
      max_end > n
    ) {
      stop(
        "max_end must be a single integer ",
        "between 1 and the length of the alignment."
      )
    }
  }
  limit <- min(c(n, max_end, position_scores$max_end), na.rm = TRUE)
  # remove alignment positions beyond the limit
  position_scores$position_scores <- scores_df[seq_len(limit), ]
  method <- match.arg(method, choices = c("cemean", "plateau", "window"))
  f <- get(paste0("find_breakpoints_", method))
  out <- f(position_scores, ...)
  class(out) <- c("breakpoints", class(out))
  return(out)
}

#' @rdname find_breakpoints
#' @param Nmax integer; the maximum number of breakpoints to look for.
#' @param seed integer; random seed.
#' @details When using `method = "cemean"` the function will estimate the number
#' of breakpoints and their locations using the Cross-Entropy Method. This
#' method will use the `CE.Normal.Mean()` function from the `breakpoint` package
#' to find the breakpoints.
#' @note When using `method == "cemean"` the speed of the function will depend
#' on the number of breakpoints to look for.
#' @examples
#' data(rbps)
#' ps <- position_scores("MN395291-1", "ON513429-1", data = rbps)
#'
#' # Find breakpoints using the Cross-Entropy Method
#' find_breakpoints(ps, method = "cemean", Nmax = 5)
find_breakpoints_cemean <- function(
  position_scores,
  Nmax = 5,
  seed = 0
) {
  if ("ps" %in% class(position_scores) == FALSE) {
    stop("position_scores must be an object created by position_scores().")
  }
  scores <- position_scores$position_scores$score
  identity <- position_scores$position_scores$identity
  core <- breakpoints_cemean_core(scores, Nmax = Nmax, seed = seed)
  pident <- mapply(
    function(s, e) round(mean(identity[s:e]), 3),
    core$start, core$end, USE.NAMES = FALSE
  )
  bpdf <- data.frame(
    pattern_id = position_scores$pattern_id,
    subject_id = position_scores$subject_id,
    start = core$start,
    end = core$end,
    mean_score = core$mean_score,
    pident = pident
  )
  return(bpdf)
}


#' Core algorithm for Cross-Entropy Method breakpoints
#' @keywords internal
#' @param scores numeric vector of alignment scores
#' @param Nmax integer; maximum number of breakpoints to look for
#' @param seed integer; random seed
#' @return data.frame with region_start, region_end, region_mean_score
#' @noRd
breakpoints_cemean_core <- function(scores, Nmax = 5L, seed = 0L) {
  Nmax <- as.integer(Nmax)
  if (!is.integer(Nmax)) stop("Nmax must be an integer.")
  if (Nmax < 1) stop("Nmax must be a positive integer.")
  seed <- as.integer(seed)
  if (!is.integer(seed)) stop("seed must be an integer")
  set.seed(seed)
  bps <- breakpoint::CE.Normal.Mean(as.data.frame(scores), Nmax = Nmax)
  if (inherits(bps, "character")){
    if (bps == "No Break-Points are Estimated") {
      return(data.frame(
        start = 1L,
        end = length(scores),
        mean_score = round(mean(scores), 3)
      ))
    } else {
      stop(bps)
    }
  }
  bps <- bps$BP.Loc
  region_start <- c(1L, bps + 1L)
  region_end <- c(bps, length(scores))
  region_mean_score <- mapply(
    function(x, y) round(mean(scores[x:y]), 3),
    region_start, region_end, USE.NAMES = FALSE
  )
  return(data.frame(
    start = region_start,
    end = region_end,
    mean_score = region_mean_score
  ))
}

#' @rdname find_breakpoints
#' @param pinit numeric; ratio of sequence length used for setting the center.
#' See Details for more information.
#' @param type character; the type of statistical process control chart to use.
#' Currently supported cards are `cusum`, `ewma` and `xbar.one`.
#' @details When using `method = "plateau"` the function will look at the first
#' few positions of the alignment (governed by `pinit`) and calculate a mean
#' score. This mean score will be used as the center: 1. it will be used as the
#' center point for a number of control charts; 2. any scores above the center
#' will be reduced to be identical with the center. Then the function will use
#' the control chart selected by `type`, and extract sections using the
#' violations.
#' @examples
#' # Find breakpoints using plateau on scores and EWMA chart
#' find_breakpoints(ps, method = "plateau", type = "ewma", lambda = 0.2)
#' # Find breakpoints using plateau on scores and CUSUM chart
#' find_breakpoints(ps, method = "plateau", type = "cusum")
find_breakpoints_plateau <- function(
    position_scores,
    pinit = 0.05,
    type = "ewma",
    ...
    ) {
  scores <- position_scores$position_scores$score
  identity <- position_scores$position_scores$identity
  core <- breakpoints_plateau_core(scores, pinit = pinit, type = type, ...)
  pident <- mapply(
    function(s, e) round(mean(identity[s:e]), 3), 
    core$start, core$end, USE.NAMES = FALSE
  )
  bpdf <- data.frame(
    pattern_id = position_scores$pattern_id,
    subject_id = position_scores$subject_id,
    start = core$start,
    end = core$end,
    mean_score = core$mean_score,
    pident = pident
  )
  return(bpdf)
}

#' Core algorithm for plateau method breakpoints
#' @keywords internal
#' @param scores numeric vector of alignment scores
#' @param pinit numeric; fraction of the sequence used to set the centre
#' @param type character; one of cusum, ewma, xbar.one
#' @param ... additional arguments forwarded to the qcc functions
#' @return data.frame with columns start, end, mean_score, pident
#' @noRd
breakpoints_plateau_core <- function(scores, pinit = 0.05, type = "ewma", ...) {
  init_length <- round(pinit * length(scores), 0)
  center <- mean(scores[1:init_length])
  scores_plateau <- ifelse(scores > center, center, scores)
  type <- match.arg(type, choices = c("cusum", "ewma", "xbar.one"))
  if (type == "xbar.one") {
    spc <- qcc::qcc(
      scores_plateau,
      type = "xbar.one",
      center = center,
      plot = FALSE,
      ...
    )
    violations <- c(spc$violations$bexond.limits, spc$violations$violating.runs)
    violations <- sort(unique(violations))
  } else if (type == "ewma") {
    spc <- qcc::ewma(
      scores_plateau,
      center = center,
      plot = FALSE,
      ...
    )
    violations <- sort(unname(spc$violations))
  } else if (type == "cusum") {
    spc <- qcc::cusum(
      scores_plateau,
      center = center,
      plot = FALSE,
      ...
    )
    violations <- sort(unique(spc$violations$lower))
  }
  pass <- !(seq_along(scores) %in% violations)
  runlen <- rle(pass)
  if (length(runlen$lengths) == 1) {
    start <- 1L
    end <- length(scores)
  } else {
    start <- c(1, cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]) + 1)
    end <- c(
      cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]),
      length(scores)
    )
  }
  data.frame(
    start = start,
    end = end,
    mean_score = mapply(
      function(s, e) round(mean(scores[s:e]), 3), 
      start, end, USE.NAMES = FALSE
    )
  )
}

#' @rdname find_breakpoints
#' @param window integer; the number of successive amino acids to look at
#' @param score_threshold numeric; mu threshold for alignment scores
#' @param pident_threshold numeric; mu threshold for amino acid identities
#' @param p_threshold numeric; p value threshold of significance
#' @details When using `method == "window"` the function will form disjunct
#' windows along the length of the alignment and look at the positions within
#' each window. For each window, it will perform one sided wilcoxon tests
#' against the `score_threshold` and the `pident_threshold` and consider the
#' window conserved if at least one of the tests passes. Finally it will
#' collapse neighboring windows that are all considered conserved or
#' non-conserved and return a data frame of sections in which conservation
#' alternates between the lines.
#' @importFrom dplyr %>%
#' @examples
#' # Find breakpoints using disjunct windows
#' find_breakpoints(ps, method = "window", window = 5)
find_breakpoints_window <- function(
    position_scores,
    window =  5,
    score_threshold = 3,
    pident_threshold = 0.4,
    p_threshold = 0.05
) {
  scores <- position_scores$position_scores$score
  identity <- position_scores$position_scores$identity
  core <- breakpoints_window_core(
    scores, 
    identity,
    window = window,
    score_threshold = score_threshold,
    pident_threshold = pident_threshold,
    p_threshold = p_threshold
  )
  identity <- mapply(
    function(s, e) round(mean(identity[s:e]), 3), 
    core$start, core$end, USE.NAMES = FALSE
  )
  out <- data.frame(
    pattern_id = position_scores$pattern_id,
    subject_id = position_scores$subject_id,
    start = core$start,
    end = core$end,
    mean_score = core$mean_score,
    pident = identity
  )
  return(out)
}

#' Core algorithm for window method breakpoints
#' @keywords internal
#' @param scores numeric vector of alignment scores
#' @param identity numeric vector of identity values (same length as scores)
#' @param window integer; size of each sliding window
#' @param score_threshold numeric; mu threshold for alignment scores
#' @param pident_threshold numeric; mu threshold for amino acid identities
#' @param p_threshold numeric; p-value threshold for significance
#' @return data.frame with columns start, end, mean_score, pident
#' @noRd
breakpoints_window_core <- function(
    scores, 
    identity = NULL,
    window = 5L,
    score_threshold = 3,
    pident_threshold = 0.4,
    p_threshold = 0.05
  ) {
  n <- length(scores)
  ngroups <- floor(n / window)
  group <- c(
    rep(1:ngroups, each = window),
    rep(ngroups + 1, times = n - window * ngroups)
  )
  ms <- data.frame(score = scores, group = group)
  if (!is.null(identity)) {
    ms$identity <- identity
  }
  ms <- ms %>%
    dplyr::group_by(!!rlang::sym("group")) %>%
    dplyr::summarise(
      window = dplyr::n(),
      mean = round(mean(!!rlang::sym("score")), 3),
      score_p = suppressWarnings(
        stats::wilcox.test(
          !!rlang::sym("score"), 
          alternative = "less", mu = score_threshold)
      )$p.value,
      identity_p = if (!is.null(identity)) {
        suppressWarnings(
          stats::wilcox.test(
            as.numeric(!!rlang::sym("identity")), 
            alternative = "less", mu = pident_threshold)
        )$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  if (is.null(identity)) {
    ms$pass <- ms$score_p >= p_threshold
  } else {
    ms$pass <- ms$score_p >= p_threshold | ms$identity_p >= p_threshold
  }
  runlen <- rle(ms$pass)
  if (length(runlen$lengths) == 1) {
    start <- 1L
    end <- n
  } else {
    start <- c(
      1, 
      window * cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]) + 1
    )
    end <- c(window * cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]), n)
  }
  data.frame(
    start = start,
    end = end,
    mean_score = mapply(
      function(s, e) round(mean(scores[s:e]), 3), 
      start, end, USE.NAMES = FALSE
    )
  )
}

