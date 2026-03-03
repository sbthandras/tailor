#' Find breakpoints in global alignment scores
#'
#' This function looks at global alignment scores by position and attempts to
#' find locations along the length of the aligned sequences where these scores
#' tend to change substantially.
#' @name find_breakpoints
#' @param position_scores an object of class `ps`, returned by
#' `position_scores()`.
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
find_breakpoints <- function(position_scores, method = "cemean", ...) {
  if (!inherits(position_scores, "ps")) {
    stop("position_scores must be an object created by position_scores().")
  }
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
  Nmax <- as.integer(Nmax)
  if (!is.integer(Nmax)) stop("Nmax must be an integer.")
  if (Nmax < 1) stop("Nmax must be a positive integer.")
  seed <- as.integer(seed)
  if (!is.integer(seed)) stop("seed must be an integer")
  set.seed(seed)
  scores <-  position_scores$position_scores
  bps <- breakpoint::CE.Normal.Mean(as.data.frame(scores$score), Nmax = Nmax)
  if (inherits(bps, "character")){
    if (bps == "No Break-Points are Estimated") {
      bpdf <- data.frame(
        pattern_id = position_scores$pattern_id,
        subject_id = position_scores$subject_id,
        start = 1,
        end = nrow(scores),
        mean_score = round(mean(scores$score), 3),
        pident = round(mean(scores$identity), 3))
      return(bpdf)
    }
    if (bps == "Error in data : single column dataframe only") {
      print(paste0(position_scores$pattern_id, " ", position_scores$subject_id))
      stop()
    }
  }
  bps <- bps$BP.Loc
  region_start <- c(1,bps+1)
  region_end <- c(bps, nrow(scores))
  region_mean_score <- mapply(
    function(x, y) round(mean(scores$score[x:y]), 3),
    region_start,
    region_end)
  region_pident <- mapply(
    function(x,y) round(mean(scores$identity[x:y]), 3),
    region_start,
    region_end)
  bpdf <- data.frame(
    pattern_id = position_scores$pattern_id,
    subject_id = position_scores$subject_id,start = region_start,
    end = region_end,
    mean_score = region_mean_score,
    pident = region_pident)
  return(bpdf)
}

#' @rdname find_breakpoints
#' @param pinit numeric; ratio of sequence length used for setting the center.
#' See Details for more information.
#' @param type character; the type of statistical process control chart to use.
#' Currently supported cards are `cusum`, `ewma` and `xbar.one`.
#' @details When using `method = "plateau"` the function will look at the first
#' few positions of the alignment (goverened by `pinit`) and calculate a mean
#' score. This mean score will be used as the center: 1. it will be used as the
#' center point for a number of control charts; 2. any scores above the center
#' will be reduced to be indentical with the center. Then the function will use
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

  init_length <- round(pinit*length(position_scores$position_scores$score), 0)
  center <- mean(position_scores$position_scores$score[1:init_length])

  type <- match.arg(type, choices = c("cusum", "ewma", "xbar.one"))

  # adjust position scores to plateau
  position_scores$position_scores$score <- ifelse(
    position_scores$position_scores$score > center,
    center,
    position_scores$position_scores$score
  )

  if (type == "xbar.one") {
    spc <- qcc::qcc(
      position_scores$position_scores$score,
      type = "xbar.one",
      center = center,
      plot = FALSE,
      ...
    )
    violations <- c(spc$violations$bexond.limits, spc$violations$violating.runs)
    violations <- sort(unique(violations))
  } else if (type == "ewma") {
    spc <- qcc::ewma(
      position_scores$position_scores$score,
      center = center,
      plot = FALSE,
      ...
    )
    violations <- sort(unname(spc$violations))
  } else if (type == "cusum") {
    spc <- qcc::cusum(
      position_scores$position_scores$score,
      center = center,
      plot = FALSE,
      ...
    )
    violations <- sort(unique(spc$violations$lower))
  }
  position_scores$position_scores$pass <- sapply(1:nrow(position_scores$position_scores), function(x) {
    if (x %in% violations) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }) |> unname()

  runlen <- rle(position_scores$position_scores$pass)
  if (length(runlen$lengths) == 1) {
    out <- data.frame(
      pattern_id = position_scores$pattern_id,
      subject_id = position_scores$subject_id,
      start = 1,
      end = nrow(position_scores$position_scores)
    )
  } else {
    out <- data.frame(
      pattern_id = position_scores$pattern_id,
      subject_id = position_scores$subject_id,
      start = c(
        1, cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]) + 1
      ),
      end = c(
        cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]),
        nrow(position_scores$position_scores)
      )
    )
  }
  out$mean_score <- mapply(
    function(x, y) round(mean(position_scores$position_scores$score[x:y]), 3),
    out$start,
    out$end)
  out$pident <- mapply(
    function(x,y) round(mean(position_scores$position_scores$identity[x:y]), 3),
    out$start,
    out$end)
  return(out)
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
  ngroups <- floor(nrow(position_scores$position_scores)/window)
  position_scores$position_scores$group <- c(
    rep(1:ngroups, each = window),
    rep(ngroups + 1, times = nrow(position_scores$position_scores) - window*ngroups)
  )
  ms <- position_scores$position_scores %>%
    dplyr::group_by(!!rlang::sym("group")) %>%
    dplyr::summarise(
      window = dplyr::n(),
      mean = round(mean(!!rlang::sym("score")), 3),
      identity_p = suppressWarnings(
        stats::wilcox.test(
          as.numeric(!!rlang::sym("identity")), 
          alternative = "less", mu = pident_threshold)
      )$p.value,
      score_p = suppressWarnings(
        stats::wilcox.test(
          !!rlang::sym("score"), 
          alternative = "less", mu = score_threshold)
      )$p.value
    )
  ms$pass <- ms$score_p >= p_threshold | ms$identity_p >= p_threshold
  runlen <- rle(ms$pass)
  if (length(runlen$lengths) == 1) {
    out <- data.frame(
      pattern_id = position_scores$pattern_id,
      subject_id = position_scores$subject_id,
      start = 1,
      end = nrow(position_scores$position_scores)
    )
  } else {
    out <- data.frame(
      pattern_id = position_scores$pattern_id,
      subject_id = position_scores$subject_id,
      start = c(
        1,
        window * cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]) + 1
      ),
      end = c(
        window * cumsum(runlen$lengths[1:(length(runlen$lengths)-1)]),
        nrow(position_scores$position_scores)
      )
    )
  }
  out$mean_score <- mapply(
    function(x, y) round(mean(position_scores$position_scores$score[x:y]), 3),
    out$start,
    out$end)
  out$pident <- mapply(
    function(x,y) round(mean(position_scores$position_scores$identity[x:y]), 3),
    out$start,
    out$end)
  return(out)
}

