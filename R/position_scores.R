#' Pairwise global alignments with position scores
#'
#' This function takes a pair of IDs within a data set, performs a pairwise 
#' global alignment on their sequences, and returns an object which contains 
#' the pattern and subject IDs and a tibble of nucleotide or amino acid 
#' identities and alignment scores by position.
#' @param pattern_id character; pattern ID
#' @param subject_id character; subject ID
#' @param data data.frame; a data frame that contains IDs and sequences
#' @param id_var character; variable within \code{data} that stores IDs.
#' @param seq_var character; variable within \code{data} that stores sequences.
#' @param submat character; the substitution matrix. Options include
#' \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"}, \code{"BLOSUM80"},
#' \code{"BLOSUM100"}, \code{"PAM30"}, \code{"PAM40"}, \code{"PAM70"},
#' \code{"PAM120"}, \code{"PAM250"}.
#' @param max_end_vars character; optional vector of variable names in 
#' \code{data} that contain numeric values representing positions which the 
#' adapters should not overlap with and terminate earlier (e.g., domain 
#' start positions). See Details for more information.
#' @param verbose logical; should verbose messages be printed to the console?
#' @return an object of class \code{ps} which is a list with five elements:
#' \code{pattern_id}, \code{subject_id}, \code{position_scores} and
#' \code{substitution_matrix}, \code{max_end}. \code{position scores} is a 
#' tibble that contains nucleotide or amino acid identities and alignment scores 
#' by position.
#' @details If you provide one or more `max_end_var` names, the function will
#' look for their values in both the pattern and subject rows, calculate their 
#' minimum  and return it in the output object. The value will then be 
#' respected by `find_breakpoints()` and the function will not look at 
#' positions beyond this value.
#' @examples
#' data(rbps)
#' position_scores("MN395291-1", "ON513429-1", data = rbps)
#' @seealso \code{find_breakpoints()}
#' @export
position_scores <- function(
    pattern_id,
    subject_id,
    data,
    id_var = "Core_ORF",
    seq_var = "translation",
    submat = "BLOSUM80",
    max_end_vars = NULL,
    verbose = getOption("verbose")
  ) {
  if (!inherits(data, "data.frame")) {
    stop("Argument 'data' must be a data.frame.")
  }
  if (!is.null(max_end_vars)) {
    if (!is.character(max_end_vars)) {
      stop("Argument 'max_end_vars' must be a character vector or NULL.")
    }
    if (!all(max_end_vars %in% names(data))) {
      stop("All variables in 'max_end_vars' must exist in 'data'.")
    }
  }
  id_var <- match.arg(id_var, names(data))
  seq_var <- match.arg(seq_var, names(data))
  if (!pattern_id %in% data[[id_var]]) {
    stop(paste0("Pattern ID not found: ", pattern_id))
  }
  if (!subject_id %in% data[[id_var]]) {
    stop(paste0("Subject ID not found: ", subject_id))
  }
  df <- data |> dplyr::filter(!!rlang::sym(id_var) %in% c(pattern_id, subject_id))
  df[[seq_var]] <- gsub("\\*", "", df[[seq_var]])
  valid_input <- validate_rbps(
    df, id_var = id_var, seq_var = seq_var, verbose = verbose
  )
  if (!valid_input) {
    stop(
      "Argument 'data' is not valid. See validate_rbps() for more information.")
  }
  pat_index_x <- which(df[[id_var]] == pattern_id)
  pat_index_y <- which(names(df) == seq_var)
  pat <- as.character(df[pat_index_x, pat_index_y])
  pat <- gsub("\\*", "", pat)
  seq_index_x <- which(df[[id_var]] == subject_id)
  seq_index_y <- which(names(df) == seq_var)
  seq <- as.character(df[seq_index_x, seq_index_y])
  seq <- gsub("\\*", "", seq)
  submat <- match.arg(submat, c(
    "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
    "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"
  ))
  if (!exists(submat)) {
    datapath <- paste0("/data/",submat,".rda")
    load(system.file(datapath,package = "pwalign"))
  }
  sm <- get(submat)
  rownames(sm) <- gsub("\\*", "-", rownames(sm))
  colnames(sm) <- gsub("\\*", "-", colnames(sm))
  ali <- pwalign::pairwiseAlignment(
    pat, seq, type = "global", substitutionMatrix = submat)
  sp <- strsplit(as.character(pwalign::pattern(ali)), split = "")[[1]]
  ss <- strsplit(as.character(pwalign::subject(ali)), split = "")[[1]]
  scs <- mapply(function(x, y){
    i <- grep(x, rownames(sm))
    j <- grep(y, colnames(sm))
    return(sm[i,j])
  }, sp, ss)
  identity <- sp == ss
  position_scores <- tibble::tibble(
    pattern = sp,
    subject = ss,
    identity = identity,
    score = scs
  )
  if (!is.null(max_end_vars)) {
    max_end <- sapply(max_end_vars, function(x) {
      c(df[[x]][pat_index_x], df[[x]][seq_index_x])
    }) |> as.numeric() |> stats::na.omit() |> unique()
    max_end <- max_end[max_end != 0]
    if (length(max_end) == 0) {
      max_end <- NULL
    } else {
      max_end <- min(max_end)
    } 
  } else {
    max_end <- NULL
  }
  if (!is.null(max_end) && max_end > nrow(position_scores)) {
    warning(
      "max_end exceeds the length of the alignment. ",
      "No constraints will be applied."
    )
    max_end <- NULL
  }
  out <- list(
    pattern_id = pattern_id,
    subject_id = subject_id,
    position_scores = position_scores,
    substitution_matrix = submat,
    max_end = max_end
  )
  class(out) <- c("ps", "list")
  return(out)
}
