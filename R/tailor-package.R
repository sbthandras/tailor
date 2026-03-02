#' Receptor Binding Proteins
#'
#' This data set contains example bacteriophage receptor binding proteins (RBPs) 
#' for demonstration purposes.
#'
#' @format A data frame with 20 rows and 2 variables:
#' \describe{
#'   \item{Core_ORF}{Unique identifier for each receptor binding protein}
#'   \item{translation}{Amino acid sequence of the protein}
#' }
"rbps"

#' Adapters
#'
#' This data set contains identified adapters from all pairwise comparisons of 
#' the RBPs in the `rbps` data set. The approach to calculate these adapters
#' was the following: first, check for adapters with the "xbar.one" method, 
#' because it seems to be highly specific (i.e. it rarely identifies adapters 
#' when there isn't one). If an adapter was found, then the "cemean" method was 
#' used to get a more accurate estimate of the adapter length.
#' 
#' @format A data frame with 90 rows and 6 variables:
#' \describe{
#'   \item{pattern_id}{pattern_id}
#'   \item{subject_id}{subject_id}
#'   \item{start}{location of the first amino acid within the adapter sequence}
#'   \item{end}{location of the last amino acid within the adapter sequence}
#'   \item{mean_score}{mean amino acid position score between 'start' and 'end'}
#'   \item{pident}{ratio of identical amino acids between 'start' and 'end'}
#' }
"adapters"


