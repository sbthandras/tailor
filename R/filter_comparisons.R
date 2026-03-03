#' Filter comparisons
#'
#' Filter comparisons (brakpoints or adapters) based on ids or properties from
#' the supporting data set.
#' @importFrom dplyr distinct
#' @param comps data.frame; the comparisons data set to be filtered. Must be of
#' class "breakpoints" or "adapter". See \code{\link{find_breakpoints}},
#' \code{\link{find_adapter}} or \code{\link{find_all_adapters}} for compatible
#' data frames.
#' @param filter_value character; the value to filter by. By default, this is
#' a pattern_id or subject id from the comparisons data frame, but can be any
#' value in the supporting data set if filter_by and data are specified.
#' @param filter_by character; name of the variable in the supporting data set
#' to filter by. By default, this is "Core_ORF", but can be any column in the
#' supporting data set.
#' @param data data.frame; the supporting data set that links ids with
#' properties. This is required if filter_by is not "Core_ORF". Must contain a
#' column with the same ids as the comparisons data frame (e.g. "Core_ORF") and
#' a column with the values to filter by (e.g. "ACL 2").
#' @examples
#' data(rbps)
#' data(adapters)
#' # convert to pident matrix and order by the original data frame
#' amat <- adapter_matrix(adapters, ids = rbps$Core_ORF)
#' # cluster RBPs into 2 to 5 clusters
#' clusters <- cluster_adapters(amat, k_min = 2, k_max = 5)
#' rbps <- dplyr::left_join(rbps, clusters, by = c("Core_ORF" = "id"))
#' # adapters with "ON513429-1"
#' adapters |> filter_comparisons("ON513429-1")
#' # adapter between "MN395291-1" and "ON513429-1"
#' adapters |> filter_comparisons("MN395291-1") |> filter_comparisons("ON513429-1")
#' # adapters in ACL 2
#' adapters |> filter_comparisons("ACL 2", filter_by = "cluster", data = rbps)
#' @export
filter_comparisons <- function(
    comps,
    filter_value,
    filter_by = "Core_ORF",
    data = NULL
  ) {
  # check that comps has the right class
  if (!inherits(comps, c("breakpoints", "adapter"))) {
    msg <- paste0(
      "'comps' must be a data frame of class 'breakpoints' or 'adapter'. ",
      "Use find_breakpoints() or find_adapter() to create a compatible data frame."
    )
    stop(msg)
  }
  if (is.null(data)) {
    out <- comps |>
      dplyr::filter(
        !!rlang::sym("pattern_id") %in% filter_value |
        !!rlang::sym("subject_id") %in% filter_value
      )
    return(out)
  }
  # check comps attributes for conversion type -> automatic ID conversion
  if (filter_by %in% names(data) == FALSE) {
    stop("Filter column not found")
  }
  if (filter_value %in% data[[filter_by]] == FALSE) {
    stop("Filter value not found")
  }
  if (!exists(as.character(substitute(comps)))) {
    stop("Conserved N-terminal database not found")
  }
  if (!exists(as.character(substitute(data)))) {
    stop("Protein database not found")
  }
  comps <- dplyr::distinct(comps) #remove once dps rows are truly unique
  if (filter_by == "Core_ORF") {
    ipat <- which(comps$pattern_id == filter_value)
    isub <- which(comps$subject_id == filter_value)
  }
  if (filter_by != "Core_ORF") {
    # breakpoint pattern and subject id-s are in Core_ORF.
    ids <- data$Core_ORF[data[[filter_by]] == filter_value]
    ipat <- which(comps$pattern_id %in% ids)
    isub <- which(comps$subject_id %in% ids)
  }
  index <- sort(unique(c(ipat, isub)))
  return(comps[index,])
}
