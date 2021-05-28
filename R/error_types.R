#' @title Calculate the number of error types in an analysis
#'
#' @description description to be written.
#'
#' @param dataset data loaded with \code{load_assoc_results()}.
#' @param alpha significance level used for hypothesis tests.
#' @param analysis_name name of analysis.
#'
#' @return A tibble with rows... columns...
#'
#' @export
calculate_error_types <- function(dataset, alpha, analysis_name) {
  stopifnot("dataset should have columns 'significant', 'causal' and 'bonferroni'" =
              all(c("significant", "causal", "bonferroni") %in% colnames(dataset)),
            "alpha needs to be a number between 0 and 1" =
              (is.numeric(alpha) && 0 < alpha && alpha < 1 &&
                 length(alpha) == 1),
            "analysis_name needs to be a string" = is.character(analysis_name))
  
  # standard significance level
  errors_tbl <- tibble::tibble(
    Analysis = paste(analysis_name,
                     "with significance",
                     alpha),
    true_positives = sum(dataset$significant & dataset$causal),
    false_positives = sum(dataset$significant & !dataset$causal),
    false_negatives = sum(!dataset$significant & dataset$causal),
    true_negatives = sum(!dataset$significant & !dataset$causal)) %>%
    tibble::add_row(
      Analysis = paste(analysis_name,
                       "with Bonferroni-correction"),
      true_positives = sum(dataset$bonferroni & dataset$causal),
      false_positives = sum(dataset$bonferroni & !dataset$causal),
      false_negatives = sum(!dataset$bonferroni & dataset$causal),
      true_negatives = sum(!dataset$bonferroni & !dataset$causal))

  return(errors_tbl)
}
