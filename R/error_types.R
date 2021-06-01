#' @title Calculate error types
#'
#' @description
#' Calculate the number of type 1 (false positive) and type 2 (false negative)
#' errors in an analysis.
#'
#' @details
#' To use this function, the user must run any test with \code{analysis_association()}
#' and use \code{load_results()} to merge this with \code{beta.txt} from the
#' simulation. Then, \code{augment_results()} needs to be run in order to obtain
#' columns for information on which SNPs are causal, significant, and significant
#' with Bonferroni correction. \cr
#' The standard significance level refers to the one used in 
#' \code{augment_results()}.
#'
#' @param causal column in dataset denoting causal SNPs.
#' @param significant column in dataset denoting which SNPs are significant.
#' @param bonferroni column in dataset denoting which SNPs are significant with 
#' Bonferroni correction.
#' @param analysis_name name of the analysis.
#'
#' @return A tibble containing the number of true positives, false positives,
#' false negatives and true negatives with each of the significance levels
#' having their own row.
#'
#' @export
calculate_error_types <- function(causal, significant, bonferroni, analysis_name) {
  stopifnot("analysis_name needs to be a string" = is.character(analysis_name))
  
  # standard significance level
  errors_tbl <- tibble::tibble(
    Analysis = paste(analysis_name,
                     "with standard significance level"),
    true_positives = sum(significant & causal),
    false_positives = sum(significant & !causal),
    false_negatives = sum(!significant & causal),
    true_negatives = sum(!significant & !causal)) %>%
    tibble::add_row(
      Analysis = paste(analysis_name,
                       "with Bonferroni-correction"),
      true_positives = sum(bonferroni & causal),
      false_positives = sum(bonferroni & !causal),
      false_negatives = sum(!bonferroni & causal),
      true_negatives = sum(!bonferroni & !causal))

  return(errors_tbl)
}
