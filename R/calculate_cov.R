#' @title Calculate covariances between liabilities
#'
#' @description Compute a sample covariance matrix. This can be used to
#' validate that covariances of data are as expected.
#'
#'
#' @param pheno dataset with phenotypes, as loaded by \code{load_phenotypes()}.
#' @param sibs \emph{optional:} number of siblings to be included in covariance
#' matrix. Defaults to using all of an individual's siblings.
#'
#' @return A covariance matrix using the data from the given dataset.
#'
#' @import stats
#'
#' @export
calculate_cov <- function(pheno, sibs) {
  stopifnot("sibs needs to be a non-negative integer" =
            (missing(sibs) ||
               (sibs >= 0 && is.numeric(sibs) && round(sibs) == sibs)))
  stopifnot("sibs is greater than the number of siblings in pheno" =
              sibs <= n_sibs(pheno))

  # If sibs isn't specified we just assign the max number of siblings in data
  if (missing(sibs)) {
    sibs <- n_sibs(pheno)
  }

  pheno[pheno == -9] <- NA
  indexes <- c(c(4:5, 8), c(seq(11, 11 + sibs * 3, by = 3)))

  # This informs the user about the distribution of the siblings
  message(noquote("Attention"))
  for (i in 0:sibs) {
    message(noquote(
      paste0("number of families with at least ", i, " siblings: ",
             sum(complete.cases(pheno[, c(seq(11, 11 + i * 3, by = 3))])))))
  }
  return(round(cov(pheno[, indexes], use = "pairwise.complete.obs"), 2))
}
