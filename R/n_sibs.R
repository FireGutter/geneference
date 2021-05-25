#' @title Find number of siblings in largest family
#'
#' @param pheno tibble loaded with load_pheno().
#'
#' @return Number of siblings in the largest family.
n_sibs <- function(pheno) {
  return(length(grep("sib", colnames(pheno))) / 3)
}
