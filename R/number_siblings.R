#' @title Title
#'
#' @param pheno_file phenofile.txt with path
#'
#' @return integer
#' @export
n_sibs <- function(pheno) {
  return(length(grep("sib", colnames(pheno))) / 3)
}