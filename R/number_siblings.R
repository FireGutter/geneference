#' @title Title
#'
#' @param pheno_file phenofile.txt with path
#'
#' @return integer
#' @export
n_sibs <- function(pheno_file) {
  return(length(grep("sib", colnames(data.table::fread(pheno_file)))) / 3)
}