
how_many_siblings_are_in_this_easily_openable_txt_file <- function(pheno_file){
  return(length(grep("sib", colnames(data.table::fread(pheno_file)))) / 3)
}