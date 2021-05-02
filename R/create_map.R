#'
#' @title create_map
#' @description This function creates a map file that is necessary for PLINK.
#'
#' Note change the first entry "data_map[, 1]" to have cromosones in the analysis.
#'
#' @param n is the amount of sequenced enteties.
#' @param path The path where the file is going to be stored.
#'
#' @import data.table
#'
#' @return This function creates a map file

create_map <- function(n, path){
  data_map <- matrix(ncol=4, nrow=n)
  data_map[, 1] <-rep(1,n) # Cromosone code
  data_map[, 2] <- seq(1:n) # Varient identifier
  data_map[, 3] <-rep(0,n) # Dummy = 0
  data_map[, 4] <- seq(1:n) # Base-pair coordiante

  fwrite(as.data.table(data_map), ## writes to locked file
         file = paste(path, "genotypes.map", sep = ""),
         quote = F,
         sep = " ",
         col.names = F)
}
