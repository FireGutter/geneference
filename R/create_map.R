#'
#' @title create_map
#' @description As PLINK requires a .map-file in order to run its analysis,
#' we create such a file and fill it with some arbitrary information, since
#' it is not actually used for the analysis.
#'
#' @param n number of sequenced genotypes.
#' @param path location where the file is going to be stored.
#'
#' @importFrom data.table as.data.table fwrite
#'
#' @return Does not return anything, but creates a file with the name
#' "genotypes.map" in the folder specified with 'path'-parameter.

create_map <- function(n, path) {
  data_map <- matrix(ncol = 4, nrow = n)
  data_map[, 1] <- rep(1, n) # Chromosome code
  data_map[, 2] <- seq(1:n) # Variant identifier
  data_map[, 3] <- rep(0, n) # Dummy with value 0
  data_map[, 4] <- seq(1:n) # Base-pair coordinate

  fwrite(as.data.table(data_map), ## writes to locked file
         file = paste0(path, "genotypes.map", sep = ""),
         quote = F,
         sep = " ",
         col.names = F)
}
