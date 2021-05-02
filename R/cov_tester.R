#'
#' @title cov_tester
#'
#' @description This function simulate data for individual enteties (the data for the parents)
#' is unknown.
#'
#' @param read the path to the file that should be read
#' @param sib is the amount of siblings
#'
#' @import data.table
#'
#' @return A covariance matrix using th data from the selected file
#'
#' @export


cov_tester <- function(read, sib = 0){
  ph <- fread(read)
  index <- c(c(4:5, 8), c(seq(11, 11 + sib*3, by = 3)))
  #index <- c(4:5, 8, 11, 14, 17, 20)
  round(100*cov(ph[, ..index]))
}

