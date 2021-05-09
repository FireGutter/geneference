#'
#' @title Test covariances between liabilities
#'
#' @description Compute a sample covariance matrix. This can be used to
#' validate that covariances of data are as expected.
#'
#'
#' @param file path to the file that should be read.
#' @param sib number of siblings. Default is 0.
#'
#' @importFrom data.table as.data.table fread
#' @importFrom stats cov
#'
#' @return A covariance matrix using the data from the selected file. All
#' entries in the matrix are multiplied by 100.
#'
#' @export


cov_tester <- function(file, sib = 0) {
  ph <- fread(file)
  ph <- as.data.frame(ph)
  indexes <- c(c(4:5, 8), c(seq(11, 11 + sib * 3, by = 3)))
  round(100 * cov(ph[, indexes]))
}
