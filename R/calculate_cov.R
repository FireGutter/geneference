#'
#' @title Calculate covariances between liabilities
#'
#' @description Compute a sample covariance matrix. This can be used to
#' validate that covariances of data are as expected.
#'
#'
#' @param file path to the file that should be read.
#' @param sibs number of siblings. optional. ÆNDRERS
#'
#' @import stats
#'
#' @return A covariance matrix using the data from the selected file. All
#' entries in the matrix are multiplied by 100.
#'
#' @export


calculate_cov <- function(file, sibs) {
  stopifnot("sib needs to be a non-negative integer" =
            (!exists("sibs") || (sibs >= 0 && class(sibs) == "numeric" && round(sibs) == sibs)))
  
  ph <- data.table::fread(file)
  ph <- as.data.frame(ph)
  
  if (missing(sibs)){  # if sibs isn't specified we just assign the max number of sibs.
    sibs <- n_sibs(ph)
  }
  
  ph[ph == -9] <- NA
  indexes <- c(c(4:5, 8), c(seq(11, 11 + sibs * 3, by = 3)))
  round(100 * cov(ph[, indexes], use = "complete.obs"))
}
