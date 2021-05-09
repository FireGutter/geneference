#'
#' @title Format conversion
#' @description Converts our format to PLINK's .bed-format.
#'
#' @param df genotypes in our simulated format.
#' @param part amount of parts any given function splits the data in.
#'
#' @return Returns a n x (6 + 2 * m) matrix so it follows PLINK's .ped format.

to_ped <- function(df, part) {
  n <- dim(df)[1]
  m <- dim(df)[2]
  m_ped_out <- m * 2 + 6
  base_pair_convert <- matrix(NA, 3, 2)
  base_pair_convert[, 1] <- matrix(c("a", "c", "c"), 1, 3)
  base_pair_convert[, 2] <- matrix(c("a", "a", "c"), 1, 3)
  out_matrix <- matrix(NA, n, m_ped_out)
  out_matrix[, 2:6] <- 1
  for (i in 1:n) {
    out_matrix[i, 1] <- i + part
    out_matrix[i, 7:m_ped_out] <- as.vector(
      t(base_pair_convert[unlist(df[i, ]) + 1, ][1:m, ]))
  }
  return(out_matrix)
}
