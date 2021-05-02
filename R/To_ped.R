#'
#' @title to_ped
#' @description turns a binary matrix into a PLINK-friendly format
#'
#' @param df input of binary matrix
#' @param part is the amount of parts, any given function splits the data in
#'
#' @return This function retruns a n x m matrix

#This function should be looked at again to make it smarter
to_ped <- function(df, part){
  n <- dim(df)[1]
  m1 <- dim(df)[2]
  m2 <- m1 * 2 + 6
  tmp <- matrix(NA, 3, 2)
  tmp[ ,1] <- matrix(c("a", "c", "c"), 1, 3)
  tmp[ ,2] <- matrix(c("a", "a", "c"), 1, 3)
  tmp2 <- matrix(NA, n, m2)
  tmp2[ ,2:6] <- 1
  for (i in 1:dim(df)[1]) {
    tmp2[i,1] <- i + part
    tmp2[i, 7:m2] <- as.vector(t(tmp[unlist(df[i, ]) + 1, ][1:m1, ]))
  }
  return(tmp2)
}
