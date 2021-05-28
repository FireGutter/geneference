#' @title Quote path
#' @noRd
#'
#' @description takes a path string and puts it in quotes such that shell
#' doesn't get confused.
#'
#' @param path a path string \code{path//to//file}.
#'
#' @return a path string in quotes \code{\"path//to//file\"}.
quote_path <- function(path) {
  if (path != "") {
    return(paste0("\"", path, "\""))
  } else {
    return(path)
  }
}
