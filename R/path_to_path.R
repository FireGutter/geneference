#' @title Quote path
#'
#' @param path a path string \code{path//to//file}
#'
#' @return a path string in quotes \code{\"path//to//file\"}
#' @export
#'
#' @examples
path_to_path <- function(path) {
  if (path != "") {
    return(paste0("\"", path, "\""))
  } else {
    return(path)
  }
}