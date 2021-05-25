#'
#' @title Check if certain files exists at a path
#'
#' @description This function check if certain files exists at the path that is
#' given to any of the simulation functions. If these files exsists, the user
#' has different options: 1. Stopping the simulation, 2. Over-writing the files
#' (deleting them and creating new ones), 3. Choose a new path for the files
#' and 4. Make a new folder at the given path, for the new files. If the user
#' does not give valid input the function will re-ask the user.
#'
#'
#' @param path for validation
#'
#' @return A choice of a valid path, new folder, over-writing the files or
#' stopping the simulation.
#'
#' @export

path_validation <- function(path){
  while(any(file.exists(paste0(path, "phenotypes.txt")), 
            file.exists(paste0(path, "beta.txt")),
            file.exists(paste0(path, "genotypes.ped")), 
            file.exists(paste0(path, "genotypes.map")),
            file.exists(paste0(path, "MAFs.txt")))) {
    
    cat("Some or more files exists on the path given.\n")
    cat("What decision would you like? \n\n")
    cat("1: Stop the simulation\n")
    cat("2: Over-write the files\n")
    cat("3: Choose new path\n")
    cat("4: Make a new folder at path for the new files\n")
    cat("\n")
    ind = readline("Enter one number or skip to stop the simulation: ")
    
    if (ind %in% c("1", "2", "3", "4", "")) {
      if (ind == "1" | ind == "") {
        stop("You chose to stop the simulation.")
      }
      else if (ind == "2") {
        cat("Over-writing the files.")
        filenames <- c("phenotypes.txt", "beta.txt", "genotypes.ped", 
                       "genotypes.map", "/MAFs.txt")
        
        for (fil in filenames) {
          if(file.exists(paste0(path, fil))) {
            unlink(paste0(path, fil))
          }
        }
      }
      else if (ind == "3") {
        path = readline("Choose your new path (full-path): ")
        while(!(dir.exists(path)) | !(substr(
          path, nchar(path), nchar(path)) %in% c("/", "\\"))){
          if (!(substr(path, nchar(path), nchar(path)) %in% c("/", "\\"))) {
            cat("Path needs to be default or a valid path ending with '/'
                or '\\\\'")
          }
          else {
            cat("Path does not exist choose new one:\n")
          }
          path = readline("Choose your new path: ")
        }
      }
      else {
        folder_name = readline("Choose your new folder name: ")
        while(dir.exists(paste0(
          path, folder_name, substr(path, nchar(path), nchar(path))))){
          
          cat("Folder does already exists. Choose a new name:\n")
          folder_name = readline("Choose new folder name: ")
        }
        dir.create(paste0(path, folder_name))
        path = paste0(path, folder_name, substr(path, nchar(path), nchar(path)))
      }
    }
    else {
      cat("Not valid input.\n")
      cat("\n")
    }
  }
  
  return(path)
}