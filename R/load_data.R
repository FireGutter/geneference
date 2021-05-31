#' @title Import phenotype data
#'
#' @description
#' Import all phenotype-related data generated in simulations.
#'
#' @details
#' Note that this function only loads data about the liabilities and phenotype
#' status of individuals and their family.
#'
#' @param pheno_file path of file with phenotypes, including file extension.
#' Simulation functions save this file as "phenotypes.txt".
#'
#' @return Returns a tibble with the data from \code{pheno_file}.
#'
#' @export
load_phenotypes <- function(pheno_file) {
  stopifnot("pheno_file needs to a valid file with extension '.txt'" =
              (file.exists(pheno_file) && file_ext(pheno_file) == "txt"))
  pheno_ds <- tibble::tibble(data.table::fread(file = pheno_file,
                                               header = TRUE))
  return(pheno_ds)
}

#' @title Import and merge results
#'
#' @description
#' Merge results from analysis performed in PLINK with the true values used by
#' the simulation.
#'
#' @details
#' The function only takes a directory as input and then prompts the user to
#' decide which files in the directory should be loaded. Therefore, simulation
#' data and results from \code{analysis_association()} should reside in the
#' same directory in order to be merged.\cr
#' The names of the columns are created by adding the name of the file as
#' suffix. E.g., the P-values from an analysis named \code{GWAS.assoc} will be
#' in the column named "P_GWAS".\cr
#' True effect sizes are stored in the column "beta". True MAFs of SNPs are
#' stored in the column "MAFs".\cr
#' Note that this function is not supposed to load phenotype information, but
#' rather loads data relating directly to the SNPs. See
#' \code{load_phenotypes()} for reading phenotype related data.\cr
#'
#' @param path path to the folder containing the data, including "/" or "\\\\"
#' at the end. The function can not take relative paths, i.e. the user always
#' needs to specify the full path to the directory.
#'
#' @return A data.table containing the different results from the files.
#'
#' @export
load_results <- function(path) {
  stopifnot("path needs to be a valid path ending with '/' or '\\\\'"
            = (dir.exists(path)
               && (substr(path, nchar(path), nchar(path)) == "/" ||
                     substr(path, nchar(path), nchar(path)) == "\\")))

  vars <- c("P", "BETA", "CHISQ", "SE")
  re <- list.files(path)
  len <- length(re)
  pick <- numeric(len)
  i <- 1
  j <- 1
  
  for(val in re){
    if (any(c(grepl("\\.assoc", val, ignore.case = T), 
              grepl("\\.qassoc", val, ignore.case = T),
              grepl("\\.txt", val, ignore.case = T))) & 
        !grepl("pheno", val, ignore.case = T)) {
      pick[i] <- re[j]
      i <- i + 1
    }
    j <- j + 1
  }
  if(i == 1){
    stop("No analysis files from PLINK found at path.")
  }
  
  pick <- pick[1:(i - 1)]
  
  while(TRUE) {
    cat("Which files would you like to load?: \n")
    nr <- 1
    for(val in pick){
      cat(nr, val)
      cat("\n")
      nr <- nr + 1
    }
    
    cat(nr, "all \n")
    cat(nr + 1, "abort \n")
    
    
    ind = readline('Enter a number or numbers seperated by " ": ')
    vals <- c(strsplit(ind, " ")[[1]])
    if(length(vals) == 1){
      if(vals %in% 1:(nr + 1)){
        break
      }
      else{
        cat("Invalid input. \n\n")
      }
    }
    else {
      for(k in vals){
        if(!(k %in% 1:(nr - 1))){
          flag <- FALSE
          break
        } else {
          flag <- TRUE
        }
      }
      if(!flag) {
        cat("Invalid input. \n\n")
      }
      else if(any(duplicated(vals) | 
                  duplicated(vals, fromLast = TRUE))){
        cat("The same number can not be repeated. \n")
      }
      else {
        break
      }
    }
  }
  
  
  
  if(length(vals) == 1 & vals[1] != nr){
    if(vals == nr + 1){
      stop("This process has been aborted.")
    }
    if (grepl("\\.txt", pick[strtoi(vals)], ignore.case = T)) {
      name <- sub("\\.txt.*", "", pick[strtoi(vals)])
      data <- data.table::fread(paste0(path, pick[strtoi(vals)]),
                                col.names = c(name))
      data <- data %>% tibble::rowid_to_column("SNP")
    }
    else {
      data <- data.table::fread(paste0(path, pick[strtoi(vals)]))
    }
  }
  else{
    if(vals[1] == nr){
      vals <- 1:(strtoi(vals) - 1)
    }
    else{
      vals <- strtoi(vals)
    }
    
    choices <- numeric(length(vals))
    for(i in seq_len(length(vals))) {
      choices[i] <- pick[vals[i]]
    }
    assoc <- grepl("\\.assoc", choices, ignore.case = T)
    qassoc <- grepl("\\.qassoc", choices, ignore.case = T)
    txt <- grepl("\\.txt", choices, ignore.case = T)
    
    assoc_val <- choices[assoc]
    qassoc_val <- choices[qassoc]
    txt_val <- choices[txt]
    
    if(any(assoc)) {
      for(i in seq_len(sum(assoc))) {
        if(i == 1) {
          tmp <- data.table::fread(paste0(path, assoc_val[i]))
          name <- sub("\\.assoc.*", "", assoc_val[i])
          
          data <- data.table::data.table("SNP" = tmp$SNP, "A1" = tmp$A1)
          
          for(v in vars){
            if(v %in% names(tmp)){
              strin <- paste0(v, "_", name)
              stup <- tmp %>% dplyr::select(dplyr::all_of(c("SNP", v)))
              data.table::setnames(stup, v, strin)
              data <- dplyr::inner_join(data, stup, by = "SNP") 
            }
          }
          
        }
        else {
          tmp <- data.table::fread(paste0(path, assoc_val[i]))
          name <- sub("\\.assoc.*", "", assoc_val[i])
          
          for(v in vars){
            if(v %in% names(tmp)){
              strin <- paste0(v, "_", name)
              stup <- tmp %>% dplyr::select(dplyr::all_of(c("SNP", v)))
              data.table::setnames(stup, v, strin)
              data <- dplyr::inner_join(data, stup, by = "SNP") 
            }
          }
        }
      }
    }
    else if(any(qassoc)) {
      for(i in seq_len(sum(qassoc))) {
        if(i == 1) {
          tmp <- data.table::fread(paste0(path, qassoc_val[i]))
          name <- sub("\\.qassoc.*", "", qassoc_val[i])
          
          data <- data.table::data.table("SNP" = tmp$SNP)
          
          for(v in vars){
            if(v %in% names(tmp)){
              strin <- paste0(v, "_", name)
              stup <- tmp %>% dplyr::select(dplyr::all_of(c("SNP", v)))
              data.table::setnames(stup, v, strin)
              data <- dplyr::inner_join(data, stup, by = "SNP") 
            }
          }
        }
        else {
          tmp <- data.table::fread(paste0(path, qassoc_val[i]))
          name <- sub("\\.qassoc.*", "", qassoc_val[i])
          
          for(v in vars){
            if(v %in% names(tmp)){
              strin <- paste0(v, "_", name)
              stup <- tmp %>% dplyr::select(dplyr::all_of(c("SNP", v)))
              data.table::setnames(stup, v, strin)
              data <- dplyr::inner_join(data, stup, by = "SNP") 
            }
          }
        }
      }
    }
    else {
      for(i in seq_len(sum(txt))) {
        if(i == 1) {
          name <- sub("\\.txt.*", "", txt_val[i])
          tmp <- data.table::fread(paste0(path, txt_val[i]), 
                                   col.names = c(name))
          data <- tmp %>% tibble::rowid_to_column("SNP")
        }
        else {
          name <- sub("\\.txt.*", "", txt_val[i])
          tmp <- data.table::fread(paste0(path, txt_val[i]),
                                   col.names = c(name))
          tmp <- tmp %>% tibble::rowid_to_column("SNP")
          data <- dplyr::inner_join(data, tmp, by = "SNP")
        }
      }
    }
    
    
    
    
    if(all(c(any(assoc), any(qassoc), any(txt)))) {
      for(i in seq_len(sum(qassoc))) {
        tmp <- data.table::fread(paste0(path, qassoc_val[i]))
        name <- sub("\\.qassoc.*", "", qassoc_val[i])
        
        for(v in vars){
          if(v %in% names(tmp)){
            strin <- paste0(v, "_", name)
            stup <- tmp %>% dplyr::select(dplyr::all_of(c("SNP", v)))
            data.table::setnames(stup, v, strin)
            data <- dplyr::inner_join(data, stup, by = "SNP") 
          }
        }
      }
      for(i in seq_len(sum(txt))) {
        name <- sub("\\.txt.*", "", txt_val[i])
        tmp <- data.table::fread(paste0(path, txt_val[i]),
                                 col.names = c(name))
        tmp <- tmp %>% tibble::rowid_to_column("SNP")
        data <- dplyr::inner_join(data, tmp, by = "SNP")
      }
    }
    else if (all(c(any(assoc), any(txt)))) {
      for(i in seq_len(sum(txt))) {
        name <- sub("\\.txt.*", "", txt_val[i])
        tmp <- data.table::fread(paste0(path, txt_val[i]),
                                 col.names = c(name))
        tmp <- tmp %>% tibble::rowid_to_column("SNP")
        data <- dplyr::inner_join(data, tmp, by = "SNP")
      }
    }
    else if (all(c(any(assoc), any(qassoc)))) {
      for(i in seq_len(sum(qassoc))) {
        tmp <- data.table::fread(paste0(path, qassoc_val[i]))
        name <- sub("\\.qassoc.*", "", qassoc_val[i])
        
        for(v in vars){
          if(v %in% names(tmp)){
            strin <- paste0(v, "_", name)
            stup <- tmp %>% dplyr::select(dplyr::all_of(c("SNP", v)))
            data.table::setnames(stup, v, strin)
            data <- dplyr::inner_join(data, stup, by = "SNP") 
          }
        }
      }
    }
    else if(all(c(any(qassoc), any(txt)))) {
      for(i in seq_len(sum(txt))) {
        name <- sub("\\.txt.*", "", txt_val[i])
        tmp <- data.table::fread(paste0(path, txt_val[i]),
                                 col.names = c(name))
        tmp <- tmp %>% tibble::rowid_to_column("SNP")
        data <- dplyr::inner_join(data, tmp, by = "SNP")
      }
    }
  }

  return(data)
}

#' @title Augment results for plotting
#'
#' @description
#' Takes data from \code{load_results()} and appends
#' columns used for plotting results.
#'
#' @details
#' This function requires that this data has been loaded with
#' \code{load_results()}. \cr
#' It generates new columns for all columns that represent p-values generated by
#' \code{analysis_association()}. One of the columns, having "_significant" as
#' its suffix, denotes whether or not the p-value is significant at the specified
#' alpha level. The other column, having "_bonferroni" as suffix, denotes whether the p-value is significant
#' with Bonferroni-corrected alpha level.\cr
#' Furthermore, the user is prompted to decide if they want to create a column
#' named "causal", which denotes whether or not a SNP is truly causal for the
#' phenotype status of an individual. The user needs to have loaded "beta.txt"
#' with \code{load_results()} in order to create this column.\cr
#' Lastly, the user is prompted to decide if they want to create a column
#' named "LTFH_transformed", which is the estimated effect sizes of a
#' regression using the LT-FH phenotype transformed to the liability scale.\cr
#' Transforming the effect sizes requires that the user has loaded results from
#' \code{analysis_association()} with both \code{pheno_name = "line_pheno"} and
#' \code{pheno_name = "LTFH_pheno"}. These need to be laoded with \code{load_results()}
#' along with \code{MAFs.txt}.
#'
#' @param data data.table from \code{load_results()}.
#' @param alpha significance level used for hypothesis tests.
#'
#' @return A data.table containing the new appended results.
#'
#' @export
augment_results <- function(data, alpha) {
  plotter <- function(){
    j <- 1
    for(name in names(data)){
      cat(j, name)
      cat("\n")
      j <- j + 1
    }
  }
  
  nr_of_names <- length(names(data))
  
  while(TRUE){
    cat("Do you want to create a causal column from the true effect sizes (by default named beta)")
    ind = readline('Enter 0 for "NO" or 1 for "YES": ')
    if(ind %in% c(0, 1)){
      if(ind == 1){
        while(TRUE) {
          cat("Which column is the beta column?\n")
          plotter()
          cat(nr_of_names + 1, "That column does not exist\n")
          val <- readline('Enter a number: ')
          if(val %in% 1:nr_of_names){
            val <- strtoi(val)
            beta_vals <- dplyr::pull(data, names(data)[val])
            data <- data %>% dplyr::mutate("causal" = beta_vals != 0)
            break
          }
          else if(val == (nr_of_names + 1)){
            cat("Since the beta column does not exist, the causal column can't be created\n")
            ind <- 0
            break
          }
          else {
            cat("Invalid input.\n")
          }
        }
        break
      }
      break
    }
    cat("Invalid input.\n")
  }
  
  m <- nrow(data)
  
  for(val in names(data)) {
    if (grepl("^P_", val, ignore.case = T)) {
      strin1 <- paste0(val, "_significant")
      strin2 <- paste0(val, "_bonferroni")
      numbers <- data %>% dplyr::select(dplyr::all_of(c(val)))
      data <- data %>% dplyr::mutate("strin1v1" = numbers < alpha, 
                                     "strin2v2" = numbers < alpha/m)
      data.table::setnames(data, "strin1v1", strin1)
      data.table::setnames(data, "strin2v2", strin2)
    }
  }
  
  while (TRUE){
    cat("Do you also want to make a transformed LTFH_beta column?")
    ind <- readline('Enter 0 for "NO" or 1 for "YES": ')
    if(ind %in% c(0, 1)){
      if(ind == 1){
        cat("Do you have standard names from the load_results() function?\n")
        sec <-  readline('Enter 0 for "NO" or 1 for "YES": ')
        if(sec %in% c(0, 1)){
          break
        }
      }
      else {
        break
      }
    }
    else {
      cat("Invalid input.\n\n")
    }
  }

  nr_of_names <- length(names(data))
  
  if(ind == 0) {
    return(data)
  }
  else if (ind == 1 & sec == 0) {
    i <- 1
    choices <- numeric(6)
    while(TRUE){
      if(i == 1){
        cat("Which column contain P-values form LTFH?\n")
        plotter()
        cat(nr_of_names + 1, "None of the columns contain P-values from LTFH.\n")
        ind <- readline('Enter a number: ')
        if(ind %in% 1:nr_of_names){
          ind <- strtoi(ind)
          choices[i] <- names(data)[ind]
          i <- i + 1
          LTFH_P <- dplyr::pull(data, names(data)[ind])
        }
        else if(ind == nr_of_names + 1){
          ind <- 0
          break
        }
        else {
          cat("Invalid input.\n")
        }
      }
      else if(i == 2){
        cat("Which column contain P-values from file 'linear'?\n")
        plotter()
        cat(nr_of_names + 1, "None of the columns contain P-values from file 'linear'.\n")
        ind <- readline('Enter a number: ')
        if(ind %in% 1:nr_of_names){
          ind <- strtoi(ind)
          choices[i] <- names(data)[ind]
          i <- i + 1
          LINE_P <- dplyr::pull(data, names(data)[ind])
        }
        else if(ind == nr_of_names + 1){
          ind <- 0
          break
        }
        else {
          cat("Invalid input.\n")
        }
      }
      else if(i == 3){
        cat("Which column contain beta-values form LTFH?\n")
        plotter()
        cat(nr_of_names + 1, "None of the columns contain beta-values from LTFH.\n")
        ind <- readline('Enter a number: ')
        if(ind %in% 1:nr_of_names){
          ind <- strtoi(ind)
          choices[i] <- names(data)[ind]
          i <- i + 1
          LTFH_BETA <- dplyr::pull(data, names(data)[ind])
        }
        else if(ind == nr_of_names + 1){
          ind <- 0
          break
        }
        else {
          cat("Invalid input.\n")
        }
      }
      else if(i == 4){
        cat("Which column contain standard-errors (SE) values from LTFH?\n")
        plotter()
        cat(nr_of_names + 1, "None of the columns contain standard-errors from LTFH.\n")
        ind <- readline('Enter a number": ')
        if(ind %in% 1:nr_of_names){
          ind <- strtoi(ind)
          choices[i] <- names(data)[ind]
          i <- i + 1
          LTFH_SE <- dplyr::pull(data, names(data)[ind])
        }
        else if(ind == nr_of_names + 1){
          ind <- 0
          break
        }
        else {
          cat("Invalid input.\n")
        }
      }
      else if(i == 5){
        cat("Which column contain the MAF values\n")
        plotter()
        cat(nr_of_names + 1, "None of the columns contain the MAF values.\n")
        ind <- readline('Enter a number: ')
        if(ind %in% 1:nr_of_names){
          ind <- strtoi(ind)
          choices[i] <- names(data)[ind]
          i <- i + 1
          MAFS <- dplyr::pull(data, names(data)[ind])
        }
        else if(ind == nr_of_names + 1){
          ind <- 0
          break
        }
        else {
          cat("Invalid input.\n")
        }
      }
      else if(i == 6) {
        ind <- readline('Choose the prevalence parameter "k": ')
        options(digits = 5)
        k <- as.double(ind)
        choices[i] <- ind
        i <- i + 1
      }
      else if(i == 7) {
        while(TRUE) {
          cat("The following data has been entered:\n")
          cat("1: P values from LTFH:", choices[1])
          cat("\n")
          cat("2: P values from file 'linear':", choices[2])
          cat("\n")
          cat("3: BETA values from LTFH:", choices[3])
          cat("\n")
          cat("4: SE from LTFH:", choices[4])
          cat("\n")
          cat("5: MAFs:", choices[5])
          cat("\n")
          cat("6: k (prevalence):", choices[6])
          cat("\n")
          ind = readline('Is this correct? Enter 7 if "YES". Enter any of the above numbers if no: ')
          if(ind %in% 1:7) {
            if(ind == "7"){
              i <- 8
              break
            }
            else{
              i <- strtoi(ind)
              break
            }
          }
          else {
            cat("Invalid input! \n")
          }
        }
      }
      else {
        break
      }
    }
  }
  else {
    stopifnot("data must have columns 'P_LTFH', 'P_linear', 'BETA_LTFH', 'SE_LTFH' and 'MAFs'" =
                all(c('P_LTFH', 'P_linear', 'BETA_LTFH', 'SE_LTFH', 'MAFs') %in% colnames(data)))
    ind <- readline('Choose the prevalence parameter "k": ')
    options(digits = 5)
    k <- as.double(ind)
    sq <- sqrt(m * median(qchisq(data$P_LTFH, 1, lower.tail = T))/median(qchisq(data$P_linear, 1, lower.tail = T)))
    obs <- (data$BETA_LTFH / (data$SE_LTFH * sq)) * sqrt(k * (1 - k) / (2 * data$MAFs * (1 - data$MAFs)))
    
    data <- data %>% dplyr::mutate("LTFH_transformed" = obs)
    
    return(data)
  }
  
  
  sq <- sqrt(m * median(qchisq(LTFH_P, 1, lower.tail = T))/median(qchisq(LINE_P, 1, lower.tail = T)))
  obs <- (LTFH_BETA / (LTFH_SE * sq)) * sqrt(k * (1 - k) / (2 * MAFS * (1 - MAFS)))
  
  data <- data %>% dplyr::mutate("LTFH_transformed" = obs)
  
  return(data)
}
