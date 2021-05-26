#'
#' @title sim_no_family
#'
#' @description This function simulate data for individual enteties (the data for the parents)
#' is unknown.
#'
#' @param n is the amount of enteties (humans)
#' @param m is the amount of SNPS per entety
#' @param q is the amount of causal SNPs
#' @param hsq heritability parameter that is squard
#' @param k is the prevalence of trait parameter
#' @param path is the file path, where the files will be stored
#'
#' @import stats
#'
#' @return This function returns five files: Three txt files: Beta, MAFs and phenotypes, a genotypes MAP file and a genotypes PED file.
#'
#' @export



sim_no_family <- function(n, m, q, hsq, k, path){
  stopifnot("n needs to be an integer greater than 0" =
              (n > 0 && class(n) == "numeric" && n == round(n)),
            "m needs to be an integer greater than 0" =
              (m > 0 && class(m) == "numeric" && m == round(m)),
            "q needs to be an integer greater than 0 and smaller than m" =
              (q > 0 && class(q) == "numeric" && q == round(q) && q <= m),
            "hsq needs to be a number between 0 and 1" =
              (hsq > 0 && hsq < 1 && class(hsq) == "numeric"),
            "k needs to be a number between 0 and 1" =
              (k > 0 && k < 1 && class(k) == "numeric"),
            "path needs to be default or a valid path ending with '/' or '\\\\'"
            = (path == "" || (dir.exists(path))
               && (substr(path, nchar(path), nchar(path)) == "/" ||
                     substr(path, nchar(path), nchar(path)) == "\\")))
  
  path = path_validation(path)

  #Function that ties everything together (one function to rule them all)

  #Look at function to determine how to best divide the overall simulation into chuncks

  #library(tidyverse)
  #library(parallel)
  #library(data.table)
  # there are many packages that allows for parallel processing in R. Future is a well implemented one in my opinion
  #library(future)
  #flock: short for file lock. Allows multiple R sessions to write to the same file. (alternatively make separate files and merge later)
  #library(flock)

  future::plan(future::multisession(workers = 2))

  SNP_simulation <- function(n, m, MAFs){
    #Function that draws from a binomial distribution and uses the MAF vector of probabilities
    #Returns a list of length n with entires m integers.
    return(lapply(1:n, function(y) rbinom(m, 2, MAFs)))
  }


  simulation <- function(n, m, MAFs){
    #Detecting the number of cores of ones computer and making one core available
    cores <- parallel::detectCores(logical = FALSE)-1

    #If the cores are creater than 1 then we check if n is a multiple of cores
    if (cores>1){
      while (n %% cores != 0){
        cores <- cores-1
      }
    } else {
      cores <- 1
    }

    #Making a cluster of the cores
    cl <- parallel::makeCluster(cores)
    #Exporting the cores to the function we want them to run on
    parallel::clusterExport(cl, c("SNP_simulation"))

    #Making the multiprocessing simulation (note: because of the way rbindlist works we have to switch rows and columns when inputting it to 'SNP_simulation'). Can maybe also use bind_cols() instead
    geno <- parallel::clusterCall(cl, function(persons, SNPs, MAF) SNP_simulation(persons, SNPs, MAFs), persons = (n/cores), SNPs = m, MAF = MAFs)

    parallel::stopCluster(cl)

    #binding it all together
    return(t(dplyr::bind_cols(geno)))
  }

  create_map(n, path)

  #Calculating the number MAF probabilities so that each core has the same probabilities
  MAFs <- runif(m, 0.01, 0.49)

  parts <- ceiling((n*m)/10000000) #varaiable to control how many parts we divide the problem into

  data.table::fwrite(data.table::as.data.table(MAFs),
         paste(path,"MAFs.txt",sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  causual_SNP <- sample.int(m, size = q, replace = F)

  beta <- matrix(0, nrow = m, ncol = 1)
  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq/q))

  data.table::fwrite(data.table::as.data.table(beta),
         paste(path,"beta.txt",sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  data.table::fwrite(data.table::as.data.table(rbind(c("FID", "IID", "pheno", "line_pheno"))),
         paste(path,"phenotypes.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  lock = tempfile()

  future.apply::future_lapply(1:parts, function(i) {
    #for (i in 1:parts) {
    if (i == parts & n%%parts != 0){
      batch_size = n%%parts
    } else {
      batch_size = ceiling(n/parts)
    }

    locked = flock::lock(lock) # locks file

    persons <- simulation(batch_size, m, MAFs)

    mu <- 2*MAFs
    sigma <- sqrt(2*MAFs*(1-MAFs))

    pliab <- sweep(sweep(persons, 2, mu, FUN = "-"), 2, sigma, FUN = "/")

    liability <- pliab %*% beta + rnorm(batch_size, 0, sqrt(1-hsq))

    critical <- qnorm(1-k)

    y <- sapply(liability, function(x) ifelse(x>critical, 2, 1))
    line_pheno <- y + 1

    id <- matrix(1:(batch_size))+(batch_size)*(i-1)

    data.table::fwrite(data.table::as.data.table(to_ped(persons, (i-1) * (batch_size))), ## writes to locked file
           paste(path,"genotypes.ped", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    data.table::fwrite(data.table::as.data.table(cbind(id,rep(1,batch_size), y, line_pheno)),
           paste(path,"phenotypes.txt",sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
