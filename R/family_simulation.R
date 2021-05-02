#'
#' @title family_simulation
#'
#' @description This function simulate data for individual enteties with specefic family history. both the liability and
#' pheno types of both siblings and parents are known to the subjects.
#'
#' @param n is the amount of enteties (humans)
#' @param m is the amount of SNPS per entety
#' @param q is the amount of causal SNPs
#' @param hsq heritability parameter that is squard
#' @param k is the prevalence of trait parameter
#' @param path is the file path, where the files will be stored. If nothing is added, the function will place the function
#' at the place of the working directory. Use getwd() to see your working directory
#' @param sib is number of siblings per individual.
#'
#' @importFrom data.table data.table as.data.table fwrite :=
#' @import future.apply
#' @import flock
#' @import dplyr
#' @import future
#' @importFrom stats rbinom
#'
#' @return This function returns five files: Three txt files: Beta, MAFs and phenotypes, a genotypes MAP file and a genotypes PED file.
#'
#' @export


family_simulation <- function(n, m, q, hsq, k, path = "", sib = 0){

  # plan(multisession(workers = availableCores(logical = F) - 1)) # Makeing parallel computational session

  parentmaker <- function(m, antal, MAFs){
    future.apply::future_sapply(1:antal, function(i){rbinom(m, 2, MAFs)}, future.seed = T)
  }
  # Make a function that created the SNPs for the parents

  parts <- ceiling((n*m)/10000000) # Here we hard code the function to simulate 10 million SNPs per session

  splits <- c(rep(ceiling(n / parts), parts - 1), n - sum(rep(ceiling(n / parts), parts - 1))) # Splitting the work on the amount of workers

  cusplits <- c(0, cumsum(splits)) # Making ids for the different cores.


  create_map(m, path) # Generate MAP-file

  MAFs <- runif(m, 0.01, 0.49) # Find MAFs for children and parents

  fwrite(as.data.table(MAFs), #Save the MAFs
         paste0(path, "MAFs.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  mu <- 2*MAFs # Find mean per SNP
  sigma <- sqrt(2*MAFs*(1-MAFs)) # Find standard diviation per SNP
  critical <- qnorm(1-k) # Calculate the critical level

  causual_SNP <- sample.int(m, size = q, replace = F) # Create causal SNPs

  beta <- matrix(0, nrow = m, ncol = 1) # Make the beta matrix

  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq/q))


  fwrite(as.data.table(beta), #Save the beta
         paste0(path, "beta.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  header <- c("FID", "IID", "pheno", "child_lg", "child_liab", "par1_pheno", "par1_lg", "par1_liab", "par2_pheno", "par2_lg", "par2_liab")

  if(sib != 0){
    sib_header <- numeric(3*sib)
    for(i in 0:(sib-1)){
      sib_header[3*i + 1] <- paste0("sib", i + 1, "_pheno", sep = "")
      sib_header[3*i + 2] <- paste0("sib", i + 1, "_lg", sep = "")
      sib_header[3*i + 3] <- paste0("sib", i + 1, "_liab", sep = "")
    }
  }
  header <- c(header, sib_header)
  # Create the header for the phenofile:
  fwrite(as.data.table(rbind(header)),
         paste0(path, "phenotypes.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  lock = tempfile()

  future.apply::future_lapply(1:parts, function(i) {

    parentmatrix <- parentmaker(m = m, antal = 2*splits[i], MAFs) # Create the parents

    child <- t(future.apply::future_sapply(seq(1, ncol(parentmatrix), 2), function(j){ # Using the parents SNPs the childrens SNPs can be calculated
      child <- rowSums(parentmatrix[,j:(j + 1), drop = FALSE])/2
      round(child + runif(m, min = -0.0001, max = 0.0001), 0)}, future.seed = T))




    if(sib != 0){

      sibtable <- data.table("begin" = numeric(splits[i]))

      for(k in 1:sib){
        sibs <- t(future.apply::future_sapply(seq(1, ncol(parentmatrix), 2), function(j){ # In the same way as the children SNPs are created we can create the siblings SNPs
          sibs <- rowSums(parentmatrix[,j:(j + 1), drop = FALSE])/2
          round(sibs + runif(m, min = -0.0001, max = 0.0001), 0)}, future.seed = T))

        sib_lg <- sweep(sweep(sibs, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta # Find the genetic liability, liability and phenotype for the siblings and add this to the "sibtable"
        sib_liab <- sib_lg + rnorm(splits[i], 0, sqrt(1 - hsq))
        sib_pheno <- sapply(sib_liab, function(x) ifelse(x > critical, 2, 1))


        sibtable <- cbind(sibtable, sib_pheno, sib_lg, sib_liab)


      }

      sibtable <- sibtable[,"begin" := NULL] # Remove the begin column in the table.

    }



    parlg <- sweep(sweep(t(parentmatrix), 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta
    parliab <- parlg + rnorm(2*splits[i], 0, sqrt(1 - hsq))

    c_lg <- sweep(sweep(child, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta

    c_liab <- c_lg + rnorm(splits[i], 0, sqrt(1 - hsq))



    # Make phenotypes:
    c_pheno <- sapply(c_liab, function(x) ifelse(x > critical, 2, 1))
    p_pheno <- sapply(parliab, function(x) ifelse(x > critical, 2, 1))


    # FID for the children/parents:

    id <- matrix(c((cusplits[i]+1):cusplits[i+1])) # Create the ID per individual


    locked = flock::lock(lock) # locks file

    fwrite(as.data.table(to_ped(child, part = cusplits[i])), # writes to locked file
           paste0(path, "genotypes.ped", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)


    fwrite(as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab, p_pheno[seq(1, 2*splits[i], 2)], parlg[seq(1, 2*splits[i], 2)], parliab[seq(1, 2*splits[i], 2)], p_pheno[seq(2, 2*splits[i], 2)], parlg[seq(2, 2*splits[i], 2)], parliab[seq(2, 2*splits[i], 2)], sibtable)),
           paste0(path, "phenotypes.txt", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)
    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
