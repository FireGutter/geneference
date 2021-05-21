#'
#' @title Simulate genotypes with family history
#' @md
#' @description Simulate genotypes for individuals with family history, where
#' each individual has an arbitrary number of siblings. Parents' genotypes are
#' simulated and used for simulating the genotypes of individuals and their
#' siblings.
#' \code{family_dist_simulaiton} makes use of parallel computation in order to
#' decrease the running time.
#'
#' @section Warning:
#' Simulating large datasets lead to very large files (eg. in the order of 100k
#' genotypes with 100k SNPs takes up approximately 40GB of space).
#' Please ensure that you have sufficient disk space at the disk where
#' \code{path} resides prior to running the simulation.
#' See FUNCTION-NAME DOCUMENTATION REFERENCE(this is a placeholder until
#' function is implemented) in order to convert the generated .ped to the
#' smaller .bed format.
#'
#' @param n number of genotypes (individuals). Given as a vector that match the
#' dist vector (see discribtion below).
#' @param m number of SNPS per genotype.
#' @param q number of causal SNPs, i.e. SNPs that effect chances of having
#' the phenotype.
#' @param hsq squared heritability parameter.
#' @param k prevalence of phenotype.
#' @param dist the distribution of siblings. Given as a vector with the same length
#' as the "n"-vector.
#' @param path directory where the files will be stored. If nothing is
#' specified, \code{family_simulation} writes its files in the current
#' working directory.
#' @details Note that dist denotes the amount of siblings matched with the amount of
#' individuals n. E.g. n = c(100, 200, 300, 400) and dist = c(0, 1, 2, 3) would mean
#' that we want a total of 100 + 200 + 300 + 400 = 1000 individuals. Here 100
#' individuals have 0 siblings, 200 has 1 sibling and so on...
#'
#' @importFrom data.table data.table as.data.table fwrite :=
#' @import future.apply
#' @import flock
#' @import dplyr
#' @import future
#' @importFrom stats rbinom rmultinom rpois start
#'
#' @return Does not return any value, but prints the following five files to
#' the \code{path} parameter specified in the function call:
#' * Three text files:
#'     * beta.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true effect of the i'th SNP.
#'     * MAFs.txt - a file of \code{m} rows with one column. The i'th row is
#'     the true Minor Allelle Frequence of the i'th SNP.
#'     * phenotypes.txt - a file of \code{n} rows, number of columns depend on
#'     number of siblings. Note that some individuals has fewer siblings as
#'     others. Thus we fill the empty value-slots with "-9". The file
#'     contains the phenotype and liability of
#'     each individual as well as information on the liabilities and phenotype
#'     status of their parents and siblings.
#' * genotypes.map - a file created such that PLINK will work with the genotype
#' data.
#' * genotypes.ped - the simulated genotypes in a PLINK-readable format.
#'
#' @export

family_dist_simulaiton <- function(n, m, q, hsq, k, dist, path = ""){
  
  # Set worker nodes:
  plan(multiprocess, workers = max(availableCores(logical = F) - 1, 1))

  parentmaker <- function(m, antal, MAFs){
    sapply(1:antal, function(i){rbinom(m, 2, MAFs)})
  }

  largest_sib <- max(dist)
  values <- numeric(length(dist))

  for (i in 1:length(dist)){
    values[i] <- ceiling((n[i]*m)/10000000)
  }


  parts <- sum(values)
  splits <- numeric(parts)
  amount_of_sibs <- splits


  lower_bound <- 1
  upper_bound <- 0

  for(i in 1:length(values)){
    x <- values[i]
    pers <- n[i]
    upper_bound <- upper_bound + x
    val <- rep(ceiling(pers / x), x - 1)
    splits[lower_bound:upper_bound] <- c(val, pers - sum(val))
    amount_of_sibs[lower_bound:upper_bound] <- rep(dist[i], upper_bound + 1 - lower_bound)
    lower_bound <- lower_bound + x
  }

  cusplits <- c(0, cumsum(splits))

  create_map(m, path) #Generate MAP-file

  MAFs <- runif(m, 0.01, 0.49) #Find MAFs for children and parents

  fwrite(as.data.table(MAFs), #Save the MAFs
         paste0(path, "MAFs.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)

  mu <- 2 * MAFs
  sigma <- sqrt(2 * MAFs * (1 - MAFs))
  critical <- qnorm(1 - k)

  causual_SNP <- sample.int(m, size = q, replace = F) #Create causal SNPs

  beta <- matrix(0, nrow = m, ncol = 1) #Make the beta matrix

  beta[causual_SNP] <- rnorm(q, 0, sqrt(hsq/q))


  fwrite(as.data.table(beta), #Save the beta
         paste0(path, "beta.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  header <- c("FID", "IID", "pheno", "child_lg", "child_liab", "par1_pheno", "par1_lg", "par1_liab", "par2_pheno", "par2_lg", "par2_liab")

  if(largest_sib != 0){
    sib_header <- numeric(3*largest_sib)
    for(i in 0:(largest_sib - 1)){
      sib_header[3 * i + 1] <- paste0("sib", i + 1, "_pheno", sep = "")
      sib_header[3 * i + 2] <- paste0("sib", i + 1, "_lg", sep = "")
      sib_header[3 * i + 3] <- paste0("sib", i + 1, "_liab", sep = "")
    }
    header <- c(header, sib_header)
  }

  #We create the header for the phenofile:
  fwrite(as.data.table(rbind(header)),
         paste0(path, "phenotypes.txt", sep = ""),
         quote = F,
         sep = " ",
         col.names = F,
         append = T)


  lock = tempfile()

  future.apply::future_lapply(1:parts, function(i) {

    parentmatrix <- parentmaker(m = m, antal = 2*splits[i], MAFs)

    child <- t(sapply(seq(1, ncol(parentmatrix), 2), function(i){
      child_mean <- rowMeans(parentmatrix[,i:(i+1), drop = FALSE])
      round_vec <- rbinom(n = m, 1, 1/2)
      snps_for_child <- rbinom(n = m, 2, 1/2)
      vals <- ifelse(child_mean == 1, 1, 0)
      child_mean[vals] <- ifelse(parentmatrix[,i][vals] == 1, snps_for_child, 1)
      return(if_else(round_vec == 1, ceiling(child_mean), floor(child_mean)))}))



    if(largest_sib != 0){

      sibtable <- data.table("start" = numeric(splits[i]))

      if (amount_of_sibs[i] != 0) {
        for(k in 1:amount_of_sibs[i]){
          sibs <- t(sapply(seq(1, ncol(parentmatrix), 2), function(i){
            sibs <- rowMeans(parentmatrix[,i:(i+1), drop = FALSE])
            round_vec <- rbinom(n = m, 1, 1/2)
            snps_for_sibs <- rbinom(n = m, 2, 1/2)
            vals <- ifelse(sibs == 1, 1, 0)
            sibs[vals] <- ifelse(parentmatrix[,i][vals] == 1, snps_for_sibs, 1)
            return(if_else(round_vec == 1, ceiling(sibs), floor(sibs)))}))

          sib_lg <- sweep(sweep(sibs, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta
          sib_liab <- sib_lg + rnorm(splits[i], 0, sqrt(1 - hsq))
          sib_pheno <- sapply(sib_liab, function(x) ifelse(x > critical, 2, 1))


          sibtable <- cbind(sibtable, sib_pheno, sib_lg, sib_liab)
        }}

      if (largest_sib != amount_of_sibs[i]) {
        Naval <- matrix(nrow = splits[i], ncol = 3 * (largest_sib - amount_of_sibs[i]),
                        data = -9)
        sibtable <- cbind(sibtable, Naval)
        #This should "fill" the empty spaces for if the amount of sibs < largest_sib
      }
      sibtable <- sibtable[,start := NULL]
    }



    parlg <- sweep(sweep(t(parentmatrix), 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta
    parliab <- parlg + rnorm(2*splits[i], 0, sqrt(1 - hsq))

    c_lg <- sweep(sweep(child, 2, mu, FUN = "-"), 2, sigma, FUN = "/") %*% beta


    c_liab <- c_lg + rnorm(splits[i], 0, sqrt(1 - hsq))



    #Make phenotypes:
    c_pheno <- sapply(c_liab, function(x) ifelse(x > critical, 2, 1))
    p_pheno <- sapply(parliab, function(x) ifelse(x > critical, 2, 1))


    #FID for the children/parents:

    id <- matrix(c((cusplits[i]+1):cusplits[i+1]))


    locked = flock::lock(lock) # locks file

    fwrite(as.data.table(to_ped(child, part = cusplits[i])), ## writes to locked file
           paste0(path, "genotypes.ped", sep = ""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    if (largest_sib != 0) {
      fwrite(as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab, p_pheno[seq(1, 2*splits[i], 2)], parlg[seq(1, 2*splits[i], 2)], parliab[seq(1, 2*splits[i], 2)], p_pheno[seq(2, 2*splits[i], 2)], parlg[seq(2, 2*splits[i], 2)], parliab[seq(2, 2*splits[i], 2)], sibtable)),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)}
    else {
      fwrite(as.data.table(cbind(id, rep(1, splits[i]), c_pheno, c_lg, c_liab, p_pheno[seq(1, 2*splits[i], 2)], parlg[seq(1, 2*splits[i], 2)], parliab[seq(1, 2*splits[i], 2)], p_pheno[seq(2, 2*splits[i], 2)], parlg[seq(2, 2*splits[i], 2)], parliab[seq(2, 2*splits[i], 2)])),
             paste0(path, "phenotypes.txt", sep = ""),
             quote = F,
             sep = " ",
             col.names = F,
             append = T)
    }



    flock::unlock(locked) #unlocks file

  }, future.seed = T)
}
