#' Title
#'
#' @param folder_path
#' @param assoc 
#' @param true_effects 
#' @param m 
#'
#' @importFrom magrittr %>% 
#' @importFrom tibble tibble
#' @return 
#' @export
#'
#' @examples
merge_results <- function(folder_path, assoc, true_effects = "beta.txt", m=10000) {
  results <- tibble(read.table(file=(paste(folder_path, assoc, sep="/")),
                               header = TRUE))
  snp_effects <- tibble(read.table(file=(paste(folder_path, true_effects, sep="/")),
                                   col.names = "true_effect")) %>% 
    rowid_to_column("SNP") %>% 
    {.}
  
  merged_tbl <- merge(results, snp_effects) %>% 
    mutate(causal = true_effect != 0,
           signif = P < 0.05,
           bonferroni = P < 0.05/m) %>%
    {.}
  
  return(merged_tbl)
}