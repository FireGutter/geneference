#' Title
#'
#' @param folder_path
#' @param assoc 
#' @param true_effects 
#' @param m 
#'
#' @importFrom magrittr %>% 
#' @import tibble
#' @import dplyr
#' @importFrom data.table fread
#' @return something
#' @export
#'

merge_results <- function(folder_path, assoc, true_effects = "beta.txt", m) {
  results <- tibble::tibble(
    data.table::fread(file = (paste(folder_path, assoc, sep = "/")),
    header = TRUE))
  snp_effects <- tibble::tibble(
    data.table::fread(file = (paste(folder_path, true_effects, sep = "/")),
               col.names = "true_effect")) %>% 
    tibble::rowid_to_column("SNP") %>% 
    {.}
  
  merged_tbl <- base::merge(results, snp_effects) %>% 
    dplyr::mutate(causal = true_effect != 0,
                  signif = P < 0.05,
                  bonferroni = P < 0.05/m) %>%
    {.}
  
  return(merged_tbl)
}