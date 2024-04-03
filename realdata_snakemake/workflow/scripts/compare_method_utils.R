# read enrich2 analysis data
extract_enrich_exp <- function(dir) {
  
  # scores
  file <- file.path(dir, "enrich2/results/tsv/simulation_exp/main_identifiers_scores.tsv")
  enrich_score <- read_tsv(file, skip = 2, col_names = FALSE)
  
  # change column name
  title <- read_tsv(file, n_max = 2, col_names = FALSE)
  colnames(enrich_score) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(enrich_score)[1] <- "var"
  rm(title)
  
  # clean the data frame
  enrich_score <- enrich_score %>% 
    pivot_longer(!var) %>%
    separate(name, into = c("cond", "stats")) %>%
    pivot_wider(names_from = stats, values_from = value)
  colnames(enrich_score)[3:5] <- str_c("enrich", colnames(enrich_score)[3:5], sep = "_")
  
  # pvalues 
  file <- file.path(dir, "enrich2/results/tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv")
  enrich_pval <- read_tsv(file, skip = 2, col_names = FALSE)
  
  # change column name
  title <- read_tsv(file, n_max = 2, col_names = FALSE)
  colnames(enrich_pval) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(enrich_pval)[1] <- "var"
  rm(title)
  
  # clean the data frame
  enrich_pval <- enrich_pval %>% 
    pivot_longer(!var) %>%
    separate(name, into = c("cond", "stats"), sep = "[.]") %>%
    pivot_wider(names_from = stats, values_from = value)
  
  colnames(enrich_pval)[3:4] <- str_c("enrich", c("pval", "z"), sep = "_")
  
  # combine score and pval
  enrich_res <- full_join(enrich_pval, enrich_score)
  
  return(enrich_res)
}
