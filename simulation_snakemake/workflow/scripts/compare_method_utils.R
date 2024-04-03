# read enrich2 analysis data
extract_enrich_exp_simdms <- function(data, mode) {
  
  # scores
  name <- paste(data, "simulation_depth_200", sep = "_")
  dir <- file.path("results", "simdms", name)
  file <- file.path(dir, "Results/tsv", 
                    paste(name, "exp", sep = "_"),
                    "main_identifiers_scores.tsv")
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
  file <- file.path(dir, "Results/tsv", 
                    paste(name, "exp", sep = "_"),
                    "main_identifiers_scores_pvalues_wt.tsv")
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
  enrich_res <- enrich_res %>% filter(cond == mode)
  
  return(enrich_res)
}

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

# correlation test for each simulation
corr_test <- function(truth, score) {
  corr_check <- data.frame(test = colnames(score),
                           num_na = -1,
                           corr_pearson = -2,
                           corr_spearman = -2)
  for (i in 1:ncol(score)) {
    corr_check[i, 2:4] <- c(sum(is.na(score[[i]])),
                            cor(truth, score[[i]], method = "pearson", use = "complete.obs"),
                            cor(truth, score[[i]], method = "spearman", use = "complete.obs"))
  }
  return(corr_check)
}

# fdr test for a method in one simulation
comp_fdr <- function(truth, test) {
  if (length(truth) != length(test)) {
    stop("length different")
  }
  
  l <- length(truth)
  TP <- sum(truth & test, na.rm = TRUE)/l
  FN <- sum(truth & (!test), na.rm = TRUE)/l
  FP <- sum((!truth) & test, na.rm = TRUE)/l
  TN <- sum((!truth) & (!test), na.rm = TRUE)/l
  
  NAR <- sum(is.na(test))/l
  if (FP + TP > 0) {
    FDR <-  FP / (FP + TP)
  } else {
    FDR <- NA
  }
  
  if (TP + FN > 0) {
    POWER <- TP / (TP + FN)
  } else {
    POWER <- NA
  }
  
  return(c(TP, FN, FP, TN, NAR, FDR, POWER))
}

# fdr test for each simulation
comp_fdr_all <- function(truth, test_df) {
  stat_fdr <- data.frame(test = colnames(test_df),
                         TP = -1, FN = -1, FP = -1, TN = -1,
                         NAR = -1, FDR = -1, POWER = -1)
  for (i in 1:ncol(test_df)) {
    stat_fdr[i, 2:ncol(stat_fdr)] <- comp_fdr(truth, test_df[, i])
  }
  
  return(stat_fdr)
}

# rank fdr test for each simulation (no fdr threshold)
comp_rankfdr <- function(effects, resolution, model, alt) {

  effect_ranked <- effects %>% 
    dplyr::select(variants, expected_test, 
                  effect = paste(model, "effects", sep = "_"),
                  test = paste(model, "tests", sep = "_")) %>%
    drop_na(effect)
  
  if (missing(alt)) {
    effect_ranked <- effect_ranked %>%
      dplyr::mutate(abs_effect = abs(effect)) %>%
      arrange(test, desc(abs_effect))
  } else if (alt == "var1") {
    effect_ranked <- effect_ranked %>% 
        mutate(test_dir = 1 - test,
               test_dir = ifelse(effect < 0, -test_dir, test_dir)) 
    effect_ranked <- effect_ranked %>% arrange(test_dir, effect)
  } else if (alt == "var2") {
    effect_ranked <- effect_ranked %>% 
        mutate(test_dir = 1 - test,
               test_dir = ifelse(effect > 0, -test_dir, test_dir))
    effect_ranked <- effect_ranked %>% arrange(test_dir, desc(effect))
  }

  rank_fdr <- data.frame(model = model, rank = (1:resolution)/resolution,
                         TP = -1, FN = -1, FP = -1, TN = -1,
                         NAR = -1, FDR = -1, Power = -1, Sensitivity = -1)
  cutoff <- floor(nrow(effect_ranked)/resolution * (1:resolution))
  for (i in 1:resolution) {
    rank_fdr[i, 3:9] <- comp_fdr(effect_ranked$expected_test[1:cutoff[i]], 
                                 rep(TRUE, cutoff[i]))
    rank_fdr[i, 10] <- sum(effect_ranked$expected_test[1:cutoff[i]])/sum(effect_ranked$expected_test)
  }

  return(tibble(rank_fdr))
}

# comp_fdr(effect_ranked$expected_test[1:cutoff],  rep(TRUE, cutoff))
# sum(effect_ranked$expected_test[1:cutoff])/sum(effect_ranked$expected_test)

# rank fdr test for each simulation
comp_rankfdr_sense <- function(effects, alt, model) {

  effects <- effects %>% 
    dplyr::select(expected_test, 
                  starts_with(paste(model, "effects", sep = "_")),
                  starts_with(paste(model, "tests", sep = "_"))) 
  
  if (alt == "var1") {
    dir_vec <- (effects[[2]] < 0)
  } else if (alt == "var2") {
    dir_vec <- (effects[[2]] > 0)
  } else {
    stop("not supporting more than two variant groups")
  }

  fdr_seq <- c(seq(0, 0.001, by = 0.00001), seq(0.002, 0.1, by = 0.001), seq(0.11, 1, by = 0.01))
  # fdr_seq <- c(0.001, 0.01, 0.05, 0.1)

  rank_fdr <- data.frame(model = model, fdr = fdr_seq,
                         TP = -1, FN = -1, FP = -1, TN = -1,
                         NAR = -1, FDR = -1, Power = -1)

  for (i in 1:nrow(rank_fdr)) {
    fdr <- rank_fdr[i, 2]
    rank_fdr[i, 3:9] <- comp_fdr(effects$expected_test, (effects[[3]] <= fdr) & dir_vec)
  }

  return(tibble(rank_fdr))
}

# dir_vec <- (effect_ranked[[3]] < 0)
# comp_fdr(effect_ranked$expected_test, (effect_ranked[[4]] <= fdr) & dir_vec)
# max(which(effect_ranked[[4]] <= fdr & dir_vec))
# sum((effect_ranked[[4]] <= fdr) & dir_vec)