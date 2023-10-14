#!/usr/bin/env Rscript
library('rosace')
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

source("workflow/scripts/compare_method_utils.R")

##### parse argument
library('optparse')
option_list <- list(
    make_option(c("-m", "--mode"), type = "character", default = NULL, 
        help = "type of experiment: growth or binding"),
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "replicate index"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "round index"),
    make_option(c("-c", "--cond"), type = "character", default = NULL, 
        help = "condition: clean or reperror"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "simulation index"),
    make_option(c("-f", "--flag"), type = "character", default = NULL, 
        help = "pos or neg. positive selection/negative selection")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
mode <- opt$mode
idx_rep <- opt$rep
idx_round <- opt$round
cond <- opt$cond
idx_sim <- opt$sim
flag <- opt$flag

# mode <- "growth"
# idx_rep <- 3
# idx_round <- 3
# cond <- "clean"
# idx_sim <- 1
# flag <- "neg"

if (flag == "pos" || flag == "posxxx") {
  alt <- "var2"
} else if (flag == "neg" || flag == "negxxx") {
  alt <- "var1"
} else {
  stop("flag can only be neg or pos.")
}

##### load data
dir <-  file.path(paste('results/sim_', flag, sep = ""), 
                  paste(mode, "_rep", idx_rep, "_rd", idx_round, "_", cond, sep = ""),
                  paste("sim", idx_sim, sep = ""))
rosace <- readRDS(file.path(dir, "rosace", "rosace_eval.rds"))
key <- "simulation"

# adjust Rosace LFSR
sc <- cbind(rosace@scores[[1]]@score, rosace@scores[[1]]@optional.score)
sc <- sc %>% select(variant, mean, sd = `stats::sd`)
sc <- sc %>% 
  rowwise() %>%
  mutate(lfsr = min(pnorm(0, mean = mean, sd = sd),
                    pnorm(0, mean = mean, sd = sd, lower.tail = FALSE))) %>%
  ungroup()
rosace@scores[[1]]@score$lfsr <- sc$lfsr
rm(sc)

##### run SLR 
rosace <- runSLR(rosace, name = key, type = "AssaySet")

##### load enrich2 score
df_enrich2 <- extract_enrich_exp(dir)
df_enrich2 <- df_enrich2 %>% dplyr::select(variant = var, estimate = enrich_score, pvalue = enrich_pval)
score_enrich2 <- CreateScoreObject(method = "ENRICH2",
                                   type = "growth",
                                   assay.name = key,
                                   score = df_enrich2)
rosace <- AddScoreData(rosace, score_enrich2)
rm(score_enrich2, df_enrich2)

##### load mutscan score except for (1, 1) scenario
if (!(idx_round == 1 && idx_rep == 1)) {
  edger_scores <- read_tsv(file.path(dir, "mutscan", "score_edgeR.tsv"))
  limma_scores <- read_tsv(file.path(dir, "mutscan", "score_limma.tsv"))

  df_edger <- data.frame(variant = rosace@var.data$variants)
  df_limma <- data.frame(variant = rosace@var.data$variants)
  df_edger <- cbind(df_edger, edger_scores %>% select(score = logFC_shrunk, fdr = FDR))
  df_limma <- cbind(df_limma, limma_scores %>% select(score = logFC, fdr = adj.P.Val))

  score_edger <- CreateScoreObject(method = "EDGER",
                                  type = "growth",
                                  assay.name = key,
                                  score = df_edger)
  score_limma <- CreateScoreObject(method = "LIMMA",
                                  type = "growth",
                                  assay.name = key,
                                  score = df_limma)
  rosace <- AddScoreData(rosace, score_edger)
  rosace <- AddScoreData(rosace, score_limma)
}

##### load dimsum score if round == 1
if (idx_round == 1) {
  load(file.path(dir, "dimsum", "DiMSum_Project", "DiMSum_Project_fitness_replicates.RData"))
  var_data <- read_tsv(file.path(dir, "dimsum", "var_data.tsv")) %>% 
    select(variant = variants, nt_seq = seq_complete) %>%
    mutate(nt_seq = tolower(nt_seq))
  all_variants <- all_variants %>% left_join(var_data) %>% filter(is.na(WT))
  df_dimsum <- data.frame(variant = all_variants$variant,
                          estimate = all_variants$fitness,
                          sd = all_variants$sigma) 
  df_dimsum <- df_dimsum %>% 
    rowwise() %>%
    mutate(z_score = estimate / sd,
           pvalue = 2 * min(pnorm(q = abs(z_score), lower.tail = FALSE),
                            pnorm(q = abs(z_score), lower.tail = TRUE)))
  score_dimsum <- CreateScoreObject(method = "DIMSUM", 
                                    type = "growth",
                                    assay.name = key,
                                    score = df_dimsum %>% 
                                      dplyr::select(variant, estimate, pvalue),
                                    optional.score = df_dimsum %>% 
                                      dplyr::select(z_score, sd))
  rosace <- AddScoreData(rosace, score_dimsum)
}

##### analysis save directory
sdir <- file.path(dir, "analysis")
if (!dir.exists(sdir)){
  dir.create(sdir)
}
save(rosace, file = file.path(sdir, "rosace_full.RData"))

##### create effects table
# load ground truth
effects <- rosace@misc$effects[, 1:3] # variants expected_effects
effects <- effects %>% mutate(expected_test = ifelse(expected_group == alt, TRUE, FALSE))
# load method results
for (score in rosace@scores) {
  df_method <- score@score[, 1:3] %>% filter(variant != "_wt")
  colnames(df_method) <- c("variants", 
                           paste(score@method, "effects", sep = "_"),
                           paste(score@method, "tests", sep = "_"))
  effects <- left_join(effects, df_method)
}
rm(df_method, score)
# p-value adjustment for frequentist method
effects <- effects %>% 
  mutate(SLR_tests = p.adjust(SLR_tests, method = "fdr"),
         ENRICH2_tests = p.adjust(ENRICH2_tests, method = "fdr"),
         ROSACE_tests = 2 * ROSACE_tests)
if (idx_round == 1) {
  effects <- effects %>%
    mutate(DIMSUM_tests = p.adjust(DIMSUM_tests, method = "fdr"))
}
write_tsv(effects, file = file.path(sdir, "effects_full.tsv"))

##### method list
method_list <- c("ROSACE", "SLR", "ENRICH2")
if (idx_round == 1) {
  method_list <- c(method_list, "DIMSUM")
}
if (!(idx_round == 1 && idx_rep == 1)) {
  method_list <- c(method_list, "EDGER", "LIMMA")
}

##### correlation
df_corr <- corr_test(truth = effects$expected_effects, 
                     score = effects[str_c(method_list, "effects", sep = "_")])
write_tsv(df_corr, file = file.path(sdir, "corr.tsv"))

##### fdr
df_fdr_all <- data.frame()
for (fdr in c(0, 0.005, 0.01, 0.05)) {
  df_fdr <- comp_fdr_all(effects$expected_test,
                         (effects[str_c(method_list, "tests", sep = "_")] <= fdr))
  df_fdr$fdr <- fdr
  df_fdr_all <- rbind(df_fdr_all, df_fdr)
}
write_tsv(df_fdr_all, file = file.path(sdir, "fdr.tsv"))

##### rankfdr
df_rankfdr <- data.frame()
df_rankfdrsense <- data.frame()
for (m in method_list) {
  df_rankfdr <- 
    rbind(df_rankfdr,
          comp_rankfdr(effect = effects, alt = alt, resolution = 1000, model = m))
  df_rankfdrsense <- 
    rbind(df_rankfdrsense,
          comp_rankfdr_sense(effect = effects, alt = alt, model = m))
}
write_tsv(df_rankfdr, file = file.path(sdir, "rankfdr.tsv"))
write_tsv(df_rankfdrsense, file = file.path(sdir, "rankfdrsense.tsv"))

##### plot
library(ggplot2)
library(cowplot)
# rankfdr plot
df_rankfdr_plot <- df_rankfdr %>% 
  select(model, rank, FDR, Sensitivity) %>%
  pivot_longer(cols = c(FDR, Sensitivity), names_to = "stats", values_to = "value")
p <- ggplot(df_rankfdr_plot, aes(x = rank, y = value)) +
  geom_step(aes(color = model), linewidth = 1, alpha = 0.5) +
  facet_wrap(~stats) +
  geom_segment(aes(x = 0, y = 0, 
                   xend = sum(effects$expected_test)/nrow(effects), yend = 1), 
               color = "black", linetype = "dashed")
save_plot(file = file.path(sdir, "rankfdr.png"), p, base_width = 10, base_height = 4)
rm(df_rankfdr_plot)
# ranksense plot
fdrsense_sig <- df_rankfdrsense %>% filter(fdr %in% c(0, 0.01, 0.05, 0.1))
p <- ggplot(df_rankfdrsense, aes(x = FDR, y = Power)) +
  geom_path(aes(color = model), alpha = 0.5) +
  labs(x = "false discovery rate", y = "power", color = NULL, shape = NULL) +
  geom_point(aes(x = FDR, y = Power, color = model, shape = as.factor(fdr)), 
             data = fdrsense_sig, size = 2, alpha = 0.8) 
save_plot(file = file.path(sdir, "rankfdrsense.png"), p, base_width = 6, base_height = 4)
