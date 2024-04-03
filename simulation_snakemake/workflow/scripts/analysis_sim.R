#!/usr/bin/env Rscript
library('rosace')
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

source("workflow/scripts/compare_method_utils.R")

# parse argument
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "simulation index"),
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "number of replicates"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "number of rounds"),
    make_option(c("-p", "--pos"), type = "character", default = NULL, 
        help = "whether is position-level model favored")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
idx_sim <- opt$sim
idx_rep <- opt$rep
idx_round <- opt$round
pos_folder <- opt$pos 

# data <- "CARD11"
# idx_sim <- 1
# idx_rep <- 3
# idx_round <- 1
# pos_folder <- "pos"

if (idx_round == 1 && idx_rep == 1) {
    method_list <- c("ENRICH2", "ROSACE", "ROSACENP", "SLR", "DIMSUM")
} else {
    method_list <- c("ENRICH2", "ROSACE", "ROSACENP", "SLR", "EDGER", "LIMMA", "DIMSUM")
}

if (data == "OCT1") {
    alt <- "var2"
} else {
    alt <- "var1"
}

##### load directory #####
dir <- file.path("results", data, pos_folder, 
                 paste("growth_rep", idx_rep, "_rd", idx_round, "_clean", sep = ""),
                 paste("sim", idx_sim, sep = ""))
sdir <- file.path(dir, "analysis")
effects <- read_tsv(file.path(sdir, "effects.tsv"))

##### correlation #####
df_corr <- corr_test(truth = effects$expected_effects, 
                     score = effects[str_c(method_list, "effects", sep = "_")])
write_tsv(df_corr, file = file.path(sdir, "corr.tsv"))

##### fdr #####
df_fdr_all <- data.frame()
for (fdr in c(0.005, 0.01, 0.05, 0.1)) {
  df_fdr <- comp_fdr_all(effects$expected_test,
                         (effects[str_c(method_list, "tests", sep = "_")] <= fdr))
  df_fdr$fdr <- fdr
  df_fdr_all <- rbind(df_fdr_all, df_fdr)
}
write_tsv(df_fdr_all, file = file.path(sdir, "fdr.tsv"))

##### rankfdr #####
df_rankfdr <- data.frame()
df_rankfdrsense <- data.frame()
for (m in method_list) {
  df_rankfdr <- 
    rbind(df_rankfdr,
          comp_rankfdr(effect = effects, resolution = 1000, model = m) %>% mutate(alt = FALSE),
          comp_rankfdr(effect = effects, resolution = 1000, model = m, alt = alt) %>% mutate(alt = TRUE))
  df_rankfdrsense <- 
    rbind(df_rankfdrsense, 
        comp_rankfdr_sense(effect = effects, alt = alt, model = m))
}
write_tsv(df_rankfdr, file = file.path(sdir, "rankfdr.tsv"))
write_tsv(df_rankfdrsense, file = file.path(sdir, "rankfdrsense.tsv"))

