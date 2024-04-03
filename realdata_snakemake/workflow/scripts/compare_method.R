#!/usr/bin/env Rscript
library('rosace')
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

source("workflow/scripts/compare_method_utils.R")

##### parse argument ##### 
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein of the experiment"),
    make_option(c("-m", "--mode"), type = "character", default = "N", 
        help = "whether analyze parsed replicate data")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
mode <- opt$mode

if (data == "OCT1") {
  key <- "1SM73"
} else if (data == "CARD11") {
  key <- "CARD11"
} else if (data == "MSH2") {
  key <- "MSH2"
} else if (data == "BRCA1") {
  key <- "BRCA1"
} else if (data == "MET") {
  key <- "MET1"
} else if (data == "BRCA1-RING") {
  key <- "BRCA1R"
} else if (data == "Cohesin") {
  key <- "COH"
}

if (mode == "N") {
  dir <- file.path("results", data)
} else {
  dir <- file.path("replicates", data, mode)
}

##### load no pos data #####
load(file.path(dir, "rosace_nopos", "rosace_eval.rda"))
score_nopos <- rosace@scores[[1]]
score_nopos@method <- "ROSACENP"

##### load data #####
load(file.path(dir, "rosace", "rosace_eval.rda"))

##### add no pos data #####
rosace <- AddScoreData(rosace, score_nopos)

##### ROSACE NP #####
# adjust ROSACE NP lfsr
sc <- rosace@scores[[1]]@score
sc <- sc %>% select(variants, mean, sd)
sc <- sc %>% 
  rowwise() %>%
  mutate(lfsr = min(pnorm(0, mean = mean, sd = sd),
                    pnorm(0, mean = mean, sd = sd, lower.tail = FALSE))) %>%
  mutate(lfsr = lfsr * 2) %>%
  ungroup()
rosace@scores[[1]]@score$lfsr <- sc$lfsr
rm(sc)

##### ROSACE #####
# adjust ROSACE lfsr
sc <- rosace@scores[[2]]@score
sc <- sc %>% select(variants, mean, sd)
sc <- sc %>% 
  rowwise() %>%
  mutate(lfsr = min(pnorm(0, mean = mean, sd = sd),
                    pnorm(0, mean = mean, sd = sd, lower.tail = FALSE))) %>%
  mutate(lfsr = lfsr * 2) %>%
  ungroup()
rosace@scores[[2]]@score$lfsr <- sc$lfsr
rm(sc)

##### NAIVE #####
rosace <- runSLR(rosace, name = key, type = "AssaySet")
# adjust p-value
rosace@scores[[3]]@score$`p-value` <- p.adjust(rosace@scores[[3]]@score$`p-value`, method = "fdr")

##### ENRICH2 #####
df_enrich2 <- extract_enrich_exp(dir)
df_enrich2 <- df_enrich2 %>% 
    dplyr::select(variants = var, estimate = enrich_score, se = enrich_SE, pvalue = enrich_pval) %>%
    mutate(pvalue = p.adjust(pvalue, method = "fdr"))
score_enrich2 <- CreateScoreObject(method = "ENRICH2",
                                   type = "growth",
                                   assay.name = key,
                                   score = df_enrich2)
rosace <- AddScoreData(rosace, score_enrich2)
rm(score_enrich2, df_enrich2)

##### load mutscan score except for (1, 1) scenario
edger_scores <- read_tsv(file.path(dir, "mutscan", "score_edgeR.tsv"))
limma_scores <- read_tsv(file.path(dir, "mutscan", "score_limma.tsv"))
df_edger <- data.frame(variants = rosace@assay.sets[[1]]@var.names)
df_limma <- data.frame(variants = rosace@assay.sets[[1]]@var.names)
df_edger <- cbind(df_edger, edger_scores %>% mutate(se = NA) %>% select(score = logFC_shrunk, se, fdr = FDR))
if (data == "Cohesin") {
  df_limma <- cbind(df_limma, limma_scores %>% mutate(se.logFC = NA) %>% select(score = logFC, se = se.logFC, fdr = adj.P.Val))
} else {    
  df_limma <- cbind(df_limma, limma_scores %>% select(score = logFC, se = se.logFC, fdr = adj.P.Val))
}
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

##### load dimsum score if round == 1
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
                            pnorm(q = abs(z_score), lower.tail = TRUE))) %>%
    mutate(pvalue = p.adjust(pvalue, method = "fdr"))
score_dimsum <- CreateScoreObject(method = "DIMSUM", 
                                    type = "growth",
                                    assay.name = key,
                                    score = df_dimsum %>% dplyr::select(variants = variant, estimate, sd, pvalue))
rosace <- AddScoreData(rosace, score_dimsum)

 
##### analysis save directory ##### 
sdir <- file.path(dir, "analysis")
if (!dir.exists(sdir)){
  dir.create(sdir)
}
save(rosace, file = file.path(sdir, "rosace_full.rda"))

##### create effects table ##### 
effects <- data.frame(rosace@var.data)
for (score in rosace@scores) {
    df_method <- score@score[, c(1, 2, 4)] %>% filter(variants != "_wt")
    colnames(df_method) <- c("variants", 
                        paste(score@method, "effects", sep = "_"),
                        paste(score@method, "tests", sep = "_"))
    effects <- left_join(effects, df_method)
}
effects <- effects %>% drop_na(ROSACE_effects, ROSACE_tests)
write_tsv(effects, file = file.path(sdir, "effects.tsv"))
