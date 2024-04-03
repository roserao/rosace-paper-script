#!/usr/bin/env Rscript
library('rosace')
library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("ggpubr") 

source("workflow/scripts/compare_method_utils.R")

# parse argument
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "binding/growth"),
    make_option(c("-m", "--mode"), type = "character", default = NULL, 
        help = "clean/reperror/jackpot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
mode <- opt$mode

# data <- "growth"
# mode <- "clean"

method_list <- c("ENRICH2", "ROSACENP", "SLR", "EDGER", "LIMMA", "DIMSUM")
model_list <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace(nopos)", "Naive")
##### loading directory #####
dir <- file.path("results", "simdms", data, mode)
sdir <- file.path(dir, "plot")
if (!dir.exists(sdir)){
  dir.create(sdir)
}

effects <- read_tsv(file.path(dir, "analysis", "effects.tsv"))

##### correlation plot #####

df_plot <- effects %>% 
    dplyr::select(variants, truth, ends_with("effects")) %>%
    pivot_longer(cols = ends_with("effects"), names_to = "method", values_to = "effect")
df_plot$method <- as.factor(df_plot$method)
levels(df_plot$method) <- model_list

p <- ggplot(df_plot, aes(truth, effect)) +
    geom_point(size = 0.1) +
    facet_wrap(vars(method)) +
    labs(title = paste(data, mode, sep = " ")) +
    theme_cowplot() +
    stat_cor(method = "pearson", label.x = -2.5, label.y = 2.5, digits = 5) +
    geom_smooth(method = "lm", se = FALSE, col = "deeppink") 
save_plot(p, file = file.path(sdir, "corr.png"), base_width = 14, base_height = 8)

# df_corr <- corr_test(truth = effects$truth, 
#                      score = effects[str_c(method_list, "effects", sep = "_")])