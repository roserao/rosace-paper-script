#!/usr/bin/env Rscript
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

##### parse argument
library('optparse')
option_list <- list(
    make_option(c("-m", "--mode"), type = "character", default = NULL, 
        help = "type of experiment: growth or binding"),
    make_option(c("-c", "--cond"), type = "character", default = NULL, 
        help = "condition: clean or reperror"),
    make_option(c("-f", "--flag"), type = "character", default = NULL, 
        help = "pos or neg. positive selection/negative selection")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
mode <- opt$mode
cond <- opt$cond
flag <- opt$flag

# mode <- "growth"
# cond <- "clean"
# flag <- "negxxx"

##### define directory
dir <-  paste('results/summary_', flag, sep = "")
dir.list <- dir(dir, pattern = paste(mode, ".*", cond, "$", sep = ""), full.names = TRUE)
sdir <- file.path(paste('results/output_', flag, sep = ""), 
                  paste(mode, cond, sep = "_"))
if (!dir.exists(sdir)) {
  dir.create(sdir)
}

df_corr <- data.frame()
df_fdr <- data.frame()
df_rankfdr <- data.frame()
df_rankfdrsense <- data.frame()
df_time <- data.frame()
for (d in dir.list) {
  
  rep <- str_extract(d, "rep[0-9]+")
  rd <- str_extract(d, "rd[0-9]+")
  
  corr_sub <- read_tsv(file.path(d, "corr.tsv")) 
  fdr_sub <- read_tsv(file.path(d, "fdr.tsv"))
  rankfdr_sub <- read_tsv(file.path(d, "rankfdr.tsv"))
  rankfdrsense_sub <- read_tsv(file.path(d, "rankfdrsense.tsv"))
  time_sub <- read_tsv(file.path(d, "benchmark.tsv"))
  
  df_corr <- rbind(df_corr, corr_sub %>% mutate(rep = rep, rd = rd))
  df_fdr <- rbind(df_fdr, fdr_sub %>% mutate(rep = rep, rd = rd))
  df_rankfdr <- rbind(df_rankfdr, rankfdr_sub %>% mutate(rep = rep, rd = rd))
  df_rankfdrsense <- rbind(df_rankfdrsense, rankfdrsense_sub %>% mutate(rep = rep, rd = rd))
  df_time <- rbind(df_time, time_sub %>% mutate(rep = rep, rd = rd))
}
rm(d, rep, rd,
   corr_sub, fdr_sub, rankfdr_sub, rankfdrsense_sub, time_sub)

write_tsv(df_corr, file.path(sdir, "corr.tsv"))
write_tsv(df_fdr, file.path(sdir, "fdr.tsv"))
write_tsv(df_rankfdr, file.path(sdir, "rankfdr.tsv"))
write_tsv(df_rankfdrsense, file.path(sdir, "rankfdrsense.tsv"))
write_tsv(df_time, file.path(sdir, "benchmark.tsv"))
