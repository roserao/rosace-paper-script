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
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "replicate index"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "round index"),
    make_option(c("-c", "--cond"), type = "character", default = NULL, 
        help = "condition: clean or reperror"),
    make_option(c("-f", "--flag"), type = "character", default = NULL, 
        help = "pos or neg. positive selection/negative selection")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
mode <- opt$mode
idx_rep <- opt$rep
idx_round <- opt$round
cond <- opt$cond
flag <- opt$flag

# mode <- "growth"
# idx_rep <- 1
# idx_round <- 1
# cond <- "clean"
# flag <- "negxxx"

##### define directory
dir <-  file.path(paste('results/sim_', flag, sep = ""), 
                  paste(mode, "_rep", idx_rep, "_rd", idx_round, "_", cond, sep = ""))
dir.list <- dir(dir, full.names = TRUE)

sdir <- file.path(paste('results/summary_', flag, sep = ""), 
                  paste(mode, "_rep", idx_rep, "_rd", idx_round, "_", cond, sep = ""))
if (!dir.exists(sdir)) {
  dir.create(sdir)
}

##### load data
df_corr <- data.frame()
df_fdr <- data.frame()
df_rankfdr <- data.frame()
df_rankfdrsense <- data.frame()
df_time <- data.frame()
for (d in dir.list) {
  
  # correlation
  corr_sub <- read_tsv(file.path(d, "analysis", "corr.tsv")) %>%
    mutate(sim = sub(".*/", "", d))
  df_corr <- rbind(df_corr, corr_sub)
  
  # fdr
  fdr_sub <- read_tsv(file.path(d, "analysis", "fdr.tsv")) %>%
    mutate(sim = sub(".*/", "", d))
  df_fdr <- rbind(df_fdr, fdr_sub)
  
  # rankfdr
  rankfdr_sub <- read_tsv(file.path(d, "analysis", "rankfdr.tsv")) %>%
    mutate(sim = sub(".*/", "", d))
  df_rankfdr <- rbind(df_rankfdr, rankfdr_sub)

  # rankfdr sense
  rankfdrsense_sub <- read_tsv(file.path(d, "analysis", "rankfdrsense.tsv")) %>%
    mutate(sim = sub(".*/", "", d))
  df_rankfdrsense <- rbind(df_rankfdrsense, rankfdrsense_sub)
  
  # computing time
  time_sub <- read_tsv(file.path(d, "rosace", "benchmark.txt")) %>%
    mutate(sim = sub(".*/", "", d))
  df_time <- rbind(df_time, time_sub)
  
}
rm(corr_sub, fdr_sub, rankfdr_sub, rankfdrsense_sub, time_sub, d)

write_tsv(df_corr, file.path(sdir, "corr.tsv"))
write_tsv(df_fdr, file.path(sdir, "fdr.tsv"))
write_tsv(df_rankfdr, file.path(sdir, "rankfdr.tsv"))
write_tsv(df_rankfdrsense, file.path(sdir, "rankfdrsense.tsv"))
write_tsv(df_time, file.path(sdir, "benchmark.tsv"))



