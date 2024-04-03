#!/usr/bin/env Rscript
library('rosace')
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

# parse argument
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data


sdir <- file.path("results", data, "summary")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}

df_corr_data <- data.frame()
df_fdr_data <- data.frame()
df_rankfdr_data <- data.frame()
df_rankfdrsense_data <- data.frame()
df_time_data <- data.frame()

for (pos in c("nopos", "pos")) {
    df_corr_all <- data.frame()
    df_fdr_all <- data.frame()
    df_rankfdr_all <- data.frame()
    df_rankfdrsense_all <- data.frame()
    df_time_all <- data.frame()

    dir <- file.path("results", data, pos)
    dir.list <- dir(dir, full.names = TRUE)
    for (cdir in dir.list) { 
        rep <- str_extract(cdir, "rep[0-9]+")
        rd <- str_extract(cdir, "rd[0-9]+")

        df_corr <- data.frame()
        df_fdr <- data.frame()
        df_rankfdr <- data.frame()
        df_rankfdrsense <- data.frame()
        df_time <- data.frame()

        cdir.list <- dir(cdir, full.names = TRUE)
        for (d in cdir.list){
            sim <- str_extract(d, "sim[0-9]+")

            # correlation
            corr_sub <- read_tsv(file.path(d, "analysis", "corr.tsv")) %>%
                mutate(sim = sim)
            df_corr <- rbind(df_corr, corr_sub)
            
            # fdr
            fdr_sub <- read_tsv(file.path(d, "analysis", "fdr.tsv")) %>%
                mutate(sim = sim)
            df_fdr <- rbind(df_fdr, fdr_sub)
            
            # rankfdr
            rankfdr_sub <- read_tsv(file.path(d, "analysis", "rankfdr.tsv")) %>%
                mutate(sim = sim)
            df_rankfdr <- rbind(df_rankfdr, rankfdr_sub)

            # rankfdr sense
            rankfdrsense_sub <- read_tsv(file.path(d, "analysis", "rankfdrsense.tsv")) %>%
                mutate(sim = sim)
            df_rankfdrsense <- rbind(df_rankfdrsense, rankfdrsense_sub)
            
            # computing time
            time_sub <- read_tsv(file.path(d, "rosace", "benchmark.txt")) %>%
                mutate(sim = sim)
            df_time <- rbind(df_time, time_sub)
        }

        df_corr$rep <- rep
        df_corr$rd <- rd
        df_fdr$rep <- rep
        df_fdr$rd <- rd
        df_rankfdr$rep <- rep
        df_rankfdr$rd <- rd     
        df_rankfdrsense$rep <- rep
        df_rankfdrsense$rd <- rd
        df_time$rep <- rep
        df_time$rd <- rd

        df_corr_all <- rbind(df_corr_all, df_corr)
        df_fdr_all <- rbind(df_fdr_all, df_fdr)
        df_rankfdr_all <- rbind(df_rankfdr_all, df_rankfdr)
        df_rankfdrsense_all <- rbind(df_rankfdrsense_all, df_rankfdrsense)
        df_time_all <- rbind(df_time_all, df_time)
    }
    df_corr_all$mode <- pos
    df_fdr_all$mode <- pos
    df_rankfdr_all$mode <- pos
    df_rankfdrsense_all$mode <- pos
    df_time_all$mode <- pos

    df_corr_data <- rbind(df_corr_data, df_corr_all)
    df_fdr_data <- rbind(df_fdr_data, df_fdr_all)
    df_rankfdr_data <- rbind(df_rankfdr_data, df_rankfdr_all)
    df_rankfdrsense_data <- rbind(df_rankfdrsense_data, df_rankfdrsense_all)
    df_time_data <- rbind(df_time_data, df_time_all)
}


write_tsv(df_corr_data, file.path(sdir, "corr.tsv"))
write_tsv(df_fdr_data, file.path(sdir, "fdr.tsv"))
write_tsv(df_rankfdr_data, file.path(sdir, "rankfdr.tsv"))
write_tsv(df_rankfdrsense_data, file.path(sdir, "rankfdrsense.tsv"))
write_tsv(df_time_data, file.path(sdir, "benchmark.tsv"))




