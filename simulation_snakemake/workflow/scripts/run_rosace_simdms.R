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
        help = "binding/growth"),
    make_option(c("-m", "--mode"), type = "character", default = NULL, 
        help = "clean/reperror/jackpot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
mode <- opt$mode

##### loading directory #####
dir <- file.path("results", "simdms", data, mode)
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))

##### run rosace #####
sdir <- file.path(dir, "rosace")
rosace <- RunRosace(object = rosace,
                name = "simulation",
                type = "AssaySet",
                savedir = sdir,
                install = FALSE)
saveRDS(rosace, file = file.path(sdir, "rosace_eval.rds"))