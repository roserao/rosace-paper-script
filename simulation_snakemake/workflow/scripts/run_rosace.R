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
        help = "target protein"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "simulation index"),
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "number of replicates"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "number of rounds"),
    make_option(c("-p", "--pos"), type = "character", default = NULL, 
        help = "whether is position-level model favored"),
    make_option(c("-m", "--mode"), type = "character", default = "P", 
        help = "position grouping or not")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
idx_sim <- opt$sim
idx_rep <- opt$rep
idx_round <- opt$round
pos_folder <- opt$pos 
mode <- opt$mode

##### loading directory #####
dir <- file.path("results", data, pos_folder, 
                 paste("growth_rep", idx_rep, "_rd", idx_round, "_clean", sep = ""),
                 paste("sim", idx_sim, sep = ""))
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))

##### preprocessing
key <- "simulation"
rosace <- FilterData(rosace, key = key, na.rmax = 0.5, min.count = 20)
rosace <- ImputeData(rosace, key = key, impute.method = "zero")
rosace <- NormalizeData(rosace, 
                        key = key,
                        normalization.method = "wt", 
                        wt.var.names = rosace@var.data[[1]][str_detect(rosace@var.data[[1]], "ctrl$")], 
                        wt.rm = FALSE)
rosace <- IntegrateData(object = rosace, key = key)


##### run Rosace
if (mode == "P") {
    sdir <- file.path(dir, "rosace")
    rosace <- RunRosace(object = rosace,
                        name = key,
                        type = "AssaySet",
                        savedir = sdir,
                        pos.col = "position",
                        ctrl.col = "mutation",
                        ctrl.name = "ctrl",
                        install = FALSE)
    saveRDS(rosace, file = file.path(sdir, "rosace_eval.rds"))
} else if (mode == "N") {
    sdir <- file.path(dir, "rosace_nopos")
    if (!dir.exists(sdir)) {
        dir.create(sdir, recursive = TRUE)
    } 
    rosace <- RunRosace(object = rosace,
                    name = key,
                    type = "AssaySet",
                    savedir = sdir,
                    # pos.col = "position",
                    # ctrl.col = "mutation",
                    # ctrl.name = "ctrl",
                    install = FALSE)
    saveRDS(rosace, file = file.path(sdir, "rosace_eval.rds"))
}
