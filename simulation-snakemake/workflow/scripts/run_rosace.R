#!/usr/bin/env Rscript
library('rosace')
library("stringr")

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
# idx_rep <- 2
# idx_round <- 2
# cond <- "clean"
# idx_sim <- 1
# flag <- "neg"

##### load data
dir <-  file.path(paste('results/sim_', flag, sep = ""), 
                  paste(mode, "_rep", idx_rep, "_rd", idx_round, "_", cond, sep = ""),
                  paste("sim", idx_sim, sep = ""))
ddir <- file.path(dir, "rosace")
rosace <- readRDS(file.path(ddir, "rosace.rds"))
key <- "simulation"

##### preprocessing
# rosace <- FilterData(rosace, key = key, na.rmax = 0.5)
rosace <- ImputeData(rosace, key = key, impute.method = "knn", na.rmax = 0.5)
rosace <- NormalizeData(rosace, 
                        key = key,
                        normalization.method = "wt", 
                        wt.var.names = rosace@var.data[[1]][str_detect(rosace@var.data[[1]], "ctrl$")], 
                        wt.rm = FALSE)
rosace <- IntegrateData(object = rosace, key = key)

##### run SLR
# rosace <- runSLR(rosace, name = key, type = "AssaySet")

##### run Rosace
rosace <- RunRosace(object = rosace,
                    name = key,
                    type = "AssaySet",
                    savedir = ddir,
                    pos.col = "position",
                    ctrl.col = "mutation",
                    ctrl.name = "ctrl",
                    install = FALSE)

##### save Rosace object
saveRDS(rosace, file = file.path(ddir, "rosace_eval.rds"))