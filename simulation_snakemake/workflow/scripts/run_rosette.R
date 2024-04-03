#!/usr/bin/env Rscript
library('rosace')

# parse argument
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "number of simulations"),
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
n_sim <- opt$sim
n_rep <- opt$rep
n_round <- opt$round
pos_folder <- opt$pos 

if (pos_folder == "pos") {
    pos <- TRUE
} else if (pos_folder == "nopos") {
    pos <- FALSE
}

##### load resette object #####
load(file.path("results", "rosette", data, "rosette.rda"))
sdir <- file.path("results", data, pos_folder)
if(!dir.exists(sdir)) {
  dir.create(sdir)
}

##### create config #####
if (data == "OCT1") {
    cfg <- CreateConfig(rosette,
                        n.sim = n_sim, save.sim = sdir, type.sim = "growth",
                        n.rep = n_rep, n.round = n_round,
                        null.var.group = 'var1', wt.effect = -2,
                        seq.shrink = 1.2, seq.depth = 100,
                        lib.shrink = 1,
                        var.shrink = 1, pos.flag = pos)
} else if (data == "MET") {
    cfg <- CreateConfig(rosette,
                        n.sim = n_sim, save.sim = sdir, type.sim = "growth",
                        n.rep = n_rep, n.round = n_round,
                        null.var.group = 'var2', wt.effect = 2,
                        seq.shrink = 1.2, seq.depth = 100,
                        lib.shrink = 1,
                        var.shrink = 1, pos.flag = pos)
} else if (data == "CARD11") {
    cfg <- CreateConfig(rosette,
                        n.sim = n_sim, save.sim = sdir, type.sim = "growth",
                        n.rep = n_rep, n.round = n_round,
                        null.var.group = 'var2', wt.effect = 2,
                        seq.shrink = 1.2, seq.depth = 100,
                        lib.shrink = 1,
                        var.shrink = 1, pos.flag = pos)
}
runRosette(config = cfg, save.tsv = TRUE, save.rosace = TRUE, save.enrich2 = TRUE)


