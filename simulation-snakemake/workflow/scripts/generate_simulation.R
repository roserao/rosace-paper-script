#!/usr/bin/env Rscript
library('rosace')

# parse argument
library('optparse')
option_list <- list(
    make_option(c("-m", "--mode"), type = "character", default = NULL, 
        help = "type of experiment: growth or binding"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "number of simulations"),
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "number of replicates"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "number of rounds"),
    make_option(c("-f", "--flag"), type = "character", default = NULL, 
        help = "pos or neg. positive selection/negative selection"),
    make_option(c("-p", "--pos"), type = "logical", default = NULL, 
        help = "whether is position-level model favored")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
mode <- opt$mode
n_sim <- opt$sim
n_rep <- opt$rep
n_round <- opt$round
flag <- opt$flag
pos <- opt$pos

# start simulation
if (flag == "pos") {
    load("data/rosette_1SM73.RData")
    if (pos) {
        cfg.clean <- CreateConfig(rosette,
                                n.sim = n_sim, save.sim = "results/sim_pos/", type.sim = mode,
                                n.rep = n_rep, n.round = n_round, 
                                null.var.group = 'var1', wt.effect = -2,
                                seq.shrink = 1.2, seq.depth = 100,
                                lib.shrink = 1,
                                var.shrink = 1, pos.flag = TRUE,
                                mode.sim = "clean")
    } else {
        cfg.clean <- CreateConfig(rosette,
                                n.sim = n_sim, save.sim = "results/sim_posxxx/", type.sim = mode,
                                n.rep = n_rep, n.round = n_round, 
                                null.var.group = 'var1', wt.effect = -2,
                                seq.shrink = 1.2, seq.depth = 100,
                                lib.shrink = 1,
                                var.shrink = 1, pos.flag = FALSE,
                                mode.sim = "clean")
    }

} else if (flag == "neg") {
    load("data/rosette_MET.RData")
    if (pos) {
        cfg.clean <- CreateConfig(rosette,
                                n.sim = n_sim, save.sim = "results/sim_neg/", type.sim = mode,
                                n.rep = n_rep, n.round = n_round, 
                                null.var.group = 'var2', wt.effect = 2,
                                seq.shrink = 1.2, seq.depth = 100,
                                lib.shrink = 1,
                                var.shrink = 1, pos.flag = TRUE,
                                mode.sim = "clean")
    } else {
        cfg.clean <- CreateConfig(rosette,
                                n.sim = n_sim, save.sim = "results/sim_negxxx/", type.sim = mode,
                                n.rep = n_rep, n.round = n_round, 
                                null.var.group = 'var2', wt.effect = 2,
                                seq.shrink = 1.2, seq.depth = 100,
                                lib.shrink = 1,
                                var.shrink = 1, pos.flag = FALSE,
                                mode.sim = "clean")
    }
} else {
    stop("Binding simulation not implemented yet.")
}
runRosette(config = cfg.clean, save.tsv = TRUE, save.rosace = TRUE, save.enrich2 = TRUE)


