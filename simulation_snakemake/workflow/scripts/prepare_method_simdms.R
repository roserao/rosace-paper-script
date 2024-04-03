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

idx_rep <- 5
idx_round <- 5

##### loading directory #####
dir <- file.path("results", "simdms", data, mode)
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))

##### mutscan ##### 

### load count table
count <- read_tsv(file.path(dir, "tsv", "counts.tsv"))
var_data <- rosace@var.data
count[is.na(count)] <- 0

### directory
source("workflow/scripts/mutscan_utils.R")
sdir <- file.path(dir, "mutscan")
if(!dir.exists(sdir)) {
  dir.create(sdir)
}

### transform to SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
    assay = list(count[, -1]),
    colData = data.frame(rep = rep(1:idx_rep, each = (idx_round + 1)),
                        time = rep(0:idx_round, idx_rep)))

### run edgeR and limma
edger_scores <- calculateRelativeFC(
    se = se,
    design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
    coef = "time",
    pseudocount = 0.5,
    WTrows = 1,  
    method = "edgeR"
)
limma_scores <- calculateRelativeFC(
    se = se,
    design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
    coef = "time",
    pseudocount = 0.5,
    WTrows = 1, 
    method = "limma"
)

### result
write_tsv(edger_scores, file = file.path(sdir, "score_edgeR.tsv"))
write_tsv(limma_scores, file = file.path(sdir, "score_limma.tsv"))


##### dimsum ##### 

### directory
sdir <- file.path(dir, "dimsum")
if(!dir.exists(sdir)) {
  dir.create(sdir)
}

##### load count data
name_col <- c("nt_seq")
for (i in 1:idx_rep) {
  name_col <- c(name_col, paste("input", i, sep = ""))
  name_col <- c(name_col, paste("output", i, "A", sep = ""))
}
seq_count <- read_tsv(file.path(dir, "tsv", "counts.tsv"))
seq_count <- seq_count[, -1]
seq_count[is.na(seq_count)] <- 0
if (idx_round > 2) {
    seq_count <- seq_count[, sort(c((idx_round+1) * (0:(idx_rep-1)) + 1, (idx_round+1) * (0:(idx_rep-1)) + 2))]
} 
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))
var_data <- rosace@var.data

# NEW: fake position and mutation
var_data$position <- rep(1:500, each = 20)[1:nrow(var_data)]
var_data$mutation <- rep(1:20, 500)[1:nrow(var_data)]

##### generate fake sequence
mut_fake <- c("ACG", "AGA", "ACA", "ATA",
  "ATT", "ACT", "ATC", "ACC",
  "AAC", "AGC", "ATG", "AAA",
  "AAG", "AGG", "CTT", "CCT",
  "CAT", "CGT", "CTC", "CCC",
  "CAC", "CGC", "CTA", "CCA")
ref_fake <- rep(mut_fake[1], max(var_data$position))
mut_seq_fake <- function(pos, ref_fake, seq) {
  ref_fake[pos] <- seq
  return(str_c(ref_fake, collapse = ""))
}
mut_fake <- data.frame(mutation = unique(var_data$mutation),
                       seq = mut_fake[2:(length(unique(var_data$mutation)) + 1)])

##### fake sequence in var data
var_data <- var_data %>% left_join(mut_fake)
var_data <- var_data %>% 
  rowwise() %>%
  mutate(seq_complete = mut_seq_fake(.data$position, ref_fake, .data$seq))

# ref_seq_fake <- str_c(ref_fake, collapse = "")
# seq_count <- rbind(colSums(seq_count[var_data$mutation == "ctrl", ]), seq_count)
# seq_count <- cbind(c(ref_seq_fake, var_data$seq_complete), 
#                    seq_count)
ref_seq_fake <- var_data$seq_complete[1]
seq_count <- cbind(var_data$seq_complete, seq_count)
colnames(seq_count) <- name_col
m
##### save ref sequence and count
write_tsv(var_data, file = file.path(sdir, "var_data.tsv"))
write_tsv(seq_count, file = file.path(sdir, "dimsum_count.tsv"))
write.table(ref_seq_fake, file = file.path(sdir, "dimsum_ref.txt"),
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

##### experiment design
exp_design <- data.frame(sample_name = name_col[-1],
                         experiment_replicate = rep(1:idx_rep, each = 2),
                         selection_id = rep(0:1, idx_rep),
                         selection_replicate = rep(c("", "1"), idx_rep),
                         technical_replicate = 1)
write_tsv(exp_design, file = file.path(sdir, "exp_design.tsv"))
