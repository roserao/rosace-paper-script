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
        help = "whether is position-level model favored")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
idx_sim <- opt$sim
idx_rep <- opt$rep
idx_round <- opt$round
pos_folder <- opt$pos 

# data <- "CARD11"
# idx_sim <- 1
# idx_rep <- 3
# idx_round <- 3
# pos_folder <- "pos"

##### loading directory #####
dir <- file.path("results", data, pos_folder, 
                 paste("growth_rep", idx_rep, "_rd", idx_round, "_clean", sep = ""),
                 paste("sim", idx_sim, sep = ""))
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))

##### mutscan ##### 

### load count table
count <- read_tsv(file.path(dir, "tsv", "sequencing_counts.tsv"))
var_data <- rosace@var.data

### directory
source("workflow/scripts/mutscan_utils.R")
sdir <- file.path(dir, "mutscan")
if(!dir.exists(sdir)) {
  dir.create(sdir)
}

### transform to SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
  assay = list(count),
  colData = data.frame(rep = rep(1:idx_rep, each = (idx_round + 1)),
                       time = rep(0:idx_round, idx_rep)))

### run edgeR and limma
if (idx_rep == 1 && idx_round == 1) {
  edger_scores <- data.frame(logFC_shrunk = rep(NA, nrow(count)), FDR = rep(NA, nrow(count)))
  limma_scores <- data.frame(logFC = rep(NA, nrow(count)), adj.P.Val = rep(NA, nrow(count)))
} else {
  edger_scores <- calculateRelativeFC(
    se = se,
    design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
    coef = "time",
    pseudocount = 0.5,
    WTrows = which(var_data$mutation == "ctrl"),  
    method = "edgeR"
  )

  limma_scores <- calculateRelativeFC(
    se = se,
    design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
    coef = "time",
    pseudocount = 0.5,
    WTrows = which(var_data$mutation == "ctrl"), #### change here!! hardcode 22
    method = "limma"
  )
}

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
seq_count <- read_tsv(file.path(dir, "tsv", "sequencing_counts.tsv"))
if (idx_round > 2) {
    seq_count <- seq_count[, sort(c((idx_round+1) * (0:(idx_rep-1)) + 1, (idx_round+1) * (0:(idx_rep-1)) + 2))]
} 
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))
var_data <- rosace@var.data

# var_data <- read_tsv(file.path(dir, "tsv", "effects.tsv")) %>% select(variants)
#   rosace@var.data <- rosace@var.data %>%
#     tidyr::separate(.data$variants, into = c("position", "mutation"), remove = FALSE) %>%
#     mutate(position = as.numeric(substr(.data$position, 4, nchar(.data$position))))

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
ref_seq_fake <- str_c(ref_fake, collapse = "")

seq_count <- rbind(colSums(seq_count[var_data$mutation == "ctrl", ]), seq_count)
seq_count <- cbind(c(ref_seq_fake, var_data$seq_complete), 
                   seq_count)
colnames(seq_count) <- name_col

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










