library("rosace")
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
idx_round <- 1
cond <- opt$cond
idx_sim <- opt$sim
flag <- opt$flag

# mode <- "growth"
# idx_rep <- 1
# cond <- "clean"
# idx_sim <- 1
# flag <- "neg"

##### define directory
dir <-  file.path(paste('results/sim_', flag, sep = ""), 
                  paste(mode, "_rep", idx_rep, "_rd", idx_round, "_", cond, sep = ""),
                  paste("sim", idx_sim, sep = ""))              
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
rosace <- readRDS(file.path(dir, "rosace", "rosace.rds"))
var_data <- rosace@var.data

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





