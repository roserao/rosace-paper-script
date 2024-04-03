library("rosace")
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

##### parse argument ##### 
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein of the experiment"),
    make_option(c("-m", "--mode"), type = "character", default = "N", 
        help = "whether analyze parsed replicate data")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
mode <- opt$mode

##### directory ##### 

### load data
if (mode == "N") {
    load(file.path("data", data, "rosace.rda"))
    sdir <- file.path("results", data)
} else {
    load(file.path("replicates", data, mode, "raw", "rosace.rda"))
    sdir <- file.path("replicates", data, mode)
}
Nrep <- length(rosace@assay.sets[[1]]@rounds) 
Nround <- (rosace@assay.sets[[1]]@rounds)[1] + 1 # TODO: Nround different for each assay

### enrich2
sdir_enrich2 <- file.path(sdir, "enrich2")
if (!dir.exists(sdir_enrich2)) {
    dir.create(sdir_enrich2, recursive = TRUE)
}
sdir_enrich2_data <- file.path(sdir_enrich2, "data")
if (!dir.exists(sdir_enrich2_data)) {
    dir.create(sdir_enrich2_data, recursive = TRUE)
}
sdir_enrich2_results <- file.path(sdir_enrich2, "results")
if (!dir.exists(sdir_enrich2_results)) {
    dir.create(sdir_enrich2_results, recursive = TRUE)
}

### mutscan
sdir_mutscan <- file.path(sdir, "mutscan")
if (!dir.exists(sdir_mutscan)) {
    dir.create(sdir_mutscan, recursive = TRUE)
}

### dimsum
sdir_dimsum <- file.path(sdir, "dimsum")
if (!dir.exists(sdir_dimsum)) {
    dir.create(sdir_dimsum, recursive = TRUE)
}

##### enrich2 ##### 

df_count <- data.frame(variants = rosace@assay.sets[[1]]@var.names,
                       rosace@assay.sets[[1]]@raw.counts)
# TODO: make sure the control label is the same 
ctrl_label <- ExtractVarAssaySet(rosace, name = names(rosace@assay.sets))$type == "synonymous" 

### output counts
for (i in 1:Nrep) {
    for (j in 1:Nround) {
        if (data == "OCT1" && i == 3 && j == 4) {
            break
        }
        if (data == "OCT1" && (mode %in% c("rep2_2", "rep2_3")) && i == 2 && j == 4) {
            break
        }
        count <- df_count[(i - 1) * Nround + j + 1]
        rownames(count) <- df_count$variants
        colnames(count) <- "count"
        count[nrow(count) + 1, 1] <- sum(count[[1]][ctrl_label], na.rm = TRUE)
        rownames(count)[nrow(count)] <- "_wt"
        utils::write.table(count, quote = FALSE, sep="\t",
                           file = file.path(sdir_enrich2_data, paste("count_rep", i, "_c", j - 1, ".tsv", sep = "")))
    }
}

### output json file
experiment <- list('name' = 'simulation', 'output directory' = sdir_enrich2_results, 'conditions' = list())
condition <- list('name' = 'clean', selections = list())
for (i in 1:Nrep) {
    selection <- list('name' = paste("rep", i, sep = ""), 'libraries' = list())
    for (j in 1:Nround) {
        if (data == "OCT1" && i == 3 && j == 4) {
            break
        }
        if (data == "OCT1" && (mode %in% c("rep2_2", "rep2_3")) && i == 2 && j == 4) {
            break
        }
        seqlib = list('counts file' = file.path(sdir_enrich2_data, paste("count_rep", i, "_c", j - 1, ".tsv", sep = "")),
                    'identifiers' = list(),
                    'name' = paste("rep", i, "_c", j - 1, sep = ""),
                    'report filtered read' = FALSE,
                    'timepoint' = j - 1)
        selection$libraries[[j]] <- seqlib
    }
    condition$selections[[i]] <- selection
}
experiment$conditions[[1]] <- condition

### save json file
jsonData <- rjson::toJSON(experiment, indent = 2)
write(jsonData, file = file.path(sdir_enrich2, "config.json"))

##### mutscan ##### 
source("workflow/scripts/mutscan_utils.R")

### load count
count <- data.frame(rosace@assay.sets[[1]]@raw.counts) %>% replace(is.na(.), 0)
ctrl_label <- ExtractVarAssaySet(rosace, name = names(rosace@assay.sets))$type == "synonymous" 

### transform to SummarizedExperiment object

if (data == "OCT1" && (mode %in% c("N", "rep3_1"))) {
    se <- SummarizedExperiment::SummarizedExperiment(
        assay = list(count),
        colData = data.frame(rep = rep(1:Nrep, each = Nround)[1:11],
                             time = rep(0:(Nround-1), Nrep)[1:11]))
} else if (data == "OCT1" && (mode %in% c("rep2_2", "rep2_3"))) {
    se <- SummarizedExperiment::SummarizedExperiment(
        assay = list(count),
        colData = data.frame(rep = rep(1:Nrep, each = Nround)[1:7],
                             time = rep(0:(Nround-1), Nrep)[1:7]))
} else {
    se <- SummarizedExperiment::SummarizedExperiment(
        assay = list(count),
        colData = data.frame(rep = rep(1:Nrep, each = Nround),
                            time = rep(0:(Nround-1), Nrep)))
}

##### run edgeR and limma
if (data != "Cohesin") {
    edger_scores <- calculateRelativeFC(
        se = se,
        design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
        coef = "time",
        pseudocount = 0.5,
        WTrows = which(ctrl_label),
        method = "edgeR"
    )

    limma_scores <- calculateRelativeFC(
        se = se,
        design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
        coef = "time",
        pseudocount = 0.5,
        WTrows = which(ctrl_label),
        method = "limma"
    )
} else {
    edger_scores <- data.frame(logFC_shrunk = rep(NA, nrow(count)), FDR = rep(NA, nrow(count)))
    limma_scores <- data.frame(logFC = rep(NA, nrow(count)), adj.P.Val = rep(NA, nrow(count)))
}
write_tsv(edger_scores, file = file.path(sdir_mutscan, "score_edgeR.tsv"))
write_tsv(limma_scores, file = file.path(sdir_mutscan, "score_limma.tsv"))

##### dimsum ##### 

### load count data
name_col <- c("nt_seq")
for (i in 1:Nrep) {
name_col <- c(name_col, paste("input", i, sep = ""))
name_col <- c(name_col, paste("output", i, "A", sep = ""))
}
seq_count <- data.frame(rosace@assay.sets[[1]]@raw.counts)
if (Nround > 2) {
    seq_count <- seq_count[, sort(c(Nround * (0:(Nrep-1)) + 1, Nround * (0:(Nrep-1)) + 2))]
} 
var_data <- ExtractVarAssaySet(rosace, name = names(rosace@assay.sets))

### generate fake sequence
mut_fake <- c("ACG", "AGA", "ACA", "ATA",
"ATT", "ACT", "ATC", "ACC",
"AAC", "AGC", "ATG", "AAA",
"AAG", "AGG", "CTT", "CCT",
"CAT", "CGT", "CTC", "CCC",
"CAC", "CGC", "CTA", "CCA")
if (max(table(var_data$position)) > 23) {
    stop("not enough fake mut codon")
}
ref_fake <- rep(mut_fake[1], max(var_data$position))
mut_seq_fake <- function(pos, ref_fake, seq) {
    ref_fake[pos] <- seq
    return(str_c(ref_fake, collapse = ""))
}
mut_fake <- data.frame(mutation = unique(var_data$mutation),
                    seq = mut_fake[2:(length(unique(var_data$mutation)) + 1)])

### fake sequence in var data
var_data <- var_data %>% left_join(mut_fake)
var_data <- var_data %>% 
    rowwise() %>%
    mutate(seq_complete = mut_seq_fake(.data$position, ref_fake, .data$seq))
    ref_seq_fake <- str_c(ref_fake, collapse = "")

seq_count <- rbind(colSums(seq_count[var_data$type == "synonymous", ], na.rm = TRUE), seq_count)
seq_count <- cbind(c(ref_seq_fake, var_data$seq_complete), 
                seq_count)
colnames(seq_count) <- name_col

### save ref sequence and count
write_tsv(var_data, file = file.path(sdir_dimsum, "var_data.tsv"))
write_tsv(seq_count, file = file.path(sdir_dimsum, "dimsum_count.tsv"))
write.table(ref_seq_fake, file = file.path(sdir_dimsum, "dimsum_ref.txt"),
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)

### experiment design
exp_design <- data.frame(sample_name = name_col[-1],
                        experiment_replicate = rep(1:Nrep, each = 2),
                        selection_id = rep(0:1, Nrep),
                        selection_replicate = rep(c("", "1"), Nrep),
                        technical_replicate = 1)
write_tsv(exp_design, file = file.path(sdir_dimsum, "exp_design.tsv"))


