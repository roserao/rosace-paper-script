library("rosace")
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

##### parse argument ##### 
library('optparse')
option_list <- list(
    make_option(c("-k", "--key"), type = "character", default = NULL, 
        help = "key of the experiment")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
key <- opt$key

##### directory ##### 

### load data
load(file.path("data", key, "rosace.rda"))
Nrep <- length(rosace@assays) 
# TODO: Nround different for each assay
Nround <- (rosace@assay.sets[[1]]@rounds)[1] + 1

### enrich2
sdir_enrich2 <- file.path("results", key, "enrich2")
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
sdir_mutscan <- file.path("results", key, "mutscan")
if (!dir.exists(sdir_mutscan)) {
    dir.create(sdir_mutscan, recursive = TRUE)
}

### dimsum
if (Nround == 2) {
    sdir_dimsum <- file.path("results", key, "dimsum")
    if (!dir.exists(sdir_dimsum)) {
        dir.create(sdir_dimsum, recursive = TRUE)
    }
}

##### enrich2 ##### 

df_count <- data.frame(variants = rosace@assay.sets[[1]]@var.names,
                       rosace@assay.sets[[1]]@raw.counts)
# TODO: make sure the control label is the same 
ctrl_label <- ExtractVarAssaySet(rosace, name = names(rosace@assay.sets))$type == "synonymous" 

### output counts
for (i in 1:Nrep) {
    for (j in 1:Nround) {
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
se <- SummarizedExperiment::SummarizedExperiment(
  assay = list(count),
  colData = data.frame(rep = rep(1:Nrep, each = Nround),
                       time = rep(0:(Nround-1), Nrep)))
##### run edgeR and limma
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

write_tsv(edger_scores, file = file.path(sdir_mutscan, "score_edgeR.tsv"))
write_tsv(limma_scores, file = file.path(sdir_mutscan, "score_limma.tsv"))

##### dimsum ##### 
if (Nround == 2) {

    ### load count data
    name_col <- c("nt_seq")
    for (i in 1:Nrep) {
    name_col <- c(name_col, paste("input", i, sep = ""))
    name_col <- c(name_col, paste("output", i, "A", sep = ""))
    }
    seq_count <- data.frame(rosace@assay.sets[[1]]@raw.counts)
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
}


