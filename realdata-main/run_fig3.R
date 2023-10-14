library(rosace)


dir <- "data/OCT1/rosace"

########## test: 123
load("data/OCT1/rosace_1SM73.RData")
sdir <- file.path(dir, "test123")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))

########## test: 1
sdir <- file.path(dir, "test1")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
load("data/OCT1/rosace_1SM73.RData")
rosace <- CreateRosaceObject(object = ExtractAssay(rosace, name = "1SM73_1"),
                   var.data = rosace@var.data)
rosace <- IntegrateData(rosace, key = "1SM73")
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))

########## test: 2
sdir <- file.path(dir, "test2")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
load("data/OCT1/rosace_1SM73.RData")
rosace <- CreateRosaceObject(object = ExtractAssay(rosace, name = "1SM73_2"),
                             var.data = rosace@var.data)
rosace <- IntegrateData(rosace, key = "1SM73")
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))

########## test: 3
sdir <- file.path(dir, "test3")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
load("data/OCT1/rosace_1SM73.RData")
rosace <- CreateRosaceObject(object = ExtractAssay(rosace, name = "1SM73_3"),
                             var.data = rosace@var.data)
rosace <- IntegrateData(rosace, key = "1SM73")
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))

########## test: 12
sdir <- file.path(dir, "test12")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
load("data/OCT1/rosace_1SM73.RData")
assay_x <- ExtractAssay(rosace, name = "1SM73_2")
rosace <- CreateRosaceObject(object = ExtractAssay(rosace, name = "1SM73_1"),
                             var.data = rosace@var.data)
rosace <- AddAssayData(object = rosace, assay = assay_x,
                       var.data = rosace@var.data)
rosace <- IntegrateData(rosace, key = "1SM73")
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))


########## test: 23
sdir <- file.path(dir, "test23")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
load("data/OCT1/rosace_1SM73.RData")
assay_x <- ExtractAssay(rosace, name = "1SM73_3")
rosace <- CreateRosaceObject(object = ExtractAssay(rosace, name = "1SM73_2"),
                             var.data = rosace@var.data)
rosace <- AddAssayData(object = rosace, assay = assay_x,
                       var.data = rosace@var.data)
rosace <- IntegrateData(rosace, key = "1SM73")
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))

########## test: 31
sdir <- file.path(dir, "test31")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}
load("data/OCT1/rosace_1SM73.RData")
assay_x <- ExtractAssay(rosace, name = "1SM73_1")
rosace <- CreateRosaceObject(object = ExtractAssay(rosace, name = "1SM73_3"),
                             var.data = rosace@var.data)
rosace <- AddAssayData(object = rosace, assay = assay_x,
                       var.data = rosace@var.data)
rosace <- IntegrateData(rosace, key = "1SM73")
rosace <- runSLR(rosace, name = "1SM73", type = "AssaySet")
save(rosace, file = file.path(sdir, "rosace.RData"))
# RunRosace(rosace,                      
#           name = "1SM73",
#           type = "AssaySet",
#           savedir = sdir,
#           pos.col = "position",
#           install = FALSE)
# save(rosace, file = file.path(sdir, "rosace_eval.RData"))



###
library(tidyverse)
dir <- "data/OCT1/enrich2/count"
count1 <- read_tsv("data/OCT1/count/1SM73/1SM73R1_sel/main_identifiers_counts_unfiltered.tsv")
rep1_c0 <- count1[, c(1,2)]
rep1_c1 <- count1[, c(1,3)]
rep1_c2 <- count1[, c(1,4)]
rep1_c3 <- count1[, c(1,5)]
colnames(rep1_c0) <- c("hgvs", "count")
colnames(rep1_c1) <- c("hgvs", "count")
colnames(rep1_c2) <- c("hgvs", "count")
colnames(rep1_c3) <- c("hgvs", "count")
write_tsv(rep1_c0, file = file.path(dir, "rep1_c_0.tsv"))
write_tsv(rep1_c1, file = file.path(dir, "rep1_c_1.tsv"))
write_tsv(rep1_c2, file = file.path(dir, "rep1_c_2.tsv"))
write_tsv(rep1_c3, file = file.path(dir, "rep1_c_3.tsv"))

count2 <- read_tsv("data/OCT1/count/1SM73/1SM73R2_sel/main_identifiers_counts_unfiltered.tsv")
rep2_c0 <- count2[, c(1,2)]
rep2_c1 <- count2[, c(1,3)]
rep2_c2 <- count2[, c(1,4)]
rep2_c3 <- count2[, c(1,5)]
colnames(rep2_c0) <- c("hgvs", "count")
colnames(rep2_c1) <- c("hgvs", "count")
colnames(rep2_c2) <- c("hgvs", "count")
colnames(rep2_c3) <- c("hgvs", "count")
write_tsv(rep2_c0, file = file.path(dir, "rep2_c_0.tsv"))
write_tsv(rep2_c1, file = file.path(dir, "rep2_c_1.tsv"))
write_tsv(rep2_c2, file = file.path(dir, "rep2_c_2.tsv"))
write_tsv(rep2_c3, file = file.path(dir, "rep2_c_3.tsv"))

count3 <- read_tsv("data/OCT1/count/1SM73/1SM73R3_sel/main_identifiers_counts_unfiltered.tsv")
rep3_c0 <- count3[, c(1,2)]
rep3_c1 <- count3[, c(1,3)]
rep3_c2 <- count3[, c(1,4)]
colnames(rep3_c0) <- c("hgvs", "count")
colnames(rep3_c1) <- c("hgvs", "count")
colnames(rep3_c2) <- c("hgvs", "count")
write_tsv(rep3_c0, file = file.path(dir, "rep3_c_0.tsv"))
write_tsv(rep3_c1, file = file.path(dir, "rep3_c_1.tsv"))
write_tsv(rep3_c2, file = file.path(dir, "rep3_c_2.tsv"))


# enrich_cmd data/OCT1/enrich2/config-test123.json WLS wt --no-plots

##### mutscan
# Add Score Object to the Rosace Object
AddScoreData <- function(object, score) {
  
  if (!isa(object, "Rosace")) {
    stop("object argument is not a 'Rosace' object")
  }
  if (!isa(score, "Score")) {
    stop("score argument is not a 'Score' object")
  }
  
  score.list <- list(score)
  names(score.list) <- names(score)
  
  if (names(score.list) %in% names(object@scores)) {
    warning(
      "Score object already exists. Update.",
      call. = FALSE,
      immediate. = TRUE
    )
    object@scores[[names(score.list)]] <- score
  } else {
    object@scores <- append(object@scores, score.list)
  }
  
  return(object)
}


calculateRelativeFC <- function(se, design, coef = NULL, contrast = NULL,
                                WTrows = NULL, selAssay = "counts",
                                pseudocount = 1, method = "edgeR",
                                normMethod = ifelse(is.null(WTrows),
                                                    "TMM", "sum")) {
  # .assertVector(x = se, type = "SummarizedExperiment")
  # .assertScalar(x = selAssay, type = "character")
  if (!(selAssay %in% SummarizedExperiment::assayNames(se))) {
    if (is.null(SummarizedExperiment::assayNames(se)) &&
        length(SummarizedExperiment::assays(se)) == 1) {
      warning("No assayNames provided in 'se', but only one ",
              "assay present - using that.")
      selAssay <- 1L
    } else {
      stop("The provided 'selAssay' not present in 'se'.")
    }
  }
  
  if (!is.null(WTrows)) {
    # .assertVector(x = WTrows, type = "character", validValues = rownames(se))
  }
  
  if (nrow(design) != ncol(se)) {
    stop("The number of rows in 'design' (", nrow(design),
         ") is not equal to the number",
         " of columns in 'se' (", ncol(se), ").")
  }
  
  # .assertScalar(x = pseudocount, type = "numeric", rngIncl = c(0, Inf))
  # 
  # .assertScalar(x = method, type = "character",
  #               validValues = c("edgeR", "limma"))
  
  if (normMethod %in% c("csaw", "TMM") && !is.null(WTrows)) {
    stop("normMethod = '", normMethod,
         "' can only be used when WTrows is NULL.")
  }
  
  if (normMethod %in% c("sum", "geomean") && is.null(WTrows)) {
    stop("normMethod = '", normMethod,
         "' can only be used when WTrows is not NULL.")
  }
  
  if (normMethod == "csaw" && method == "limma") {
    stop("normMethod = 'csaw' can only be used with method = 'edgeR'.")
  }
  
  # .assertScalar(x = normMethod, type = "character",
  #               validValues = c("csaw", "TMM", "geomean", "sum"))
  
  if (is.null(coef) && is.null(contrast)) {
    stop("'coef' and 'contrast' can not both be NULL.")
  }
  
  if (!is.null(contrast) && !is.null(dim(contrast))) {
    stop("'contrast' must be a vector.")
  }
  
  ## Create DGEList from SummarizedExperiment
  dge <- edgeR::DGEList(counts = as.matrix(SummarizedExperiment::assay(se, selAssay)),
                        samples = SummarizedExperiment::colData(se))
  if (normMethod == "csaw") {
    ## csaw normalization - also calculate normalization factors since
    ## aveLogCPM does not use provided offsets
    ## In this case, we know that WTrows is NULL, so all features
    ## will be used for the normalization
    dge <- edgeR::calcNormFactors(dge)
    dge <- csaw::normOffsets(dge)
  } else if (normMethod == "TMM") {
    ## TMM normalization, with all features
    dge <- edgeR::calcNormFactors(dge)
  } else if (normMethod == "geomean") {
    ## Use size factors (offsets) derived from the geometric mean
    ## of the WT rows
    tmp0 <- dge$counts[WTrows, , drop = FALSE]
    tmp0 <- tmp0[apply(tmp0, 1, min) > 0, , drop = FALSE]
    logoffsets <- apply(tmp0, 2, function(s) mean(log(s)))
    dge <- edgeR::scaleOffset(dge, logoffsets)
  } else if (normMethod == "sum") {
    ## Use size factors (offsets) derived from the sum of the
    ## WT rows
    tmp0 <- dge$counts[WTrows, , drop = FALSE]
    dge <- edgeR::scaleOffset(dge, log(colSums(tmp0)))
  }
  
  ## Fit model and perform test
  if (method == "edgeR") {
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmQLFit(dge, design = design)
    qlf <- edgeR::glmQLFTest(fit, coef = coef, contrast = contrast)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    ## Calculate shrunken fold changes. Only when testing
    ## a single coefficient or contrast
    predfc <- edgeR::predFC(dge, design = design, prior.count = pseudocount)
    if (length(coef) == 1 && is.null(contrast)) {
      tt$logFC_shrunk <- predfc[, coef]
    } else if (!is.null(contrast)) {
      tt$logFC_shrunk <- c(predfc %*% cbind(contrast))
    }
    tt$df.total <- qlf$df.total
    tt$df.prior <- qlf$df.prior
    tt$df.test <- qlf$df.test
  } else if (method == "limma") {
    if (!is.null(dge$offset)) {
      vm <- limma::voom(dge, design = design, lib.size = exp(dge$offset))
    } else {
      vm <- limma::voom(dge, design = design,
                        lib.size = edgeR::getNormLibSizes(dge))
    }
    fit <- limma::lmFit(vm, design = design)
    if (!is.null(contrast)) {
      fit <- limma::contrasts.fit(fit, contrasts = contrast)
      coef <- 1
    }
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = coef,
                          confint = TRUE, number = Inf, sort.by = "none")
    if (length(coef) == 1) {
      tt$se.logFC <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
    }
    tt$df.total <- fit$df.total
    tt$df.prior <- fit$df.prior
  }
  tt
}


library(tidyverse)
dir <- "data/OCT1/rosace"

# test1, test2, test3
sdir <- file.path(dir, "test3") # change here
load(file.path(sdir, "rosace.RData"))
count <- read_tsv("data/OCT1/count/1SM73/1SM73R3_sel/main_identifiers_counts_unfiltered.tsv") # change here


count_mat <- as.matrix(count[, -1])
count_mat[is.na(count_mat)] <- 0

se <- SummarizedExperiment::SummarizedExperiment(
  assay = list(count_mat),
  colData = data.frame(rep = rep(1, 3),
                       time = rep(0:2, 1)))
edger_scores <- calculateRelativeFC(
  se = se,
  design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
  coef = "time",
  pseudocount = 0.5,
  WTrows = "1",
  method = "edgeR"
)
edger_scores <- edger_scores[-1, ]
limma_scores <- calculateRelativeFC(
  se = se,
  design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
  coef = "time",
  pseudocount = 0.5,
  WTrows = "1", 
  method = "limma"
)
limma_scores <- limma_scores[-1, ]
df_edger <- data.frame(variant = count$hgvs[-1])
df_limma <- data.frame(variant = count$hgvs[-1])
df_edger <- cbind(df_edger, edger_scores %>% select(score = logFC_shrunk, fdr = FDR))
df_limma <- cbind(df_limma, limma_scores %>% select(score = logFC, fdr = adj.P.Val))
score_edger <- CreateScoreObject(method = "EDGER",
                                 type = "growth",
                                 assay.name = "1SM73",
                                 score = df_edger)
score_limma <- CreateScoreObject(method = "LIMMA",
                                 type = "growth",
                                 assay.name = "1SM73",
                                 score = df_limma)
rosace <- AddScoreData(rosace, score_edger)
rosace <- AddScoreData(rosace, score_limma)
save(rosace, file = file.path(sdir, "rosace_freq.RData"))


# test123
sdir <- file.path(dir, "test123")
load(file.path(sdir, "rosace.RData"))

count1 <- read_tsv("data/OCT1/count/1SM73/1SM73R1_sel/main_identifiers_counts_unfiltered.tsv")
colnames(count1) <- c("hgvs", "rep1.c0", "rep1.c1", "rep1.c2", "rep1.c3")
count2 <- read_tsv("data/OCT1/count/1SM73/1SM73R2_sel/main_identifiers_counts_unfiltered.tsv")
colnames(count2) <- c("hgvs", "rep2.c0", "rep2.c1", "rep2.c2", "rep2.c3")
count3 <- read_tsv("data/OCT1/count/1SM73/1SM73R3_sel/main_identifiers_counts_unfiltered.tsv")
colnames(count3) <- c("hgvs", "rep3.c0", "rep3.c1", "rep3.c2", "rep3.c3")

count <- inner_join(count1, count2) %>% inner_join(count3)
count_mat <- as.matrix(count[, -1])
count_mat[is.na(count_mat)] <- 0

se <- SummarizedExperiment::SummarizedExperiment(
  assay = list(count_mat),
  colData = data.frame(rep = c(rep(1, 4), rep(2, 4), rep(3, 3)),
                       time = c(rep(0:3, 2), 0:2))
  )
edger_scores <- calculateRelativeFC(
  se = se,
  design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
  coef = "time",
  pseudocount = 0.5,
  WTrows = "1",
  method = "edgeR"
)
edger_scores <- edger_scores[-1, ]
limma_scores <- calculateRelativeFC(
  se = se,
  design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
  coef = "time",
  pseudocount = 0.5,
  WTrows = "1", 
  method = "limma"
)
limma_scores <- limma_scores[-1, ]
df_edger <- data.frame(variant = count$hgvs[-1])
df_limma <- data.frame(variant = count$hgvs[-1])
df_edger <- cbind(df_edger, edger_scores %>% select(score = logFC_shrunk, fdr = FDR))
df_limma <- cbind(df_limma, limma_scores %>% select(score = logFC, fdr = adj.P.Val))

score_edger <- CreateScoreObject(method = "EDGER",
                                 type = "growth",
                                 assay.name = "1SM73",
                                 score = df_edger)
score_limma <- CreateScoreObject(method = "LIMMA",
                                 type = "growth",
                                 assay.name = "1SM73",
                                 score = df_limma)
rosace <- AddScoreData(rosace, score_edger)
rosace <- AddScoreData(rosace, score_limma)
save(rosace, file = file.path(sdir, "rosace_freq.RData"))

# test12
count1 <- read_tsv("data/OCT1/count/1SM73/1SM73R1_sel/main_identifiers_counts_unfiltered.tsv")
colnames(count1) <- c("hgvs", "rep1.c0", "rep1.c1", "rep1.c2", "rep1.c3")
count2 <- read_tsv("data/OCT1/count/1SM73/1SM73R2_sel/main_identifiers_counts_unfiltered.tsv")
colnames(count2) <- c("hgvs", "rep2.c0", "rep2.c1", "rep2.c2", "rep2.c3")
count3 <- read_tsv("data/OCT1/count/1SM73/1SM73R3_sel/main_identifiers_counts_unfiltered.tsv")
colnames(count3) <- c("hgvs", "rep3.c0", "rep3.c1", "rep3.c2")


sdir <- file.path(dir, "test12")
load(file.path(sdir, "rosace.RData"))
count <- inner_join(count1, count2)
count_mat <- as.matrix(count[, -1])
count_mat[is.na(count_mat)] <- 0

se <- SummarizedExperiment::SummarizedExperiment(
  assay = list(count_mat),
  colData = data.frame(rep = c(rep(1, 4), rep(2, 4)),
                       time = c(rep(0:3, 1), 0:3))
)
edger_scores <- calculateRelativeFC(
  se = se,
  design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
  coef = "time",
  pseudocount = 0.5,
  WTrows = "1",
  method = "edgeR"
)
edger_scores <- edger_scores[-1, ]
limma_scores <- calculateRelativeFC(
  se = se,
  design = model.matrix(~time, data = SummarizedExperiment::colData(se)),
  coef = "time",
  pseudocount = 0.5,
  WTrows = "1", 
  method = "limma"
)
limma_scores <- limma_scores[-1, ]
df_edger <- data.frame(variant = count$hgvs[-1])
df_limma <- data.frame(variant = count$hgvs[-1])
df_edger <- cbind(df_edger, edger_scores %>% select(score = logFC_shrunk, fdr = FDR))
df_limma <- cbind(df_limma, limma_scores %>% select(score = logFC, fdr = adj.P.Val))

score_edger <- CreateScoreObject(method = "EDGER",
                                 type = "growth",
                                 assay.name = "1SM73",
                                 score = df_edger)
score_limma <- CreateScoreObject(method = "LIMMA",
                                 type = "growth",
                                 assay.name = "1SM73",
                                 score = df_limma)
rosace <- AddScoreData(rosace, score_edger)
rosace <- AddScoreData(rosace, score_limma)
save(rosace, file = file.path(sdir, "rosace_freq.RData"))



     