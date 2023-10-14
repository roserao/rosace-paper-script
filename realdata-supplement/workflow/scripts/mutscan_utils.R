#!/usr/bin/env Rscript

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