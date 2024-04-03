#!/usr/bin/env Rscript
library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")

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

if (data == "OCT1") {
  eve <- FALSE
  aa <- "single"
  stop <- NA
  lof <- "pos"
} else if (data == "CARD11") {
  eve <- TRUE
  aa <- "triple"
  stop <- "Ter" # to * in the aa mapping
  lof <- "neg"
} else if (data == "MSH2") {
  eve <- TRUE
  aa <- "single"
  stop <- "*"
  lof <- "pos"
} else if (data == "BRCA1") {
  eve <- TRUE
  aa <- "single"
  stop <- "*"
  lof <- "neg"
} else if (data == "MET") {
  eve <- TRUE
  aa <- "single"
  stop <- "X"
  lof <- "neg"
} else if (data == "BRCA1-RING") {
  eve <- TRUE
  aa <- "triple"
  stop <- "Ter"
  lof <- "neg"
} else if (data == "Cohesin") {
  eve <- FALSE
  aa <- "triple"
  stop <- "Ter"
  lof <- "neg"
}

##### hyperparameter #####
model_vec <- c("ENRICH2", "ROSACE", "ROSACENP", "SLR", "EDGER", "LIMMA", "DIMSUM")
method_vec <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Rosace(nopos)", "Naive")
# color_vec <- c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#F0E442", "#56B4E9")
# linetype_vec <- c(2, 3, 1, 6, 4, 4, 5)
color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#ff7f0e', '#e377c2')
linetype_vec <- c(4, 6, 3, 4, 1, 2, 5)

##### plot save directory ##### 
if (mode == "N") {
    dir <- file.path("results", data)
} else {
    dir <- file.path("replicates", data, mode)
}
sdir <- file.path(dir, "plot")
if (!dir.exists(sdir)){
  dir.create(sdir)
}

##### effects dataframe ##### 
df <- read_tsv(file = file.path(dir, "analysis", "effects.tsv"))
if (aa == "triple") {
    load("data/sysdata.rda")
    df <- df %>% 
        left_join(aa_table, c("wildtype" = "triple")) %>%
        select(-wildtype) %>%
        rename(wildtype = single) %>% 
        left_join(aa_table, c("mutation" = "triple")) %>%
        select(-mutation) %>%
        rename(mutation = single) 
} else {
    if (!is.na(stop) && stop != "*") {
        df$mutation[df$mutation == stop] <- "*"
    }
}
if (data == "MET") {
    df$position <- 1058 + df$position
}

if (data != "Cohesin") {
    df_gnomAD  <- read_tsv(file = file.path("data", data, "gnomAD_cleaned.tsv")) %>%
        select(position, wildtype, mutation, clinvar) # clinvar
    df <- df %>% left_join(df_gnomAD)  

    df_alphaM <-  read_tsv(file = file.path("data", data, "AlphaM_cleaned.tsv")) %>%
        select(position, wildtype, mutation, alpham)
    df <- df %>% left_join(df_alphaM)
}

if (eve) {
    df_eve  <- read_tsv(file = file.path("data", data, "EVE_cleaned.tsv")) # EVE.1
    df <- df %>% left_join(df_eve)
}

##### synonymous mutation rank plot ##### 

# rank fdr test for each simulation
comp_rankfdr <- function(df, resolution, model, rank_by_effect = TRUE) {
    df_ranked <- df %>% 
        select(variants, type,
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) %>%
        mutate(effect = abs(effect)) 

    if (rank_by_effect){
        df_ranked <- df_ranked %>% arrange(desc(effect), test) 
    }  else {
        df_ranked <- df_ranked %>% arrange(test, desc(effect)) 
    }
    
    n_syn <- sum(df_ranked$type == "synonymous", na.rm = TRUE)
    cutoff <- floor(nrow(df_ranked)/resolution * (1:resolution))
    result <- data.frame(rank = 1:resolution, fdr = -0.1)
    for (i in 1:resolution) {
        df_ranked_sig <- df_ranked[1:cutoff[i], ] 
        result[i, 2] <- sum(df_ranked_sig$type == "synonymous", na.rm = TRUE)/n_syn
    }
    result$method <- model
    return(result)
}

# RANK BY EFFECT: run function on each method
rank_pform <- data.frame()
for (mod in model_vec) {
  res <- comp_rankfdr(df, 500, mod)
  rank_pform <- rbind(rank_pform, res)
}
rank_pform$method <- as.factor(rank_pform$method)
levels(rank_pform$method) <- method_vec
write_tsv(rank_pform, file = file.path(dir, "analysis", "rank_syn_fdr.tsv"))

# RANK BY EFFECT: plot
p1a <- ggplot(rank_pform, aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL, title = data)
cowplot::save_plot(p1a, file = file.path(sdir, "rank_fdr_syn.png"), base_width = 6, base_height = 4)
rm(rank_pform)

# RANK BY TEST STATISTICS: run function on each method
rank_pform <- data.frame()
for (mod in model_vec) {
  res <- comp_rankfdr(df, 500, mod, rank_by_effect = FALSE)
  rank_pform <- rbind(rank_pform, res)
}
rank_pform$method <- as.factor(rank_pform$method)
levels(rank_pform$method) <- method_vec
write_tsv(rank_pform, file = file.path(dir, "analysis", "rank_syn_fdr_test.tsv"))

# RANK BY TEST STATISTICS: plot
p1b <- ggplot(rank_pform, aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL, title = data)
cowplot::save_plot(p1b, file = file.path(sdir, "rank_fdr_syn_test.png"), base_width = 6, base_height = 4)
rm(rank_pform)

##### TP/TN test new version for EVE ##### 
if (eve) {
    eve_rank_ROC <- function(df, resolution, model, rank_by_effect = TRUE) {
        df_ranked <- df %>% 
            select(variants, type, EVE.1, 
                effect = starts_with(paste(model, "effects", sep = "_")),
                test = starts_with(paste(model, "tests", sep = "_"))) 
        #    mutate(effect = abs(effect))

        # if (rank_by_effect){
        #     df_ranked <- df_ranked %>% arrange(desc(effect), test) 
        # }  else {
        #     df_ranked <- df_ranked %>% arrange(test, desc(effect)) 
        # }
        if (lof == "neg") {
            df_ranked <- df_ranked %>% 
                mutate(test_dir = ifelse(effect < 0, -(1 - test), test), 
                    effect_dir = effect) 
        } else if (lof == "pos") {
            df_ranked <- df_ranked %>% 
                mutate(test_dir = ifelse(effect > 0, -(1 - test), test), 
                    effect_dir = -effect) 
        }
        if (rank_by_effect) {
            df_ranked <- df_ranked %>% arrange(effect_dir, test_dir) 
        } else {
            df_ranked <- df_ranked %>% arrange(test_dir, effect_dir) 
        }

        n_ctrl <- sum(df$EVE.1 == "Benign", na.rm = TRUE)
        n_sig <- sum(df$EVE.1 == "Pathogenic", na.rm = TRUE)
        
        cutoff <- floor(nrow(df_ranked)/resolution * (1:resolution))
        result <- data.frame(rank = (1:resolution)/resolution, FP = -1, TP = -1, 
                            n_ctrl = n_ctrl, n_sig = n_sig)
        for (i in 1:resolution) {
            df_ranked_sub <- df_ranked[1:cutoff[i], ] 
            result[i, 2] <- sum(df_ranked_sub$EVE.1 == "Benign", na.rm = TRUE)/n_ctrl
            if (n_sig != 0) {
                result[i, 3] <- sum(df_ranked_sub$EVE.1 == "Pathogenic", na.rm = TRUE)/n_sig
            }
        }
        result$method <- model
        return(result)
    }
    
    # RANK BY EFFECT: run function on each method
    eve_roc <- data.frame()
    for (mod in model_vec) {
        res <- eve_rank_ROC(df, 500, mod)
        eve_roc <- rbind(eve_roc, res)
    }
    eve_roc$method <- as.factor(eve_roc$method)
    levels(eve_roc$method) <- method_vec
    write_tsv(eve_roc, file = file.path(dir, "analysis", "eve_roc_rank.tsv"))

    # RANK BY EFFECT: plot
    p2a <- ggplot(eve_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (EVE10pct Benign)", y = "TP (EVE10pct Pathogenic)", color = NULL, shape = "rank proportion", 
            title = data, subtitle = paste("nBenign = ", eve_roc$n_ctrl[1], " nPathogenic = ", eve_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = eve_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p2a, file = file.path(sdir, "eve_roc_rank.png"), base_width = 6, base_height = 4)
    rm(eve_roc)

    # RANK BY TEST STATISTICS: run function on each method
    eve_roc <- data.frame()
    for (mod in model_vec) {
        res <- eve_rank_ROC(df, 500, mod, rank_by_effect = FALSE)
        eve_roc <- rbind(eve_roc, res)
    }
    eve_roc$method <- as.factor(eve_roc$method)
    levels(eve_roc$method) <- method_vec
    write_tsv(eve_roc, file = file.path(dir, "analysis", "eve_roc_rank_test.tsv"))

    # RANK BY TEST STATISTICS: plot
    p2b <- ggplot(eve_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (EVE10pct Benign)", y = "TP (EVE10pct Pathogenic)", color = NULL, shape = "rank proportion", 
            title = data, subtitle = paste("nBenign = ", eve_roc$n_ctrl[1], " nPathogenic = ", eve_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = eve_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p2b, file = file.path(sdir, "eve_roc_rank_test.png"), base_width = 6, base_height = 4)
}


##### TP/TN test for EVE ##### 
if (eve) {
    # EVE fdr test for each simulation
    eve_ROC <- function(df, model) {
        df <- df %>% 
            select(variants, type, EVE.1,
                effect = starts_with(paste(model, "effects", sep = "_")),
                test = starts_with(paste(model, "tests", sep = "_"))) %>%
            filter(!is.na(EVE.1)) 

        df_ctrl <- df %>% filter(EVE.1 == "Benign")
        df_sig <- df %>% filter(EVE.1 == "Pathogenic")
        n_ctrl <- sum(df$EVE.1 == "Benign", na.rm = TRUE)
        n_sig <- sum(df$EVE.1 == "Pathogenic", na.rm = TRUE)

        cutoff <- seq(0, 0.2, by = 0.01)
        result <- data.frame(test_cutoff = cutoff, FP = -1, TP = -1, 
                             n_ctrl = n_ctrl, n_sig = n_sig)
        for (i in 1:length(cutoff)) {
            result[i, 2] <- sum(df_ctrl$test <= cutoff[i], na.rm = TRUE)/n_ctrl
            if (n_sig != 0) {
                result[i, 3] <- sum(df_sig$test <= cutoff[i], na.rm = TRUE)/n_sig
            }
        }
        result$method <- model
        return(result)
    }

    # run function on each method
    eve_roc <- data.frame()
    for (mod in model_vec) {
        res <- eve_ROC(df, mod)
        eve_roc <- rbind(eve_roc, res)
    }
    eve_roc$method <- as.factor(eve_roc$method)
    levels(eve_roc$method) <- method_vec
    write_tsv(eve_roc, file = file.path(dir, "analysis", "eve_roc.tsv"))

    p2c <- ggplot(eve_roc %>% filter(test_cutoff != 0), aes(x = FP, y = TP)) +
    geom_path(aes(color = method), alpha = 0.7) +
    scale_color_manual(values = color_vec) +
    scale_linetype_manual(values = linetype_vec) +
    labs(x = "FP (EVE10pct Benign)", y = "TP (EVE10pct Pathogenic)", color = NULL, shape = "test cutoff", 
         title = data, subtitle = paste("nBenign = ", eve_roc$n_ctrl[1], " nPathogenic = ", eve_roc$n_sig[1], sep = "")) +
    cowplot::theme_cowplot() +
    geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(test_cutoff)), 
                data = eve_roc %>% filter(test_cutoff %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p2c, file = file.path(sdir, "eve_roc.png"), base_width = 6, base_height = 4)
}

##### TP/TN test rank for alpham ##### 
alpham_rank_ROC <- function(df, resolution, model, rank_by_effect = TRUE) {
    df_ranked <- df %>% 
        select(variants, type, alpham, 
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) 
       # mutate(effect_abs = abs(effect))
    
    # if (rank_by_effect){
    #     df_ranked <- df_ranked %>% arrange(desc(effect_abs), test) 
    # }  else {
    #     df_ranked <- df_ranked %>% arrange(test, desc(effect_abs))
    # }
    
    if (lof == "neg") {
        df_ranked <- df_ranked %>% 
            mutate(test_dir = ifelse(effect < 0, -(1 - test), test), 
                   effect_dir = effect) 
    } else if (lof == "pos") {
        df_ranked <- df_ranked %>% 
            mutate(test_dir = ifelse(effect > 0, -(1 - test), test), 
                   effect_dir = -effect) 
    }
    if (rank_by_effect) {
        df_ranked <- df_ranked %>% arrange(effect_dir, test_dir) 
    } else {
        df_ranked <- df_ranked %>% arrange(test_dir, effect_dir) 
    }

    n_ctrl <- sum(df$alpham %in% c("likely_benign"), na.rm = TRUE)
    n_sig <- sum(df$alpham %in% c("likely_pathogenic"), na.rm = TRUE)
    
    cutoff <- floor(nrow(df_ranked)/resolution * (1:resolution))
    result <- data.frame(rank = (1:resolution)/resolution, FP = -1, TP = -1, 
                        n_ctrl = n_ctrl, n_sig = n_sig)
    for (i in 1:resolution) {
        df_ranked_sub <- df_ranked[1:cutoff[i], ] 
        result[i, 2] <- sum(df_ranked_sub$alpham %in% c("likely_benign"), na.rm = TRUE)/n_ctrl
        if (n_sig != 0) {
            result[i, 3] <- sum(df_ranked_sub$alpham %in% c("likely_pathogenic"), na.rm = TRUE)/n_sig
        }
    }
    result$method <- model
    return(result)
}

if (data != "Cohesin") {
    # RANK BY EFFECT: run function on each method
    alpham_roc <- data.frame()
    for (mod in model_vec) {
        res <- alpham_rank_ROC(df, 500, mod, rank_by_effect = TRUE)
        alpham_roc <- rbind(alpham_roc, res)
    }
    alpham_roc$method <- as.factor(alpham_roc$method)
    levels(alpham_roc$method) <- method_vec
    write_tsv(alpham_roc, file = file.path(dir, "analysis", "alpham_roc_rank.tsv"))

    # RANK BY EFFECT: plot
    p3a <- ggplot(alpham_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (alphaM Benign)", y = "TP (alphaM Pathogenic)", color = NULL, shape = "rank", 
            title = data, subtitle = paste("nBenign = ", alpham_roc$n_ctrl[1], " nPathogenic = ", alpham_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                    data = alpham_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3a, file = file.path(sdir, "alpham_roc_rank.png"), base_width = 6, base_height = 4)

    # RANK BY TEST STATISTICS: run function on each method
    alpham_roc <- data.frame()
    for (mod in model_vec) {
        res <- alpham_rank_ROC(df, 500, mod, rank_by_effect = FALSE)
        alpham_roc <- rbind(alpham_roc, res)
    }
    alpham_roc$method <- as.factor(alpham_roc$method)
    levels(alpham_roc$method) <- method_vec
    write_tsv(alpham_roc, file = file.path(dir, "analysis", "alpham_roc_rank_test.tsv"))

    # RANK BY TEST STATISTICS: plot
    p3b <- ggplot(alpham_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (alphaM Benign)", y = "TP (alphaM Pathogenic)", color = NULL, shape = "rank", 
            title = data, subtitle = paste("nBenign = ", alpham_roc$n_ctrl[1], " nPathogenic = ", alpham_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                    data = alpham_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3b, file = file.path(sdir, "alpham_roc_rank_test.png"), base_width = 6, base_height = 4)

}

##### TP/TN test for AlphaM ##### 
# AlphaM fdr test for each simulation
alpham_ROC <- function(df, model) {
    df <- df %>% 
        select(variants, alpham,
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) %>%
        filter(!is.na(alpham)) 

    df_ctrl <- df %>% filter(alpham == "likely_benign")
    df_sig <- df %>% filter(alpham == "likely_pathogenic")
    n_ctrl <- sum(df$alpham == "likely_benign", na.rm = TRUE)
    n_sig <- sum(df$alpham == "likely_pathogenic", na.rm = TRUE)

    cutoff <- seq(0, 0.2, by = 0.01)
    result <- data.frame(test_cutoff = cutoff, FP = -1, TP = -1, 
                            n_ctrl = n_ctrl, n_sig = n_sig)
    for (i in 1:length(cutoff)) {
        result[i, 2] <- sum(df_ctrl$test <= cutoff[i], na.rm = TRUE)/n_ctrl
        if (n_sig != 0) {
            result[i, 3] <- sum(df_sig$test <= cutoff[i], na.rm = TRUE)/n_sig
        }
    }
    result$method <- model
    return(result)
}

# run function on each method
if (data != "Cohesin") {
    alpham_roc <- data.frame()
    for (mod in model_vec) {
        res <- alpham_ROC(df, mod)
        alpham_roc <- rbind(alpham_roc, res)
    }
    alpham_roc$method <- as.factor(alpham_roc$method)
    levels(alpham_roc$method) <- method_vec
    write_tsv(alpham_roc, file = file.path(dir, "analysis", "alpham_roc.tsv"))

    p2c <- ggplot(alpham_roc %>% filter(test_cutoff != 0), aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (alphaM Benign)", y = "TP (alphaM Pathogenic)", color = NULL, shape = "test cutoff", 
                title = data, subtitle = paste("nBenign = ", alpham_roc$n_ctrl[1], " nPathogenic = ", alpham_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(test_cutoff)), 
                    data = alpham_roc %>% filter(test_cutoff %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p2c, file = file.path(sdir, "alpham_roc.png"), base_width = 6, base_height = 4)
}

##### TP/TN test rank for genomAD ##### 
gnomad_rank_ROC <- function(df, resolution, model, rank_by_effect = TRUE) {
    df_ranked <- df %>% 
        select(variants, type, clinvar, 
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) 
       # mutate(effect_abs = abs(effect))
    
    # if (rank_by_effect){
    #     df_ranked <- df_ranked %>% arrange(desc(effect_abs), test) 
    # }  else {
    #     df_ranked <- df_ranked %>% arrange(test, desc(effect_abs))
    # }
    
    if (lof == "neg") {
        df_ranked <- df_ranked %>% 
            mutate(test_dir = ifelse(effect < 0, -(1 - test), test), 
                   effect_dir = effect) 
    } else if (lof == "pos") {
        df_ranked <- df_ranked %>% 
            mutate(test_dir = ifelse(effect > 0, -(1 - test), test), 
                   effect_dir = -effect) 
    }
    if (rank_by_effect) {
        df_ranked <- df_ranked %>% arrange(effect_dir, test_dir) 
    } else {
        df_ranked <- df_ranked %>% arrange(test_dir, effect_dir) 
    }

    n_ctrl <- sum(df$clinvar %in% c("Benign", "Likely benign"), na.rm = TRUE)
    n_sig <- sum(df$clinvar %in% c("Pathogenic", "Likely pathogenic"), na.rm = TRUE)
    
    cutoff <- floor(nrow(df_ranked)/resolution * (1:resolution))
    result <- data.frame(rank = (1:resolution)/resolution, FP = -1, TP = -1, 
                        n_ctrl = n_ctrl, n_sig = n_sig)
    for (i in 1:resolution) {
        df_ranked_sub <- df_ranked[1:cutoff[i], ] 
        result[i, 2] <- sum(df_ranked_sub$clinvar %in% c("Benign", "Likely benign"), na.rm = TRUE)/n_ctrl
        if (n_sig != 0) {
            result[i, 3] <- sum(df_ranked_sub$clinvar %in% c("Pathogenic", "Likely pathogenic"), na.rm = TRUE)/n_sig
        }
    }
    result$method <- model
    return(result)
}

if (data != "Cohesin") {
    # RANK BY EFFECT: run function on each method
    gnomad_roc <- data.frame()
    for (mod in model_vec) {
    res <- gnomad_rank_ROC(df, 500, mod, rank_by_effect = TRUE)
    gnomad_roc <- rbind(gnomad_roc, res)
    }
    gnomad_roc$method <- as.factor(gnomad_roc$method)
    levels(gnomad_roc$method) <- method_vec
    write_tsv(gnomad_roc, file = file.path(dir, "analysis", "gnomad_roc_rank.tsv"))

    # RANK BY EFFECT: plot
    p3a <- ggplot(gnomad_roc, aes(x = FP, y = TP)) +
    geom_path(aes(color = method), alpha = 0.7) +
    scale_color_manual(values = color_vec) +
    scale_linetype_manual(values = linetype_vec) +
    labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "rank", 
        title = data, subtitle = paste("nBenign = ", gnomad_roc$n_ctrl[1], " nPathogenic = ", gnomad_roc$n_sig[1], sep = "")) +
    cowplot::theme_cowplot() +
    geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = gnomad_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3a, file = file.path(sdir, "gnomad_roc_rank.png"), base_width = 6, base_height = 4)

    # RANK BY TEST STATISTICS: run function on each method
    gnomad_roc <- data.frame()
    for (mod in model_vec) {
    res <- gnomad_rank_ROC(df, 500, mod, rank_by_effect = FALSE)
    gnomad_roc <- rbind(gnomad_roc, res)
    }
    gnomad_roc$method <- as.factor(gnomad_roc$method)
    levels(gnomad_roc$method) <- method_vec
    write_tsv(gnomad_roc, file = file.path(dir, "analysis", "gnomad_roc_rank_test.tsv"))

    # RANK BY TEST STATISTICS: plot
    p3b <- ggplot(gnomad_roc, aes(x = FP, y = TP)) +
    geom_path(aes(color = method), alpha = 0.7) +
    scale_color_manual(values = color_vec) +
    scale_linetype_manual(values = linetype_vec) +
    labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "rank", 
        title = data, subtitle = paste("nBenign = ", gnomad_roc$n_ctrl[1], " nPathogenic = ", gnomad_roc$n_sig[1], sep = "")) +
    cowplot::theme_cowplot() +
    geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = gnomad_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3b, file = file.path(sdir, "gnomad_roc_rank_test.png"), base_width = 6, base_height = 4)

}


##### TP/TN test for genomAD ##### 

# gnomad fdr test for each simulation
gnomad_ROC <- function(df, model) {
    df <- df %>% 
        select(variants, type, clinvar,
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) %>%
        filter(!is.na(clinvar)) %>%
        drop_na()
    df_ctrl <- df %>% filter(clinvar %in% c("Benign", "Likely benign"))
    df_sig <- df %>% filter(clinvar %in% c("Pathogenic", "Likely pathogenic"))

    cutoff <- seq(0, 0.2, by = 0.01)
    result <- data.frame(test_cutoff = cutoff, FP = -1, TP = -1, 
                         n_ctrl = nrow(df_ctrl), n_sig = nrow(df_sig))
    for (i in 1:length(cutoff)) {
        result[i, 2] <- sum(df_ctrl$test <= cutoff[i], na.rm = TRUE)/nrow(df_ctrl)
        if (nrow(df_sig) != 0) {
            result[i, 3] <- sum(df_sig$test <= cutoff[i], na.rm = TRUE)/nrow(df_sig)
        }
    }
    result$method <- model
    return(result)
}

if (data != "Cohesin") {
    # run function on each method
    gnomad_roc <- data.frame()
    for (mod in model_vec) {
    res <- gnomad_ROC(df, mod)
    gnomad_roc <- rbind(gnomad_roc, res)
    }
    gnomad_roc$method <- as.factor(gnomad_roc$method)
    levels(gnomad_roc$method) <- method_vec
    write_tsv(gnomad_roc, file = file.path(dir, "analysis", "gnomad_roc.tsv"))

    p3c <- ggplot(gnomad_roc %>% filter(test_cutoff != 0), aes(x = FP, y = TP)) +
    geom_path(aes(color = method), alpha = 0.7) +
    scale_color_manual(values = color_vec) +
    scale_linetype_manual(values = linetype_vec) +
    labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "test cutoff", 
        title = data, subtitle = paste("nBenign = ", gnomad_roc$n_ctrl[1], " nPathogenic = ", gnomad_roc$n_sig[1], sep = "")) +
    cowplot::theme_cowplot() +
    geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(test_cutoff)), 
                data = gnomad_roc %>% filter(test_cutoff %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3c, file = file.path(sdir, "gnomad_roc.png"), base_width = 6, base_height = 4)
}


##### correlation plot of ranking of synonymous and stop codon ##### 
# rosace vs enrich2
df_plot <- df %>% 
    filter(variants != "_wt") %>%
    select(type, ROSACE_effects, ENRICH2_effects) %>% rowwise() %>%
    mutate(syn = ifelse(type == "synonymous", "synonymous", "others")) %>%
    mutate(stop_indel = ifelse(type != "synonymous" & type != "missense", "stop/indel", "others")) 

cor_syn <- cor(df_plot$ROSACE_effects[df_plot$syn == "synonymous"], df_plot$ENRICH2_effects[df_plot$syn == "synonymous"], 
    use = "complete.obs", method = "spearman")
cor_nonsyn <- cor(df_plot$ROSACE_effects[df_plot$syn != "synonymous"], df_plot$ENRICH2_effects[df_plot$syn != "synonymous"], 
    use = "complete.obs", method = "spearman")
cor_all <- cor(df_plot$ROSACE_effects, df_plot$ENRICH2_effects, 
    use = "complete.obs", method = "spearman")

p4 <- ggplot(df_plot, aes(ROSACE_effects, ENRICH2_effects)) +
    geom_point(aes(color = syn, alpha = syn, size = syn)) +
    scale_color_manual(values = c("grey", "red")) +
    scale_size_manual(values = c(0.1, 0.1)) +
    scale_alpha_manual(values = c(0.2, 0.9)) +
    labs(title = data) +
    guides(alpha = "none", size = "none") +
    theme_classic() +
    annotate("text", x = min(df_plot$ROSACE_effects, na.rm = TRUE)+2, 
        y = max(df_plot$ENRICH2_effects, na.rm = TRUE)-1, size = 3,
        label = paste("spearman corr\nw/ and w/o syn:\n", round(cor_all, 5), " and ", round(cor_nonsyn, 5), sep = ""))
cowplot::save_plot(p4, file = file.path(sdir, "syn_corr.png"), base_width = 6, base_height = 4)

p5 <- ggplot(df_plot[-1, ], aes(ROSACE_effects, ENRICH2_effects)) +
    geom_point(aes(color = stop_indel, alpha = stop_indel, size = stop_indel)) +
    scale_color_manual(values = c("grey", "blue")) +
    scale_size_manual(values = c(0.1, 0.1)) +
    scale_alpha_manual(values = c(0.2, 0.9)) +
    labs(title = data) +
    guides(alpha = "none", size = "none") +
    theme_classic()
cowplot::save_plot(p5, file = file.path(sdir, "stop_corr.png"), base_width = 6, base_height = 4)
rm(df_plot)

##### MET stop codon illustration #####
if (data == "MET" && mode == "N") {
    df_stop <- df %>% 
        filter(variants != "_wt") %>%
        filter(position %in% df$position[df$mutation == "*"]) %>%
        filter(position <= 1130) %>%
        mutate(stop = ifelse(mutation == "*", "stop", "others"))
    
    p6 <- ggplot(df_stop, aes(ROSACE_effects, ENRICH2_effects)) +
        geom_point(aes(color = as.factor(position), shape = stop), size = 2) +
        labs(title = data, color = "position", shape = "stop") +
        cowplot::theme_cowplot() +
        facet_wrap(vars(position), nrow = 2)

    cowplot::save_plot(p6, file = file.path(sdir, "stop_position.png"), base_width = 8, base_height = 6)
}

##### positive control set #####

exp_roc_replicate <- function(df, resolution, model, pos, rank_by_effect = TRUE) {
    df_ranked <- df %>% 
        select(variants, type,
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) %>%
        mutate(effect = abs(effect)) 

    if (rank_by_effect){
        df_ranked <- df_ranked %>% arrange(desc(effect), test) 
    }  else {
        df_ranked <- df_ranked %>% arrange(test, desc(effect)) 
    }
    
    n_syn <- sum(df_ranked$type == "synonymous", na.rm = TRUE)
    cutoff <- floor(nrow(df_ranked)/resolution * (1:resolution))
    result <- data.frame(rank = 1:resolution, FP = -0.1, TP = -0.1)
    for (i in 1:resolution) {
        df_ranked_sig <- df_ranked[1:cutoff[i], ] 
        result[i, 2] <- sum(df_ranked_sig$type == "synonymous", na.rm = TRUE)/n_syn
        result[i, 3] <- sum(df_ranked_sig$variants %in% pos, na.rm = TRUE)/length(pos)
    }
    result$method <- model
    return(result)
}

exp_power_replicate <- function(df, model, pos) {
    df <- df %>% 
        select(variants, type, 
            effect = starts_with(paste(model, "effects", sep = "_")),
            test = starts_with(paste(model, "tests", sep = "_"))) %>%
        filter(variants %in% pos) %>%
        drop_na()
    
    cutoff <- seq(0, 0.2, by = 0.01)
    result <- data.frame(test_cutoff = cutoff, power = -1)
    for (i in 1:length(cutoff)) {
        result[i, 2] <- sum(df$test <= cutoff[i], na.rm = TRUE)/length(pos)
    }
    result$method <- model
    return(result)
}


if (data == "OCT1") {
    pos <- c("p.(E284A)", "p.(W64R)", "p.(R175H)", "p.(D149N)",
        "p.(D149R)", "p.(R486W)", "p.(E386K)", 
        "p.(D303E)", "p.(D303G)", "p.(Q241N)")
} 

if (mode != "N" && data == "OCT1") {

    # RANK BY EFFECT: run function on each method
    exp_roc <- data.frame()
    for (mod in model_vec) {
        res <- exp_roc_replicate(df, 500, mod, pos = pos, rank_by_effect = TRUE)
        exp_roc <- rbind(exp_roc, res)
    }
    exp_roc$method <- as.factor(exp_roc$method)
    levels(exp_roc$method) <- method_vec
    write_tsv(exp_roc, file = file.path(dir, "analysis", "exp_roc_rank.tsv"))

    # RANK BY EFFECT: plot
    p7a <- ggplot(exp_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (synonymous)", y = "TP (experiment)", color = NULL, shape = "rank proportion", 
            title = data) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = exp_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p7a, file = file.path(sdir, "exp_roc_rank.png"), base_width = 6, base_height = 4)
    rm(exp_roc)

    # RANK BY TEST STATISTICS: run function on each method
    exp_roc <- data.frame()
    for (mod in model_vec) {
        res <- exp_roc_replicate(df, 500, mod, pos = pos, rank_by_effect = FALSE)
        exp_roc <- rbind(exp_roc, res)
    }
    exp_roc$method <- as.factor(exp_roc$method)
    levels(exp_roc$method) <- method_vec
    write_tsv(exp_roc, file = file.path(dir, "analysis", "exp_roc_rank_test.tsv"))

    # RANK BY TEST STATISTICS: plot
    p7b <- ggplot(exp_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (synonymous)", y = "TP (experiment)", color = NULL, shape = "rank proportion", 
            title = data) +
        cowplot::theme_cowplot() +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = exp_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p7b, file = file.path(sdir, "exp_roc_rank_test.png"), base_width = 6, base_height = 4)
    rm(exp_roc)

    # power by experiment validation
    exp_power <- data.frame()
    for (mod in model_vec) {
        res <- exp_power_replicate(df, mod, pos = pos)
        exp_power <- rbind(exp_power, res)
    }
    exp_power$method <- as.factor(exp_power$method)
    levels(exp_power$method) <- method_vec
    write_tsv(exp_power, file = file.path(dir, "analysis", "exp_power.tsv"))
}



