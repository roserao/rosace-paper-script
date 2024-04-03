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
        help = "target protein of the experiment")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data

v_mode <- c('rep1_1', 'rep1_2', 'rep1_3', 'rep2_1', 'rep2_2', 'rep2_3', 'rep3_1')
v_rep <- c(1, 1, 1, 2, 2, 3)
v_rep_text <- c("number of replicate = 1", "number of replicate = 2", "number of replicate = 3")
method_vec <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Rosace(nopos)", "Naive")
color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#ff7f0e', '#e377c2')
linetype_vec <- c(4, 6, 3, 4, 1, 2, 5)

#####  directory ##### 
dir <- file.path("replicates", data)
sdir <- file.path("replicates", data, "summary")
if (!dir.exists(sdir)){
  dir.create(sdir)
}

### Both MET and OCT1

##### plot1A: replicate level rank_syn_fdr  # 3 figures 1 row
# file: rank_syn_fdr.tsv
rank_syn_fdr <- data.frame()
for (i in 1:6) {
    df <- read_tsv(file.path(dir, v_mode[i], "analysis", "rank_syn_fdr.tsv"))
    df$rep <- v_rep[i]
    rank_syn_fdr <- rbind(rank_syn_fdr, df)
}
rank_syn_fdr <- rank_syn_fdr %>% group_by(rank, rep, method) %>%
    summarise(fdr = mean(fdr, na.rm = TRUE))

rank_syn_fdr$rep <- as.factor(rank_syn_fdr$rep)
levels(rank_syn_fdr$rep) <- v_rep_text
rank_syn_fdr$method <-  factor(rank_syn_fdr$method, levels = method_vec)
write_tsv(rank_syn_fdr, file = file.path(sdir, "rank_syn_fdr.tsv"))

p1a <- ggplot(rank_syn_fdr, aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  facet_wrap(vars(rep)) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL, title = data) 
cowplot::save_plot(p1a, file = file.path(sdir, "rank_fdr_syn.png"), base_width = 15, base_height = 4.5)
rm(rank_syn_fdr)

##### plot1B: replicate level rank_syn_fdr_test (same as figure 3A) # 3 figures 1 row
# file: rank_syn_fdr_test.tsv
rank_syn_fdr <- data.frame()
for (i in 1:6) {
    df <- read_tsv(file.path(dir, v_mode[i], "analysis", "rank_syn_fdr_test.tsv"))
    df$rep <- v_rep[i]
    rank_syn_fdr <- rbind(rank_syn_fdr, df)
}
rank_syn_fdr <- rank_syn_fdr %>% group_by(rank, rep, method) %>%
    summarise(fdr = mean(fdr, na.rm = TRUE))

rank_syn_fdr$rep <- as.factor(rank_syn_fdr$rep)
levels(rank_syn_fdr$rep) <- v_rep_text
rank_syn_fdr$method <-  factor(rank_syn_fdr$method, levels = method_vec)
write_tsv(rank_syn_fdr, file = file.path(sdir, "rank_syn_fdr_test.tsv"))

p1a <- ggplot(rank_syn_fdr, aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  facet_wrap(vars(rep)) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL, title = data) 
cowplot::save_plot(p1a, file = file.path(sdir, "rank_fdr_syn_test.png"), base_width = 15, base_height = 4.5)
rm(rank_syn_fdr)

##### plot2A: replicate level gnomad_roc_rank # 3 figures 1 row
# file: gnomad_roc_rank.tsv
gnomad_roc <- data.frame()
for (i in 1:6) {
    df <- read_tsv(file.path(dir, v_mode[i], "analysis", "gnomad_roc_rank.tsv"))
    df$rep <- v_rep[i]
    gnomad_roc <- rbind(gnomad_roc, df)
}
gnomad_roc <- gnomad_roc %>% group_by(rank, rep, method) %>%
    summarise(FP = mean(FP, na.rm = TRUE),
              TP = mean(TP, na.rm = TRUE),
              n_ctrl = max(n_ctrl, na.rm = TRUE), 
              n_sig = max(n_sig, na.rm = TRUE))
gnomad_roc$rep <- as.factor(gnomad_roc$rep)
levels(gnomad_roc$rep) <- v_rep_text
gnomad_roc$method <-  factor(gnomad_roc$method, levels = method_vec)
write_tsv(gnomad_roc, file = file.path(sdir, "gnomad_roc_rank.tsv"))

p2a <- ggplot(gnomad_roc, aes(x = FP, y = TP)) +
  geom_path(aes(color = method), alpha = 0.7) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "rank", 
       title = data, subtitle = paste("nBenign = ", gnomad_roc$n_ctrl[1], " nPathogenic = ", gnomad_roc$n_sig[1], sep = "")) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(rep)) +
  geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
             data = gnomad_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
cowplot::save_plot(p2a, file = file.path(sdir, "gnomad_roc_rank.png"), base_width = 15, base_height = 4.5)
rm(gnomad_roc)

##### plot2B: replicate level gnomad_roc_rank_test # 3 figures 1 row
# file: gnomad_roc_rank_test.tsv
gnomad_roc <- data.frame()
for (i in 1:6) {
    df <- read_tsv(file.path(dir, v_mode[i], "analysis", "gnomad_roc_rank_test.tsv"))
    df$rep <- v_rep[i]
    gnomad_roc <- rbind(gnomad_roc, df)
}
gnomad_roc <- gnomad_roc %>% group_by(rank, rep, method) %>%
    summarise(FP = mean(FP, na.rm = TRUE),
              TP = mean(TP, na.rm = TRUE),
              n_ctrl = max(n_ctrl, na.rm = TRUE), 
              n_sig = max(n_sig, na.rm = TRUE))
gnomad_roc$rep <- as.factor(gnomad_roc$rep)
levels(gnomad_roc$rep) <- v_rep_text
gnomad_roc$method <-  factor(gnomad_roc$method, levels = method_vec)
write_tsv(gnomad_roc, file = file.path(sdir, "gnomad_roc_rank_test.tsv"))

p2b <- ggplot(gnomad_roc, aes(x = FP, y = TP)) +
  geom_path(aes(color = method), alpha = 0.7) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "rank", 
       title = data, subtitle = paste("nBenign = ", gnomad_roc$n_ctrl[1], " nPathogenic = ", gnomad_roc$n_sig[1], sep = "")) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(rep)) +
  geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
             data = gnomad_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
cowplot::save_plot(p2b, file = file.path(sdir, "gnomad_roc_rank_test.png"), base_width = 15, base_height = 4.5)
rm(gnomad_roc)

##### plotxA: replicate level alpham_roc_rank # 3 figures 1 row
# file: alpham_roc_rank.tsv
alpham_roc <- data.frame()
for (i in 1:6) {
    df <- read_tsv(file.path(dir, v_mode[i], "analysis", "alpham_roc_rank.tsv"))
    df$rep <- v_rep[i]
    alpham_roc <- rbind(alpham_roc, df)
}
alpham_roc <- alpham_roc %>% group_by(rank, rep, method) %>%
    summarise(FP = mean(FP, na.rm = TRUE),
              TP = mean(TP, na.rm = TRUE),
              n_ctrl = max(n_ctrl, na.rm = TRUE), 
              n_sig = max(n_sig, na.rm = TRUE))
alpham_roc$rep <- as.factor(alpham_roc$rep)
levels(alpham_roc$rep) <- v_rep_text
alpham_roc$method <-  factor(alpham_roc$method, levels = method_vec)
write_tsv(alpham_roc, file = file.path(sdir, "alpham_roc_rank.tsv"))

p2a <- ggplot(alpham_roc, aes(x = FP, y = TP)) +
  geom_path(aes(color = method), alpha = 0.7) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "rank", 
       title = data, subtitle = paste("nBenign = ", alpham_roc$n_ctrl[1], " nPathogenic = ", alpham_roc$n_sig[1], sep = "")) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(rep)) +
  geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
             data = alpham_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
cowplot::save_plot(p2a, file = file.path(sdir, "alpham_roc_rank.png"), base_width = 15, base_height = 4.5)
rm(alpham_roc)

##### plotxB: replicate level alpham_roc_rank_test # 3 figures 1 row
# file: alpham_roc_rank_test.tsv
alpham_roc <- data.frame()
for (i in 1:6) {
    df <- read_tsv(file.path(dir, v_mode[i], "analysis", "alpham_roc_rank_test.tsv"))
    df$rep <- v_rep[i]
    alpham_roc <- rbind(alpham_roc, df)
}
alpham_roc <- alpham_roc %>% group_by(rank, rep, method) %>%
    summarise(FP = mean(FP, na.rm = TRUE),
              TP = mean(TP, na.rm = TRUE),
              n_ctrl = max(n_ctrl, na.rm = TRUE), 
              n_sig = max(n_sig, na.rm = TRUE))
alpham_roc$rep <- as.factor(alpham_roc$rep)
levels(alpham_roc$rep) <- v_rep_text
alpham_roc$method <-  factor(alpham_roc$method, levels = method_vec)
write_tsv(alpham_roc, file = file.path(sdir, "alpham_roc_rank_test.tsv"))

p2b <- ggplot(alpham_roc, aes(x = FP, y = TP)) +
  geom_path(aes(color = method), alpha = 0.7) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "FP (clinvar Benign)", y = "TP (clinvar Pathogenic)", color = NULL, shape = "rank", 
       title = data, subtitle = paste("nBenign = ", alpham_roc$n_ctrl[1], " nPathogenic = ", alpham_roc$n_sig[1], sep = "")) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(rep)) +
  geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
             data = alpham_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
cowplot::save_plot(p2b, file = file.path(sdir, "alpham_roc_rank_test.png"), base_width = 15, base_height = 4.5)
rm(alpham_roc)


### MET only
if (data == "MET") {
    ##### plot3A: replicate level eve_roc_rank # 3 figures 1 row
    # file: eve_roc_rank.tsv
    eve_roc <- data.frame()
    for (i in 1:6) {
        df <- read_tsv(file.path(dir, v_mode[i], "analysis", "eve_roc_rank.tsv"))
        df$rep <- v_rep[i]
        eve_roc <- rbind(eve_roc, df)
    }
    eve_roc <- eve_roc %>% group_by(rank, rep, method) %>%
        summarise(FP = mean(FP, na.rm = TRUE),
                TP = mean(TP, na.rm = TRUE),
                n_ctrl = max(n_ctrl, na.rm = TRUE), 
                n_sig = max(n_sig, na.rm = TRUE))
    eve_roc$rep <- as.factor(eve_roc$rep)
    levels(eve_roc$rep) <- v_rep_text 
    eve_roc$method <-  factor(eve_roc$method, levels = method_vec)
    write_tsv(eve_roc, file = file.path(sdir, "eve_roc_rank.tsv"))

    p3a <- ggplot(eve_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (EVE10pct Benign)", y = "TP (EVE10pct Pathogenic)", color = NULL, shape = "rank proportion", 
            title = data, subtitle = paste("nBenign = ", eve_roc$n_ctrl[1], " nPathogenic = ", eve_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        facet_wrap(vars(rep)) +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = eve_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3a, file = file.path(sdir, "eve_roc_rank.png"), base_width = 15, base_height = 4.5)
    rm(eve_roc)

    ##### plot3B: replicate level eve_roc_rank_test # 3 figures 1 row
    # file: eve_roc_rank_test.tsv
    eve_roc <- data.frame()
    for (i in 1:6) {
        df <- read_tsv(file.path(dir, v_mode[i], "analysis", "eve_roc_rank_test.tsv"))
        df$rep <- v_rep[i]
        eve_roc <- rbind(eve_roc, df)
    }
    eve_roc <- eve_roc %>% group_by(rank, rep, method) %>%
        summarise(FP = mean(FP, na.rm = TRUE),
                TP = mean(TP, na.rm = TRUE),
                n_ctrl = max(n_ctrl, na.rm = TRUE), 
                n_sig = max(n_sig, na.rm = TRUE))
    eve_roc$rep <- as.factor(eve_roc$rep)
    levels(eve_roc$rep) <- v_rep_text 
    eve_roc$method <-  factor(eve_roc$method, levels = method_vec)
    write_tsv(eve_roc, file = file.path(sdir, "eve_roc_rank_test.tsv"))

    p3b <- ggplot(eve_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (EVE10pct Benign)", y = "TP (EVE10pct Pathogenic)", color = NULL, shape = "rank proportion", 
            title = data, subtitle = paste("nBenign = ", eve_roc$n_ctrl[1], " nPathogenic = ", eve_roc$n_sig[1], sep = "")) +
        cowplot::theme_cowplot() +
        facet_wrap(vars(rep)) +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = eve_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p3b, file = file.path(sdir, "eve_roc_rank_test.png"), base_width = 15, base_height = 4.5)
    rm(eve_roc)

}

### OCT1 only

if (data == "OCT1") {

    ##### plot4A: replicate level exp_roc_rank # 3 figures 1 row
    # file: exp_roc_rank.tsv
    exp_roc <- data.frame()
    for (i in 1:6) {
        df <- read_tsv(file.path(dir, v_mode[i], "analysis", "exp_roc_rank.tsv"))
        df$rep <- v_rep[i]
        exp_roc <- rbind(exp_roc, df)
    }
    exp_roc <- exp_roc %>% group_by(rank, rep, method) %>%
        summarise(FP = mean(FP, na.rm = TRUE),
                TP = mean(TP, na.rm = TRUE))
    exp_roc$rep <- as.factor(exp_roc$rep)
    levels(exp_roc$rep) <- v_rep_text
    exp_roc$method <-  factor(exp_roc$method, levels = method_vec)
    write_tsv(exp_roc, file = file.path(sdir, "exp_roc_rank.tsv"))

    p4a <- ggplot(exp_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (synonymous)", y = "TP (experiment)", color = NULL, shape = "rank proportion", 
            title = data) +
        cowplot::theme_cowplot() +
        facet_wrap(vars(rep)) +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = exp_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p4a, file = file.path(sdir, "exp_roc_rank.png"), base_width = 15, base_height = 4.5)
    rm(exp_roc)

    ##### plot4B: replicate level exp_roc_rank_test # 3 figures 1 row
    # file: exp_roc_rank_test.tsv
    exp_roc <- data.frame()
    for (i in 1:6) {
        df <- read_tsv(file.path(dir, v_mode[i], "analysis", "exp_roc_rank_test.tsv"))
        df$rep <- v_rep[i]
        exp_roc <- rbind(exp_roc, df)
    }
    exp_roc <- exp_roc %>% group_by(rank, rep, method) %>%
        summarise(FP = mean(FP, na.rm = TRUE),
                TP = mean(TP, na.rm = TRUE))
    exp_roc$rep <- as.factor(exp_roc$rep)
    levels(exp_roc$rep) <- v_rep_text
    exp_roc$method <-  factor(exp_roc$method, levels = method_vec)
    write_tsv(exp_roc, file = file.path(sdir, "exp_roc_rank_test.tsv"))

    p4b <- ggplot(exp_roc, aes(x = FP, y = TP)) +
        geom_path(aes(color = method), alpha = 0.7) +
        scale_color_manual(values = color_vec) +
        scale_linetype_manual(values = linetype_vec) +
        labs(x = "FP (synonymous)", y = "TP (experiment)", color = NULL, shape = "rank proportion", 
            title = data) +
        cowplot::theme_cowplot() +
        facet_wrap(vars(rep)) +
        geom_point(aes(x = FP, y = TP, color = method, shape = as.factor(rank)), 
                data = exp_roc %>% filter(rank %in% c(0.01, 0.05, 0.1, 0.2)), size = 2, alpha = 0.7) 
    cowplot::save_plot(p4b, file = file.path(sdir, "exp_roc_rank_test.png"), base_width = 15, base_height = 4.5)
    # rm(exp_roc)

    ##### plot5: replicate level power (same as figure 3B) # 1 figure
    # file: exp_power.tsv
    exp_power <- data.frame()
    for (i in 1:6) {
        df <- read_tsv(file.path(dir, v_mode[i], "analysis", "exp_power.tsv"))
        df$rep <- v_rep[i]
        exp_power <- rbind(exp_power, df)
    }
    exp_power$method <-  factor(exp_power$method, levels = method_vec)
    write_tsv(exp_power, file = file.path(sdir, "exp_roc_rank_test.tsv"))

    p5 <- ggplot(exp_power %>% filter(test_cutoff %in% c(0.05)), aes(as.factor(rep), power*10)) +
    geom_boxplot(aes(fill = method), color = "grey", alpha = 0.3) +
    geom_jitter(aes(fill = method), size = 2.5, color = "black", shape = 21,
                position = position_jitterdodge(jitter.width = 0.2)) +
    #  facet_wrap(vars(test_cutoff)) +
    cowplot::theme_cowplot() +
    scale_y_continuous(breaks = 0:10) +
    scale_fill_manual(values = color_vec) +
    scale_color_manual(values = color_vec) +
    labs(y = "# of validated variants called", x = "# of replicates")
    cowplot::save_plot(p5, file = file.path(sdir, "exp_power.png"), base_width = 6, base_height = 4)


}

