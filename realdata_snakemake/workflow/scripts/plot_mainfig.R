#!/usr/bin/env Rscript
library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")

data <- "OCT1"

##### hyperparameter #####
v_rep_text <- c("number of replicates = 1", "number of replicates = 2", "number of replicates = 3")
method_vec <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Rosace(nopos)", "Naive")
# color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#ff7f0e', '#e377c2')
# linetype_vec <- c(4, 6, 3, 4, 1, 2, 5)
color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#e377c2')
linetype_vec <- c(4, 6, 3, 4, 1, 5)

##### loading directory #####
dir <- file.path("replicates", data, "summary")
sdir <- file.path("plot_mainfig")
if (!dir.exists(sdir)){
  dir.create(sdir)
}

##### Figure 3B power plot #####
exp_power <- read_tsv(file.path(dir, "exp_roc_rank_test.tsv"))
exp_power$method <- factor(exp_power$method, levels = method_vec)
p5 <- ggplot(exp_power %>% filter(test_cutoff %in% c(0.05)) %>%
                filter(method != "Rosace(nopos)"), 
        aes(as.factor(rep), power*10)) +
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


##### Figure 3A FDR plot #####
rank_syn_fdr <- read_tsv(file = file.path(dir, "rank_syn_fdr_test.tsv"))
rank_syn_fdr$method <- factor(rank_syn_fdr$method, levels = method_vec)
rank_syn_fdr$rep <- as.factor(rank_syn_fdr$rep)
levels(rank_syn_fdr$rep) <- v_rep_text 

p1a <- ggplot(rank_syn_fdr %>% filter(method != "Rosace(nopos)") %>%
        filter(rep != "number of replicates = 2"), 
    aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  facet_wrap(vars(rep)) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL, title = data) 
cowplot::save_plot(p1a, file = file.path(sdir, "rank_fdr_syn_test.png"), base_width = 9, base_height = 4.5)


