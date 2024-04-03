#!/usr/bin/env Rscript
library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")

##### hyper-parameter #####
method_vec <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Rosace(nopos)", "Naive")
# color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#ff7f0e', '#e377c2')
# linetype_vec <- c(4, 6, 3, 4, 1, 2, 5)
color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#e377c2')
linetype_vec <- c(4, 6, 3, 4, 1, 5)
pos_vec <- c("position not favored (Rosette)", "position favored (Rosette mod.)")
data_vec <- c("positive selection (OCT1)", "negative selection (MET)")

##### directory #####
oct1_dir <- file.path("results", "OCT1", "summary")
met_dir <- file.path("results", "MET", "summary")
sdir <- file.path("plot_mainfig")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}

##### rank-false discovery plot #####

func_rankfdr_plot <- function(dir) {
    df_rankfdr <- read_tsv(file.path(dir, "rankfdr.tsv"))
    rankfdr_plot <- df_rankfdr %>% 
        group_by(model, rank, rep, rd, mode, alt) %>%
        summarise(FDR = mean(FDR), Sensitivity = mean(Sensitivity))
    rankfdr_plot$model <- as.factor(rankfdr_plot$model)
    levels(rankfdr_plot$model) <- method_vec
    rankfdr_plot$mode <- as.factor(rankfdr_plot$mode)
    levels(rankfdr_plot$mode) <- pos_vec
    rankfdr_plot <- rankfdr_plot %>% unite("cond", c(rd, rep), sep = "_")
    rankfdr_plot <- rankfdr_plot %>%
        filter(alt == FALSE, cond != "rd1_rep1", model != "Rosace(nopos)")
    rankfdr_plot$cond <- as.factor(rankfdr_plot$cond)
    levels(rankfdr_plot$cond) <- c("R=3 T=1", "R=1 T=3", "R=3 T=3")

    return(rankfdr_plot)
}

rankfdr_plot <- rbind(
    func_rankfdr_plot(oct1_dir) %>% mutate(data = data_vec[1]),
    func_rankfdr_plot(met_dir) %>% mutate(data = data_vec[2])
)

p <- ggplot(rankfdr_plot %>%
    filter(mode == "position not favored (Rosette)"), aes(x = rank, y = FDR)) +
  geom_step(aes(color = model, linetype = model), linewidth = 0.5, alpha = 0.8) +
  facet_grid(vars(data), vars(cond)) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "ranked variants by hypothesis testing", color = NULL, y = "false discovery rate",
       color = NULL, linetype = NULL) +
  cowplot::theme_cowplot() 
save_plot(p, file = file.path(sdir, "rankfdr.png"), base_width = 12.5, base_height = 6)


##### false discovery rate vs sensitivity plot #####

func_fdrsense_plot <- function(dir) {
    df_rankfdrsense <- read_tsv(file.path(dir, "rankfdrsense.tsv"))
    fdrsense_plot <- df_rankfdrsense %>% 
        group_by(model, fdr, rep, rd, mode) %>%
        summarise(FDR = mean(FDR), Sensitivity = mean(Power))
    fdrsense_plot$model <- as.factor(fdrsense_plot$model)
    levels(fdrsense_plot$model) <- method_vec
    fdrsense_plot$mode <- as.factor(fdrsense_plot$mode)
    levels(fdrsense_plot$mode) <- pos_vec

    fdrsense_plot <- fdrsense_plot %>% unite("cond", c(rd, rep), sep = "_")
    fdrsense_plot <- fdrsense_plot %>%
        filter(cond != "rd1_rep1", cond != "rd3_rep3", model != "Rosace(nopos)")
    fdrsense_plot$cond <- as.factor(fdrsense_plot$cond)
    levels(fdrsense_plot$cond) <- c("R=3 T=1", "R=1 T=3")

    return(fdrsense_plot)
}

fdrsense_plot <- rbind(
    func_fdrsense_plot(oct1_dir) %>% mutate(data = data_vec[1]),
    func_fdrsense_plot(met_dir) %>% mutate(data = data_vec[2])
)
fdrsense_plot <- fdrsense_plot %>% replace(is.na(.), 0) 
fdrsense_sig <- fdrsense_plot %>% filter(fdr %in% c(0.001, 0.01, 0.05, 0.1))

p <- ggplot(fdrsense_plot %>% filter(fdr != 0), aes(x = FDR, y = Sensitivity)) +
  geom_path(aes(color = model), alpha = 0.7) +
  facet_grid(mode ~ cond + data) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "false discovery rate", y = "sensitivity", color = NULL, shape = NULL) +
  cowplot::theme_cowplot() +
  geom_point(aes(x = FDR, y = Sensitivity, color = model, shape = as.factor(fdr)), 
             data = fdrsense_sig, size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0.001, linetype = 2, alpha = 0.2) +
  geom_vline(xintercept = 0.01, linetype = 2, alpha = 0.2) +
  geom_vline(xintercept = 0.05, linetype = 2, alpha = 0.2) + 
  geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.2) 
save_plot(p, file = file.path(sdir, "fdrsensitivity.png"), base_width = 15, base_height = 7)

