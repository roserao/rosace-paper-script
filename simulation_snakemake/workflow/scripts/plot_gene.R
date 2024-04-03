#!/usr/bin/env Rscript
library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")

# parse argument
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data

# data <- "OCT1"

##### directory #####
dir <- file.path("results", data, "summary")
sdir <- file.path("results", data, "plot")
if (!dir.exists(sdir)) {
  dir.create(sdir)
}

##### hyperparameter #####
model_vec <- c("ENRICH2", "ROSACE", "ROSACENP", "SLR", "EDGER", "LIMMA", "DIMSUM")
method_vec <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Rosace(nopos)", "Naive")
color_vec <- c('#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#d62728', '#ff7f0e', '#e377c2')
linetype_vec <- c(4, 6, 3, 4, 1, 2, 5)
pos_vec <- c("pos not favored (Rosette)", "pos favored (Rosette mod.)")

##### load file #####
df_corr <- read_tsv(file.path(dir, "corr.tsv"))
df_fdr <- read_tsv(file.path(dir, "fdr.tsv"))
df_rankfdr <- read_tsv(file.path(dir, "rankfdr.tsv"))
df_rankfdrsense <- read_tsv(file.path(dir, "rankfdrsense.tsv"))
df_time <- read_tsv(file.path(dir, "benchmark.tsv"))

##### correlation plot #####
corr_plot <- df_corr %>% 
  pivot_longer(cols = starts_with("corr"), names_to = c("corr_type"), values_to = "corr") %>%
  mutate(corr_type = substr(corr_type, 6, nchar(corr_type))) 
corr_plot$test <- as.factor(corr_plot$test)
levels(corr_plot$test) <- method_vec
corr_plot <- corr_plot %>% unite("cond", c(rd, rep), sep = "_")
corr_plot$mode <- as.factor(corr_plot$mode)
levels(corr_plot$mode) <- pos_vec

p <- ggplot(corr_plot, aes(x = corr_type, y = corr)) +
  geom_boxplot(aes(color = test), alpha = 1) +
  scale_color_manual(values = color_vec) +
  facet_grid(mode ~ cond)  +
  labs(y = "correlation score", x = "correlation type", color = NULL, title = data) +
  theme_cowplot()
save_plot(p, file = file.path(sdir, "corr.png"), base_width = 14, base_height = 6)

##### time plot #####
time_plot <- df_time %>% 
  select(s, mean_load, rd, rep, mode) %>%
  pivot_longer(cols = c(s, mean_load), names_to = "stats", values_to = "value")
time_plot <- time_plot %>% unite("cond", c(rd, rep), sep = "_")
time_plot$mode <- as.factor(time_plot$mode)
levels(time_plot$mode) <- pos_vec

p <- ggplot(time_plot, aes(y = value, x = cond)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, size = 0.5) +
  facet_grid(stats ~ mode, scales = "free") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = NULL, y = "seconds/-", title = data) +
  theme_cowplot()
save_plot(p, file = file.path(sdir, "time.png"), base_width = 6, base_height = 5)

##### fdr plot #####
fdr_plot <- df_fdr %>% 
  select(test, FDR, POWER, NAR, fdr, sim, rep, rd, mode) %>%
  pivot_longer(cols = c("FDR", "POWER"), names_to = c("statistics"), values_to = "value")
fdr_plot$test <- as.factor(fdr_plot$test)
levels(fdr_plot$test) <- method_vec
fdr_plot <- fdr_plot %>% unite("cond", c(rd, rep), sep = "_")
fdr_plot$mode <- as.factor(fdr_plot$mode)
levels(fdr_plot$mode) <- pos_vec

p <- ggplot(fdr_plot, aes(x = as.factor(fdr), y = value)) +
  geom_boxplot(aes(color = test), alpha = 1) +
  scale_color_manual(values = color_vec) +
  facet_grid(cond ~ mode + statistics)  +
  labs(y = "value", x = "fdr test threshold", color = NULL, title = data) +
  theme_cowplot()
save_plot(p, file = file.path(sdir, "fdr.png"), base_width = 14, base_height = 10)


# df <- rankfdr_plot %>% filter(rep == "rep3", rd == "rd1", mode == "pos", alt == TRUE, model == "ROSACE", Sensitivity > 0.9, rank < 0.04)
# df <- df_rankfdrsense %>% filter(rep == "rep3", rd == "rd1", mode == "pos", model == "ROSACE")
##### rankfdr plot #####
rankfdr_plot <- df_rankfdr %>% 
    group_by(model, rank, rep, rd, mode, alt) %>%
    summarise(FDR = mean(FDR), Sensitivity = mean(Sensitivity))
rankfdr_plot <- rankfdr_plot %>% unite("cond", c(rd, rep), sep = "_")
rankfdr_plot$model <- as.factor(rankfdr_plot$model)
levels(rankfdr_plot$model) <- method_vec
rankfdr_plot$mode <- as.factor(rankfdr_plot$mode)
levels(rankfdr_plot$mode) <- pos_vec
# rankfdr_plot$cond <- as.factor(rankfdr_plot$cond)
# levels(rankfdr_plot$cond) <- c("R=3 T=1", "R=1 T=3", "R=3 T=3")

p <- ggplot(rankfdr_plot %>% filter(!alt), aes(x = rank, y = FDR)) +
  geom_step(aes(color = model, linetype = model), linewidth = 0.5, alpha = 0.8) +
  facet_grid(vars(mode), vars(cond)) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "ranked variants by hypothesis testing", color = NULL, y = "false discovery rate",
       color = NULL, linetype = NULL, title = data) +
  cowplot::theme_cowplot() 
save_plot(p, file = file.path(sdir, "rankfdr.png"), base_width = 14, base_height = 6)

p <- ggplot(rankfdr_plot %>% filter(!alt), aes(x = rank, y = Sensitivity)) +
  geom_step(aes(color = model, linetype = model), linewidth = 0.5, alpha = 0.8) +
  facet_grid(vars(mode), vars(cond)) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "ranked variants by hypothesis testing", color = NULL, y = "sensitivity",
       color = NULL, linetype = NULL, title = data) +
  cowplot::theme_cowplot() 
save_plot(p, file = file.path(sdir, "ranksensitivity.png"), base_width = 14, base_height = 6)


##### fdrsense plot #####
fdrsense_plot <- df_rankfdrsense %>% 
    group_by(model, fdr, rep, rd, mode) %>%
    summarise(FDR = mean(FDR), Sensitivity = mean(Power))
fdrsense_plot <- fdrsense_plot %>% unite("cond", c(rd, rep), sep = "_")
fdrsense_plot$mode <- as.factor(fdrsense_plot$mode)
levels(fdrsense_plot$mode) <- pos_vec
fdrsense_plot$model <- as.factor(fdrsense_plot$model)
levels(fdrsense_plot$model) <- method_vec

fdrsense_plot <- fdrsense_plot %>% replace(is.na(.), 0) 
fdrsense_sig <- fdrsense_plot %>% filter(fdr %in% c(0.001, 0.01, 0.05, 0.1))

p <- ggplot(rankfdr_plot %>% filter(alt) %>% arrange(model, cond, mode, rank), aes(x = FDR, y = Sensitivity)) +
  geom_path(aes(color = model, linetype = model), alpha = 0.7) +
  facet_grid(vars(mode), vars(cond)) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "false discovery rate", y = "sensitivity", color = NULL, linetype = NULL, shape = NULL, title = data)  +
  cowplot::theme_cowplot() 
  # geom_point(aes(x = FDR, y = Sensitivity, color = model, shape = as.factor(fdr)), 
  #            data = fdrsense_plot, size = 2, alpha = 0.7) +
  # geom_vline(xintercept = 0.01, linetype = 2, alpha = 0.2) +
  # geom_vline(xintercept = 0.05, linetype = 2, alpha = 0.2) + 
  # geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.2) 
save_plot(p, file = file.path(sdir, "rankfdrsensitivity.png"), base_width = 14, base_height = 6)


p <- ggplot(fdrsense_plot %>% filter(fdr != 0), aes(x = FDR, y = Sensitivity)) +
  geom_path(aes(color = model), alpha = 0.7) +
  facet_grid(vars(mode), vars(cond)) +
  scale_color_manual(values = color_vec) +
  scale_linetype_manual(values = linetype_vec) +
  labs(x = "false discovery rate", y = "sensitivity", color = NULL, shape = NULL, title = data) +
  cowplot::theme_cowplot() +
  geom_point(aes(x = FDR, y = Sensitivity, color = model, shape = as.factor(fdr)), 
             data = fdrsense_sig, size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0.001, linetype = 2, alpha = 0.2) +
  geom_vline(xintercept = 0.01, linetype = 2, alpha = 0.2) +
  geom_vline(xintercept = 0.05, linetype = 2, alpha = 0.2) + 
  geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.2) 
save_plot(p, file = file.path(sdir, "fdrsensitivity.png"), base_width = 14, base_height = 6)