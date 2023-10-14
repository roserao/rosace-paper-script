library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

dir1 <- "data/output_neg/growth_clean"
df_corr1 <- read_tsv(file.path(dir1, "corr.tsv"))
dir2 <- "data/output_negxxx/growth_clean"
df_corr2 <- read_tsv(file.path(dir2, "corr.tsv"))
dir3 <- "data/output_pos/growth_clean"
df_corr3 <- read_tsv(file.path(dir3, "corr.tsv"))
dir4 <- "data/output_posxxx/growth_clean"
df_corr4 <- read_tsv(file.path(dir4, "corr.tsv"))

df_corr <- rbind(
  df_corr1 %>% mutate(type = "position favored", sel = "negative selection"),
  df_corr2 %>% mutate(type = "position not favored", sel = "negative selection"),
  df_corr3 %>% mutate(type = "position favored", sel = "positive selection"),
  df_corr4 %>% mutate(type = "position not favored", sel = "positive selection")
)
rm(df_corr1, df_corr2, df_corr3, df_corr4,
   dir1, dir2, dir3, dir4)

corr_plot <- df_corr %>% 
  pivot_longer(cols = starts_with("corr"), names_to = c("corr_type"), values_to = "corr") %>%
  mutate(corr_type = substr(corr_type, 6, nchar(type))) 
corr_plot$test <- as.factor(corr_plot$test)
levels(corr_plot$test) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")
corr_plot <- corr_plot %>% unite("cond", c(rd, rep), sep = "_")

p <- ggplot(corr_plot, aes(x = corr_type, y = corr)) +
  geom_boxplot(aes(color = test), alpha = 1) +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  facet_grid(type + sel ~ cond)  +
  labs(y = "correlation score", x = "correlation type", color = NULL)

save_plot(p, file = file.path("data", "plot_sim", "corr.png"),
          base_width = 14, base_height = 10)
rm(p, corr_plot, df_corr)

### Time Plot

dir1 <- "data/output_neg/growth_clean"
df_time1 <- read_tsv(file.path(dir1, "benchmark.tsv"))
dir2 <- "data/output_negxxx/growth_clean"
df_time2 <- read_tsv(file.path(dir2, "benchmark.tsv"))
dir3 <- "data/output_pos/growth_clean"
df_time3 <- read_tsv(file.path(dir3, "benchmark.tsv"))
dir4 <- "data/output_posxxx/growth_clean"
df_time4 <- read_tsv(file.path(dir4, "benchmark.tsv"))

df_time <- rbind(
  df_time1 %>% mutate(type = "position favored", sel = "negative selection"),
  df_time2 %>% mutate(type = "position not favored", sel = "negative selection"),
  df_time3 %>% mutate(type = "position favored", sel = "positive selection"),
  df_time4 %>% mutate(type = "position not favored", sel = "positive selection")
)
rm(df_time1, df_time2, df_time3, df_time4,
   dir1, dir2, dir3, dir4)

time_plot <- df_time %>% 
  select(s, mean_load, rd, rep, type, sel) %>%
  pivot_longer(cols = c(s, mean_load), names_to = "stats", values_to = "value")
time_plot$test <- as.factor(time_plot$test)
time_plot <- time_plot %>% unite("cond", c(rd, rep), sep = "_")

p <- ggplot(time_plot, aes(y = value, x = cond)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, size = 0.5) +
  facet_grid(stats ~ type + sel, scales = "free") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = NULL, y = "seconds/-")
save_plot(p, file = file.path("data", "plot_sim", "time.png"),
          base_width = 13, base_height = 6)
rm(p, time_plot, df_time)
