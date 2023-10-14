library("tidyverse")
library("cowplot")

dir_pos <- "data/output_posxxx/growth_clean"
dir_neg <- "data/output_negxxx/growth_clean"
type <- "plot_xxx"

# dir_pos <- "data/output_pos/growth_clean"
# dir_neg <- "data/output_neg/growth_clean"
# type <- "plot_ooo"

# df_fdr_pos <- read_tsv(file.path(dir_pos, "fdr.tsv"))
# df_fdr_neg <- read_tsv(file.path(dir_neg, "fdr.tsv"))
df_rankfdr_pos <- read_tsv(file.path(dir_pos, "rankfdr.tsv"))
df_rankfdr_neg <- read_tsv(file.path(dir_neg, "rankfdr.tsv"))
df_fdrsense_pos <- read_tsv(file.path(dir_pos, "rankfdrsense.tsv"))
df_fdrsense_neg <- read_tsv(file.path(dir_neg, "rankfdrsense.tsv"))


##### FDR plot #####
# fdr_plot_pos <- df_fdr_pos %>% 
#   select(test, fdr, FDR, POWER, sim, rep, rd) %>%
#   mutate(test = sub("\\..*", "", test)) %>%
#   mutate(screen = "positive selection")
# fdr_plot_neg <- df_fdr_neg %>% 
#   select(test, fdr, FDR, POWER, sim, rep, rd) %>%
#   mutate(test = sub("\\..*", "", test)) %>%
#   mutate(screen = "negative selection")
# 
# fdr_plot <- rbind(fdr_plot_pos, fdr_plot_neg)
# fdr_plot <- fdr_plot %>%
#   filter(rd %in% c("rd1", "rd3"), rep %in% c("rep1", "rep3")) %>%
#   filter(rd != "rd1" | rep != "rep1")
# fdr_plot <- fdr_plot %>% unite("cond", c(rd, rep), sep = "_")
# # efdr = 0.05
# 
# 
# # plot 1A
# fdr_plot$test <- as.factor(fdr_plot$test)
# levels(fdr_plot$test) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")
# fdr_plot$cond <- as.factor(fdr_plot$cond)
# levels(fdr_plot$cond) <- c("R=3 T=1", "R=1 T=3", "R=3 T=3")

# ggplot(fdr_plot %>% filter(fdr == 0.05), aes(x = cond, y = FDR)) +
#   geom_boxplot(aes(color = test)) +
#   geom_jitter(aes(color = test),
#               position = position_jitterdodge(jitter.width = 0.1),
#               alpha = 0.5, size = 0.5) +
#   cowplot::theme_cowplot() +
#   facet_wrap(vars(screen), ncol = 1) +
#   labs(x = "condition", color = NULL, y = "false discovery rate") +
#   scale_color_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
#   scale_fill_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
#   geom_hline(yintercept = 0.05, linetype = 2)



##### FDR plot #####

rankfdr_plot_pos <- df_rankfdr_pos %>% 
  select(model, rank, FDR, Sensitivity, sim, rep, rd) %>%
  group_by(model, rank, rep, rd) %>%
  summarise(FDR = mean(FDR), Sensitivity = mean(Sensitivity)) %>%
  mutate(screen = "positive selection") 
rankfdr_plot_neg <- df_rankfdr_neg %>% 
  select(model, rank, FDR, Sensitivity, sim, rep, rd) %>%
  group_by(model, rank, rep, rd) %>%
  summarise(FDR = mean(FDR), Sensitivity = mean(Sensitivity)) %>%
  mutate(screen = "negative selection") 

rankfdr_plot <- rbind(rankfdr_plot_pos, rankfdr_plot_neg)
rankfdr_plot <- rankfdr_plot %>%
  filter(rd %in% c("rd1", "rd3"), rep %in% c("rep1", "rep3")) %>%
  filter(rd != "rd1" | rep != "rep1")
rankfdr_plot <- rankfdr_plot %>% unite("cond", c(rd, rep), sep = "_")

rankfdr_plot$model <- as.factor(rankfdr_plot$model)
levels(rankfdr_plot$model) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")
rankfdr_plot$cond <- as.factor(rankfdr_plot$cond)
levels(rankfdr_plot$cond) <- c("R=3 T=1", "R=1 T=3", "R=3 T=3")

p1 <- ggplot(rankfdr_plot, aes(x = rank, y = FDR)) +
  geom_step(aes(color = model, linetype = model), linewidth = 0.5, alpha = 0.8) +
  facet_grid(vars(screen), vars(cond)) +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  scale_linetype_manual(values=c(2, 3, 1, 6, 4, 5)) +
  labs(x = "ranked variants by hypothesis testing", color = NULL, y = "false discovery rate",
       color = NULL, linetype = NULL) +
  cowplot::theme_cowplot() 


##### Power plot #####

# ggplot(fdr_plot, aes(x = cond, y = POWER)) +
#   geom_boxplot(aes(color = test)) +
#   geom_jitter(aes(color = test),
#               position = position_jitterdodge(jitter.width = 0.1),
#               alpha = 0.5, size = 0.5) +
#   cowplot::theme_cowplot() +
#   facet_wrap(vars(screen), ncol = 1) +
#   labs(x = "condition", color = NULL, y = "false discovery rate") +
#   scale_fill_manual(values=c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")) +
#   scale_color_manual(values=c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")) +
#   geom_hline(yintercept = 0.05, linetype = 2)
# 

p2 <- ggplot(rankfdr_plot, aes(x = rank, y = Sensitivity)) +
  geom_step(aes(color = model, linetype = model), linewidth = 0.5, alpha = 0.8) +
  facet_grid(vars(screen), vars(cond)) +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  scale_linetype_manual(values=c(2, 3, 1, 6, 4, 5)) +
  labs(x = "ranked variants by hypothesis testing", color = NULL, y = "sensitivity",
       color = NULL, linetype = NULL) +
  cowplot::theme_cowplot()



fdrsense_plot_pos <- df_fdrsense_pos %>% 
  select(model, fdr, FDR, Power, sim, rep, rd) %>%
  group_by(model, fdr, rep, rd) %>%
  summarise(FDR = mean(FDR), Sensitivity = mean(Power)) %>%
  mutate(screen = "positive selection") 
fdrsense_plot_neg <- df_fdrsense_neg %>% 
  select(model, fdr, FDR, Power, sim, rep, rd) %>%
  group_by(model, fdr, rep, rd) %>%
  summarise(FDR = mean(FDR), Sensitivity = mean(Power)) %>%
  mutate(screen = "negative selection") 

fdrsense_plot <- rbind(fdrsense_plot_pos, fdrsense_plot_neg)
fdrsense_plot <- fdrsense_plot %>%
  filter(rd %in% c("rd1", "rd3"), rep %in% c("rep1", "rep3")) %>%
  filter(rd != "rd1" | rep != "rep1") 
fdrsense_plot <- fdrsense_plot %>% unite("cond", c(rd, rep), sep = "_")

fdrsense_plot$model <- as.factor(fdrsense_plot$model)
levels(fdrsense_plot$model) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")
fdrsense_plot$cond <- as.factor(fdrsense_plot$cond)
levels(fdrsense_plot$cond) <-  c("R=3 T=1", "R=1 T=3", "R=3 T=3")
fdrsense_plot <- fdrsense_plot %>% replace(is.na(.), 0) 
fdrsense_sig <- fdrsense_plot %>% filter(fdr %in% c(0.001, 0.01, 0.05, 0.1))


p3 <- ggplot(fdrsense_plot %>% filter(fdr != 0), aes(x = FDR, y = Sensitivity)) +
  geom_path(aes(color = model), alpha = 0.7) +
  facet_grid(vars(screen), vars(cond)) +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  scale_linetype_manual(values=c(2, 3, 1, 6, 4, 5)) +
  labs(x = "false discovery rate", y = "sensitivity", color = NULL, shape = NULL) +
  cowplot::theme_cowplot() +
  geom_point(aes(x = FDR, y = Sensitivity, color = model, shape = as.factor(fdr)), 
             data = fdrsense_sig, size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0.001, linetype = 2, alpha = 0.2) +
  geom_vline(xintercept = 0.01, linetype = 2, alpha = 0.2) +
  geom_vline(xintercept = 0.05, linetype = 2, alpha = 0.2) + 
  geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.2) 


save_plot(p1, file = file.path("data", type, "rankfdr.png"),
          base_width = 12, base_height = 6)
save_plot(p2, file = file.path("data", type, "ranksense.png"),
          base_width = 12, base_height = 6)
save_plot(p3, file = file.path("data", type, "rankfdrsense.png"),
          base_width = 12, base_height = 6)

