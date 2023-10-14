library("tidyverse")
library("cowplot")

fdrsense_plot_list <- function(dir_pos, dir_neg) {
  df_fdrsense_pos <- read_tsv(file.path(dir_pos, "rankfdrsense.tsv"))
  df_fdrsense_neg <- read_tsv(file.path(dir_neg, "rankfdrsense.tsv"))
  
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
    filter(rd != "rd1" | rep != "rep1")  %>% 
    filter(rd != "rd3" | rep != "rep3")
  fdrsense_plot <- fdrsense_plot %>% unite("cond", c(rd, rep), sep = "_")
  
  fdrsense_plot$model <- as.factor(fdrsense_plot$model)
  levels(fdrsense_plot$model) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")
  fdrsense_plot$cond <- as.factor(fdrsense_plot$cond)
  levels(fdrsense_plot$cond) <-  c("R=3 T=1", "R=1 T=3")
  fdrsense_plot <- fdrsense_plot %>% replace(is.na(.), 0) 
  fdrsense_sig <- fdrsense_plot %>% filter(fdr %in% c(0.001, 0.01, 0.05, 0.1))
  rm(df_fdrsense_neg, df_fdrsense_pos, fdrsense_plot_neg, fdrsense_plot_pos)
  
  return(list(fdrsense_plot = fdrsense_plot,
              fdrsense_sig = fdrsense_sig))
}

list_xxx <- fdrsense_plot_list(dir_pos = "data/output_posxxx/growth_clean",
                               dir_neg = "data/output_negxxx/growth_clean")
list_ooo <- fdrsense_plot_list(dir_pos = "data/output_pos/growth_clean",
                               dir_neg = "data/output_neg/growth_clean")

fdrsense_plot <- rbind(list_xxx$fdrsense_plot %>% mutate(sim = "position not favored (Rosette)"),
                       list_ooo$fdrsense_plot %>% mutate(sim = "position favored (Rosette mod.)"))
fdrsense_sig <-rbind(list_xxx$fdrsense_sig %>% mutate(sim = "position not favored (Rosette)"),
                     list_ooo$fdrsense_sig %>% mutate(sim = "position favored (Rosette mod.)"))


p <- ggplot(fdrsense_plot %>% filter(fdr != 0), aes(x = FDR, y = Sensitivity)) +
  geom_path(aes(color = model), alpha = 0.7) +
  facet_grid(sim ~ screen + cond) +
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
save_plot(p, file = file.path("data", "plot_sim", "ranksense.png"),
          base_width = 14, base_height = 6)

