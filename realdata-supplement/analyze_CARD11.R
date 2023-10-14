library(rosace)
library(tidyverse)
library(cowplot)

sdir <- "results/CARD11/analysis"
load(file = file.path(sdir, "rosace_full.RData"))
effects <- read_tsv(file.path(sdir, "effects.tsv"))

label_variant <- function(test, effect) {
  if (is.na(test) || is.na(effect)) {
    return("NA")
  } else if (test < 0.05 && effect > 0) {
    return("GOF")
  } else if (test < 0.05 && effect < 0) {
    return("LOF")
  } else {
    return("FUNC")
  }
}

effects <- effects %>% rowwise() %>%
  mutate(ROSACE_alt = label_variant(ROSACE_tests, ROSACE_effects),
         SLR_alt = label_variant(SLR_tests, SLR_effects),
         LIMMA_alt = label_variant(LIMMA_tests, LIMMA_effects),
         EDGER_alt = label_variant(EDGER_tests, EDGER_effects),
         DIMSUM_alt = label_variant(DIMSUM_tests, DIMSUM_effects),
         ENRICH2_alt = label_variant(ENRICH2_tests, ENRICH2_effects)) 


##### positive control #####
pos_ctrl <- c("p.Phe115Ile", "p.Thr117Pro", "p.Gly126Arg", "p.Phe130Cys", 
              "p.Cys49Tyr", "p.Cys49Glu", "p.Cys49Phe", "p.Cys49Asn")
df_pos_ctrl <- effects %>% filter(variants %in% pos_ctrl)

df_pos_ctrl <- data.frame(
  method = c("ROSACE", "SLR", "LIMMA", "EDGER", "DIMSUM", "ENRICH2"),
  sensitivity = colSums(df_pos_ctrl %>% select(ends_with("alt")) == "GOF") / 8)
rownames(df_pos_ctrl) <- NULL

df_pos_ctrl$method <- as.factor(df_pos_ctrl$method)
levels(df_pos_ctrl$method) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")

##### negative control #####

# df_neg_ctrl <- df_rosace %>% filter(type == "synonymous")
rank_pform <- data.frame()
for (mod in c("ENRICH2", "ROSACE", "SLR", "EDGER", "LIMMA", "DIMSUM")) {
  res <- comp_rankfdr(effects, 500, mod)
  rank_pform <- rbind(rank_pform, res)
}

# rank fdr test for each simulation
comp_rankfdr <- function(effects, resolution, model) {
  
  effects <- effects %>% mutate(ENRICH2_effects = abs(ENRICH2_effects),
                                ROSACE_effects = abs(ROSACE_effects),
                                SLR_effects = abs(SLR_effects),
                                EDGER_effects = abs(EDGER_effects),
                                LIMMA_effects = abs(LIMMA_effects),
                                DIMSUM_effects = abs(DIMSUM_effects))
  effect_ranked <- effects %>% 
    dplyr::select(variants, type,
                  starts_with(paste(model, "effects", sep = "_")),
                  starts_with(paste(model, "tests", sep = "_"))) %>%
    arrange(across(ends_with("tests")), desc(across(ends_with("effects"))))
  
  effect_syn <- effect_ranked %>% filter(type == "synonymous")
  n_syn <- sum(!is.na(effect_syn[, 3]))
  
  cutoff <- floor(nrow(effect_ranked)/resolution * (1:resolution))
  result <- data.frame(rank = 1:resolution, sense = -0.1, fdr = -0.1)
  for (i in 1:resolution) {
    effects_ranked_sig <- effect_ranked[1:cutoff[i], ]
    effects_ranked_sig <- effects_ranked_sig %>% drop_na()
    result[i, 2] <- length(intersect(effects_ranked_sig$variants, pos_ctrl))
    result[i, 3] <- sum(effects_ranked_sig$type == "synonymous", na.rm = T)/n_syn
  }
  result$method <- model
  return(result)
}
rank_pform$method <- as.factor(rank_pform$method)
levels(rank_pform$method) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")
write_tsv(rank_pform, file = file.path(sdir, "rank_pform.tsv"))

p1 <- ggplot(rank_pform, aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  cowplot::theme_cowplot() +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  scale_linetype_manual(values=c(2, 3, 1, 6, 4, 5)) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL, title = "CARD11")
save_plot(p1, file = file.path(sdir, "rank_fdr.png"), 
          base_width = 6, base_height = 4)

p2 <- ggplot(df_pos_ctrl, aes(x = method, y = sensitivity)) +
  geom_bar(aes(fill = method), stat = "identity", width = .4) +
  theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  labs(title = "CARD11", 
       y = "sensitivity (8 GOF control)")
save_plot(p2, file = file.path(sdir, "sensitivity.png"), 
          base_width = 7, base_height = 4)
  
