library(rosace)
library(tidyverse)
library(cowplot)

sdir <- "results/MSH2/analysis"
load(file = file.path(sdir, "rosace_full.RData"))
effects <- read_tsv(file.path(sdir, "effects.tsv"))

label_variant <- function(test, effect) {
  if (is.na(test) || is.na(effect)) {
    return("NA")
  } else if (test < 0.05 && effect > 0) {
    return("Deleterious")
  } else if (test < 0.05 && effect < 0) {
    return("Neutral")
  } else {
    return("Neutral")
  }
}
effects <- effects %>% rowwise() %>%
  mutate(ROSACE_alt = label_variant(ROSACE_tests, ROSACE_effects),
         SLR_alt = label_variant(SLR_tests, SLR_effects),
         LIMMA_alt = label_variant(LIMMA_tests, LIMMA_effects),
         EDGER_alt = label_variant(EDGER_tests, EDGER_effects),
         DIMSUM_alt = label_variant(DIMSUM_tests, DIMSUM_effects),
         ENRICH2_alt = label_variant(ENRICH2_tests, ENRICH2_effects)) 
effects$variants <- str_replace_all(effects$variants, "[^[:alnum:]]", "")

##### Positive and Negative Control ##### 
df_val <- read_excel("data/MSH2/raw/mmc2.xlsx", sheet = 4, skip = 1)
df_val <- df_val %>% filter(`Conflicts among previous reports?` == FALSE,
                            SpliceAI <= 0.20)
df_val <- df_val %>% 
  select(variants = Variant, truth = `Final functional classification`,
         score_paper = `LOF score`) %>%
  mutate(test_paper = ifelse(score_paper > 0, "Deleterious", "Neutral")) 
sum(df_val$truth == df_val$test_paper, na.rm = TRUE) # 154/168
df_val <- df_val %>% left_join(effects)

df_val_pos <- df_val %>% filter(truth == "Deleterious")
df_val_pos <- data.frame(
  method = c("ROSACE", "SLR", "LIMMA", "EDGER", "DIMSUM", "ENRICH2"),
  sensitivity = colSums(df_val_pos %>% select(ends_with("alt")) == df_val_pos$truth) / nrow(df_val_pos))
rownames(df_val_pos) <- NULL
df_val_pos

df_val_neg <- df_val %>% filter(truth == "Neutral")
df_val_neg <- data.frame(
  method = c("ROSACE", "SLR", "LIMMA", "EDGER", "DIMSUM", "ENRICH2"),
  specificity = colSums(df_val_neg %>% select(ends_with("alt")) == df_val_neg$truth) / nrow(df_val_neg))
rownames(df_val_neg) <- NULL
df_val_neg

df_val_result <- df_val_pos %>% left_join(df_val_neg)
df_val_result$method <- as.factor(df_val_result$method)
levels(df_val_result$method) <- c("DiMSum", "mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")

p2 <- ggplot(df_val_result, aes(x = method, y = sensitivity)) +
  geom_bar(aes(fill = method), stat = "identity", width = .4) +
  theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  labs(title = "MSH2", 
       y = "sensitivity (77 LOF control)")
save_plot(p2, file = file.path(sdir, "sensitivity.png"), 
          base_width = 7, base_height = 4)

p3 <- ggplot(df_val_result, aes(x = method, y = sensitivity)) +
  geom_bar(aes(fill = method), stat = "identity", width = .4) +
  theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  labs(title = "MSH2", 
       y = "specificity (91 Neutral control)")
save_plot(p3, file = file.path(sdir, "specificity.png"), 
          base_width = 7, base_height = 4)

##### Synonymous mutation ##### 


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
  result <- data.frame(rank = 1:resolution, fdr = -0.1)
  for (i in 1:resolution) {
    effects_ranked_sig <- effect_ranked[1:cutoff[i], ]
    effects_ranked_sig <- effects_ranked_sig %>% drop_na()
    result[i, 2] <- sum(effects_ranked_sig$type == "synonymous", na.rm = T)/n_syn
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
       color = NULL, linetype = NULL, title = "MSH2")
print(p1)
save_plot(p1, file = file.path(sdir, "rank_fdr.png"), 
          base_width = 6, base_height = 4)
