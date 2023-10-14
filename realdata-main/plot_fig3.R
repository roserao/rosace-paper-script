library("tidyverse")
library("rosace")

dir <- "data/OCT1/rosace/"
key_list <- c("test1", "test2", "test3",
              "test12", "test23", "test31",
              "test123")

##### load enrich2 score
# for (key in key_list) {
#   load(file.path(dir, key, "rosace.RData")) # Rosace
#   df_enrich2 <- extract_enrich_exp(key)
#   df_enrich2 <- df_enrich2 %>% dplyr::select(variant = var, estimate = enrich_score, pvalue = enrich_pval)
#   score_enrich2 <- CreateScoreObject(method = "ENRICH2",
#                                      type = "growth",
#                                      assay.name = "1SM73",
#                                      score = df_enrich2)
#   rosace <- AddScoreData(rosace, score_enrich2)
#   rm(df_enrich2, score_enrich2)
#   save(rosace, file = file.path(dir, key, "rosace.RData"))
# }
# rm(rosace, key, AddScoreData, extract_enrich_exp)

##### create effects table
for (key in key_list) {
  load(file.path(dir, key, "modv1", "rosace.RData"))
  # load method results
  effects <- rosace@var.data
  for (score in rosace@scores) {
    df_method <- score@score[, 1:3] %>% filter(variant != "_wt")
    colnames(df_method) <- c("variants", 
                             paste(score@method, "effects", sep = "_"),
                             paste(score@method, "tests", sep = "_"))
    effects <- left_join(effects, df_method)
  }
  rm(df_method, score)
  # p-value adjustment for frequentist method
  effects <- effects %>% 
    mutate(SLR_tests = p.adjust(SLR_tests, method = "fdr"),
           ROSACE_tests = ROSACE_tests * 2,
           ENRICH2_tests = p.adjust(ENRICH2_tests, method = "fdr"))
  # fdr threshold: 0.1
  effects <- effects %>% 
    mutate(SLR_alt = (SLR_tests <= 0.05),
           ENRICH2_alt = (ENRICH2_tests <= 0.05),
           ROSACE_alt = (ROSACE_tests <= 0.05),
           LIMMA_alt = (LIMMA_tests <= 0.05),
           EDGER_alt = (EDGER_tests <= 0.05))
  
  write_tsv(effects, file = file.path(dir, key, "effects.tsv"))
}
rm(rosace, effects, key)

##### analyze effects: positive selection, negative selection

# positive selection list from Willow
pos <- c("p.(E284A)", "p.(W64R)", "p.(R175H)", "p.(D149N)",
         "p.(D149R)", "p.(R486W)", "p.(E386K)", 
         "p.(D303E)", "p.(D303G)", "p.(Q241N)")


# val <- read_csv("data/OCT1/validation/oct1_scores_validation.csv")
# val <- val %>% mutate(upper = uptake_mean + 1 * uptake_sd,
#                       lower = uptake_mean - 1 * uptake_sd,
#                       test.wt = !((lower > 100) | (upper < 100)),
#                       label = ifelse(uptake_mean > 100, "GOF", "LOF"))
# pos_ctrl <- val %>% filter(!test.wt) %>% select(hgvs, label, uptake_mean, uptake_sd)
# rm(val)
# negative selection list: synonymous mutations


# organize summary table

rank_pform <- data.frame()
for (key in key_list) {
  effects <- read_tsv(file = file.path(dir, key, "effects.tsv"))
  for (mod in c("ENRICH2", "ROSACE", "SLR", "EDGER", "LIMMA")) {
    res <- comp_rankfdr(effects, 500, mod)
    res$test <- key
    rank_pform <- rbind(rank_pform, res)
  }
}
write_tsv(rank_pform, file = "data/OCT1/rank_performance.tsv")

# rank fdr test for each simulation
comp_rankfdr <- function(effects, resolution, model) {
  
  effects <- effects %>% mutate(ENRICH2_effects = abs(ENRICH2_effects),
                                ROSACE_effects = abs(ROSACE_effects),
                                SLR_effects = abs(SLR_effects),
                                EDGER_effects = abs(EDGER_effects),
                                LIMMA_effects = abs(LIMMA_effects))
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
    result[i, 2] <- length(intersect(effects_ranked_sig$variants, pos))
    result[i, 3] <- sum(effects_ranked_sig$type == "synonymous", na.rm = T)/n_syn
  }
  result$method <- model
  return(result)
}

pform_table_all <- data.frame()
for (key in key_list) {
  effects <- read_tsv(file = file.path(dir, key, "effects.tsv"))
  pform_table <- data.frame()
  for (fdr in seq(0.001, 0.1, by = 0.001)) {
    effects <- effects %>%
      mutate(SLR_alt = (SLR_tests <= fdr),
             ENRICH2_alt = (ENRICH2_tests <= fdr),
             LIMMA_alt = (LIMMA_tests <= fdr),
             EDGER_alt = (EDGER_tests <= fdr),
             ROSACE_alt = (ROSACE_tests <= fdr))
    expected_pos <- effects %>% filter(variants %in% pos)
    expected_neg <- effects %>% filter(type == "synonymous")
    row <- data.frame(key = key,
                      SLR_power = sum(expected_pos$SLR_alt, na.rm = TRUE)/10,
                      ENRICH2_power = sum(expected_pos$ENRICH2_alt, na.rm = TRUE)/10,
                      LIMMA_power = sum(expected_pos$LIMMA_alt, na.rm = TRUE)/10,
                      EDGER_power = sum(expected_pos$EDGER_alt, na.rm = TRUE)/10,
                      ROSACE_power = sum(expected_pos$ROSACE_alt, na.rm = TRUE)/10,
               #       ROSACEv2_power = sum(expected_pos$ROSACEv2_alt, na.rm = TRUE)/10,
                      SLR_fdr = sum(expected_neg$SLR_alt, na.rm = TRUE)/sum(!is.na(expected_neg$SLR_alt)),
                      ENRICH2_fdr = sum(expected_neg$ENRICH2_alt, na.rm = TRUE)/sum(!is.na(expected_neg$ENRICH2_alt)),
                      LIMMA_fdr = sum(expected_neg$LIMMA_alt, na.rm = TRUE)/sum(!is.na(expected_neg$LIMMA_alt)),
                      EDGER_fdr = sum(expected_neg$EDGER_alt, na.rm = TRUE)/sum(!is.na(expected_neg$EDGER_alt)),
                      ROSACE_fdr = sum(expected_neg$ROSACE_alt, na.rm = TRUE)/sum(!is.na(expected_neg$ROSACE_alt)),
               #       ROSACEv2_fdr = sum(expected_neg$ROSACEv2_alt, na.rm = TRUE)/sum(!is.na(expected_neg$ROSACEv2_alt)),
                      fdr = fdr)
    pform_table <- rbind(pform_table, row)
  }
  pform_table_all <- rbind(pform_table_all, pform_table)
}
rm(row, effects, expected_pos, expected_neg, key, fdr)
write_tsv(pform_table_all, file = "data/OCT1/performance.tsv")


##### plot performance table
library("tidyverse")
library("cowplot")
rank_pform <- read_tsv(file = "data/OCT1/rank_performance.tsv")

label <- data.frame(test = c("test1", "test2", "test3",
                             "test12", "test23", "test31",
                             "test123"),
                    nrep = c(rep("number of replicates = 1", 3), 
                             rep("number of replicates = 2", 3), 
                             "number of replicates = 3"))
rank_pform <- rank_pform %>% left_join(label)
rank_pform_mean <- rank_pform %>% 
  select(-test) %>% 
  group_by(nrep, rank, method) %>%
  summarise(sense = mean(sense), fdr = mean(fdr)) %>%
  ungroup()


# test 3, 31, 123
# rank fdr plot
# ggplot(rank_pform_mean, aes(x = rank, y = sense/12)) +
#   geom_line(aes(color = method)) +
#   facet_wrap(vars(nrep))

# ggplot(rank_pform_mean,
#        aes(x = rank/500, y = fdr)) +
#   geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.8) +
#   facet_wrap(vars(nrep)) +
#   cowplot::theme_cowplot() +
#   scale_color_manual(values=c("#0077BB", "#CC3311", "#009988", "orange", "purple")) +
#   scale_linetype_manual(values=c(3, 1, 6, 4, 5, 2)) +
#   labs(x = "ranked variants by hypothesis testing",
#        y = "% of synonymous called",
#        color = NULL, linetype = NULL)

rank_pform_mean$method <- as.factor(rank_pform_mean$method)
levels(rank_pform_mean$method) <- c("mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")

p1 <- ggplot(rank_pform_mean %>%
         filter(nrep != "number of replicates = 2"), 
       aes(x = rank/500, y = fdr)) +
  geom_step(aes(color = method, linetype = method), size = 0.5, alpha = 0.9) +
  facet_wrap(vars(nrep)) +
  cowplot::theme_cowplot() +
  scale_color_manual(values=c("#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  scale_linetype_manual(values=c(3, 1, 6, 4, 5)) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL)
print(p1)


# box power plot
pform_table <- read_tsv(file = "data/OCT1/performance.tsv")
df_power <- pform_table %>% 
  filter(fdr == 0.05) %>%
  select(key, ends_with("power")) 
df_power$nrep <- c(1, 1, 1, 2, 2, 2, 3)
df_power <- df_power %>%
  pivot_longer(cols = ends_with("power"),
               names_to = "method", values_to = "power") %>%
  mutate(method = sub("\\_.*", "", method))
df_power$method <- as.factor(df_power$method)
levels(df_power$method) <- c("mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")

df_power2 <- pform_table %>% 
  filter(fdr == 0.01) %>%
  select(key, ends_with("power")) 
df_power2$nrep <- c(1, 1, 1, 2, 2, 2, 3)
df_power2 <- df_power2 %>%
  pivot_longer(cols = ends_with("power"),
               names_to = "method", values_to = "power") %>%
  mutate(method = sub("\\_.*", "", method))
df_power2$method <- as.factor(df_power2$method)
levels(df_power2$method) <- c("mutscan(edgeR)", "Enrich2", "mutscan(limma)", "Rosace", "Naive")


df_power <- rbind(df_power %>% filter(method != "Rosace"),
                  df_power2 %>% filter(method == "Rosace"))

p2 <- ggplot(df_power, aes(as.factor(nrep), power*10)) +
  geom_boxplot(aes(fill = method), color = "grey", alpha = 0.3) +
  geom_jitter(aes(fill = method), size = 2.5, color = "black", shape = 21,
              position = position_jitterdodge(jitter.width = 0.2)) +
  cowplot::theme_cowplot() +
  scale_y_continuous(breaks = 0:10) +
  scale_fill_manual(values=c("#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  scale_color_manual(values=c("#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")) +
  labs(y = "# of validated variants called", x = "# of replicates")
print(p2)
 

save_plot(p1, file = file.path("data", "plot_oct1", "fdr.png"),
          base_width = 10, base_height = 4)
save_plot(p2, file = file.path("data", "plot_oct1", "power.png"),
          base_width = 6, base_height = 4)






##### function
extract_enrich_exp <- function(dir) {
  
  # scores
  file <- paste("data/OCT1/enrich2/result/", 
                dir, "/tsv/", dir, "_exp/main_identifiers_scores.tsv",
                sep = "")
  enrich_score <- read_tsv(file, skip = 2, col_names = FALSE)
  
  # change column name
  title <- read_tsv(file, n_max = 2, col_names = FALSE)
  colnames(enrich_score) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(enrich_score)[1] <- "var"
  rm(title)
  
  # clean the data frame
  enrich_score <- enrich_score %>% 
    pivot_longer(!var) %>%
    separate(name, into = c("cond", "stats")) %>%
    pivot_wider(names_from = stats, values_from = value)
  colnames(enrich_score)[3:5] <- str_c("enrich", colnames(enrich_score)[3:5], sep = "_")
  
  # pvalues 
  file <- paste("data/OCT1/enrich2/result/", 
                dir, "/tsv/", dir, "_exp/main_identifiers_scores_pvalues_wt.tsv",
                sep = "")
  enrich_pval <- read_tsv(file, skip = 2, col_names = FALSE)
  
  # change column name
  title <- read_tsv(file, n_max = 2, col_names = FALSE)
  colnames(enrich_pval) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(enrich_pval)[1] <- "var"
  rm(title)
  
  # clean the data frame
  enrich_pval <- enrich_pval %>% 
    pivot_longer(!var) %>%
    separate(name, into = c("cond", "stats"), sep = "[.]") %>%
    pivot_wider(names_from = stats, values_from = value)
  
  colnames(enrich_pval)[3:4] <- str_c("enrich", c("pval", "z"), sep = "_")
  
  # combine score and pval
  enrich_res <- full_join(enrich_pval, enrich_score)
  
  return(enrich_res)
}


#### for SookWah
rosace@scores[["1SM73_ROSACE"]]@score <-
  cbind(ExtractScore(rosace, name = "1SM73_ROSACE")@score, 
        ExtractScore(rosace, name = "1SM73_ROSACE")@optional.score %>% select(sd = `stats::sd`))
df <- OutputScore(rosace, name = "1SM73_ROSACE")
df <- df %>% 
  rowwise() %>%
  mutate(lfsr = min(pnorm(0, mean = mean, sd = sd),
                    pnorm(0, mean = mean, sd = sd, lower.tail = FALSE))) %>%
  ungroup()

df <- df %>% 
  arrange(across(lfsr), desc(across(mean)))

which(df$variants == "p.(M420del)")
which(df$variants == "p.(L160F)")
