library(rosace)
library(tidyverse)
library(cowplot)
library(readxl)

##### Load Data: Gene CARD11 #####
df <- read_xlsx("data/CARD11/raw/raw_Ibrutinib.xlsx", skip = 1)

key <- "CARD11"
type <- "growth"

df <- df %>% select(variants = `aa change`, ref_aa = `ref aa`, 
                    rep1_c0 = t11_c_0, rep1_c1 = t11_c_1, 
                    rep2_c0 = t12_c_0, rep2_c1 = t12_c_1, 
                    rep3_c0 = t13_c_0, rep3_c1 = t13_c_1, 
                    rep4_c0 = t14_c_0, rep4_c1 = t14_c_1, 
                    rep5_c0 = t9_c_0, rep5_c1 = t9_c_1, 
                    score_paper = log2_score, se_paper = SE_log, test_paper = score)

# only leave synonymous mutation and missense mutation
df <- df %>% rowwise() %>% 
  filter(!is.na(ref_aa)) %>%
  mutate(variants = ifelse(variants == "p.=", str_c("p.", ref_aa, substr(ref_aa, 1, 3), sep = ""), variants)) %>% 
  mutate(position = as.numeric(substr(ref_aa, 4, nchar(ref_aa))),
         wildtype = substr(variants, 3, 5),
         mutation = substr(variants, nchar(variants)-2, nchar(variants)),
         type = ifelse(mutation == wildtype, "synonymous", "missense"),
         type = ifelse(mutation == "Ter", "nonsense", type)) 
 

# test
pos_ctrl <- c("p.Phe115Ile", "p.Thr117Pro", "p.Gly126Arg", "p.Phe130Cys", 
              "p.Cys49Tyr", "p.Cys49Glu", "p.Cys49Phe", "p.Cys49Asn")
df_pos_ctrl <- df %>% filter(variants %in% pos_ctrl)
df_neg_ctrl <- df %>% filter(type == "synonymous")


##### Create Rosace Object #####

df_input <- df %>% select(1, 16:19, 3:12) %>%
  mutate(position = as.numeric(substr(variants, 6, nchar(variants)-3))) %>%
  group_by(variants, position, wildtype, mutation, type) %>% 
  summarise(rep1_c0 = sum(rep1_c0, na.rm = TRUE), rep1_c1 = sum(rep1_c1, na.rm = TRUE),
            rep2_c0 = sum(rep2_c0, na.rm = TRUE), rep2_c1 = sum(rep2_c1, na.rm = TRUE),
            rep3_c0 = sum(rep3_c0, na.rm = TRUE), rep3_c1 = sum(rep3_c1, na.rm = TRUE),
            rep4_c0 = sum(rep4_c0, na.rm = TRUE), rep4_c1 = sum(rep4_c1, na.rm = TRUE),
            rep5_c0 = sum(rep5_c0, na.rm = TRUE), rep5_c1 = sum(rep5_c1, na.rm = TRUE)) %>%
  ungroup()
df_rep1 <- df_input %>% select(variants, rep1_c0, rep1_c1) %>%
  replace(is.na(.), 0) %>%
  filter(!(rep1_c0 == 0 & rep1_c1 == 0))
df_rep2 <- df_input %>% select(variants, rep2_c0, rep2_c1) %>%
  replace(is.na(.), 0) %>%
  filter(!(rep2_c0 == 0 & rep2_c1 == 0))
df_rep3 <- df_input %>% select(variants, rep3_c0, rep3_c1) %>%
  replace(is.na(.), 0) %>%
  filter(!(rep3_c0 == 0 & rep3_c1 == 0))
df_rep4 <- df_input %>% select(variants, rep4_c0, rep4_c1) %>%
  replace(is.na(.), 0) %>%
  filter(!(rep4_c0 == 0 & rep4_c1 == 0))
df_rep5 <- df_input %>% select(variants, rep5_c0, rep5_c1) %>%
  replace(is.na(.), 0) %>%
  filter(!(rep5_c0 == 0 & rep5_c1 == 0))

assay1 <- CreateAssayObject(counts = as.matrix(df_rep1[, 2:3]),
                            var.names = df_rep1$variants,
                            key = key, rep = 1, type = type)
assay2 <- CreateAssayObject(counts = as.matrix(df_rep2[, 2:3]),
                            var.names = df_rep2$variants,
                            key = key, rep = 2, type = type)
assay3 <- CreateAssayObject(counts = as.matrix(df_rep3[, 2:3]),
                            var.names = df_rep3$variants,
                            key = key, rep = 3, type = type)
assay4 <- CreateAssayObject(counts = as.matrix(df_rep4[, 2:3]),
                            var.names = df_rep4$variants,
                            key = key, rep = 4, type = type)
assay5 <- CreateAssayObject(counts = as.matrix(df_rep5[, 2:3]),
                            var.names = df_rep5$variants,
                            key = key, rep = 5, type = type)
rosace <- CreateRosaceObject(object = assay1)
rosace <- AddAssayData(object = rosace, assay = assay2)
rosace <- AddAssayData(object = rosace, assay = assay3)
rosace <- AddAssayData(object = rosace, assay = assay4)
rosace <- AddAssayData(object = rosace, assay = assay5)

rosace@var.data <- rosace@var.data %>% left_join(df_input[, 1:5])

rosace <- NormalizeData(rosace, key = key,
                        normalization.method = "wt",
                        wt.var.names = (rosace@var.data %>% filter(type == "synonymous"))$variants,
                        wt.rm = FALSE)
rosace <- IntegrateData(object = rosace, key = key)

rosace <- runSLR(rosace, name = key, type = "AssaySet")
df_score <- OutputScore(rosace, name = paste(key, "SLR", sep = "_"))
ggplot(df_score, aes(estimate)) +
  geom_histogram(aes(fill = type), bins = 50, color = "grey") +
  facet_wrap(vars(type), ncol = 1, scale = "free_y") +
  theme_cowplot()

save(rosace, file = "data/CARD11/rosace/rosace.rda")

##### Run Rosace #####

rosace <- RunRosace(object = rosace,
                    name = "CARD11",
                    type = "AssaySet",
                    savedir = "results/CARD11/rosace/result",
                    pos.col = "position",
                    ctrl.col = "type",
                    ctrl.name = "synonymous")
save(rosace, file = "results/CARD11/rosace/result/rosace_eval.rda")
