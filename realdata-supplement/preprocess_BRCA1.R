library(rosace)
library(tidyverse)
library(cowplot)
library(readxl)

##### Load Data: Gene MSH2, Cell Line HAP1 #####
df <- read_xlsx("data/BRCA1/raw/table1.xlsx", skip = 2)
key <- "BRCA1"
type <- "growth"
df <- df %>% 
  select(position_dna = `position (hg19)`, wildtype_dna = reference, mutant_dna = alt, 
         position_aa = aa_pos, wildtype_aa = aa_ref, mutant_aa = aa_alt, type = consequence,
         truth = clinvar_simple,
         score_paper = function.score.mean, test_paper = func.class, pval_paper = p.nonfunctional,
         library, d5.r1, d11.r1, d5.r2, d11.r2) %>%
  mutate(pval_paper = 1 - pval_paper)

##### label type 
df <- df %>% unite("variants", wildtype_dna, position_dna, mutant_dna, remove = FALSE)

##### Create Rosace Object #####

var.data <- df %>% 
  select(variants = variants, position_dna = position_dna, position = position_aa, type = type) %>%
  mutate(position = ifelse(position == "NA", position_dna, position),
         type = tolower(type)) 
table(var.data$position)

key <- "BRCA1"
type <- "growth"
assay1 <- CreateAssayObject(counts = as.matrix(df %>% select(library, d5.r1, d11.r1)),
                            var.names = df$variants,
                            key = key, rep = 1, type = type)
assay2 <- CreateAssayObject(counts = as.matrix(df %>% select(library, d5.r2, d11.r2)),
                            var.names = df$variants,
                            key = key, rep = 2, type = type)
rosace <- CreateRosaceObject(object = assay1)
rosace <- AddAssayData(object = rosace, assay = assay2)
GetAssayName(rosace)

rosace@var.data <- rosace@var.data %>% left_join(var.data)
#rosace <- NormalizeData(rosace, key = key, normalization.method = "total")
rosace <- NormalizeData(rosace, key = key,
                        normalization.method = "wt",
                        wt.var.names = (rosace@var.data %>% filter(type == "synonymous"))$variants,
                        wt.rm = FALSE)
rosace <- IntegrateData(object = rosace, key = key)
# rosace <- runSLR(rosace, name = key, type = "AssaySet")
# df_score <- OutputScore(rosace, name = paste(key, "SLR", sep = "_"))
# ggplot(df_score, aes(estimate)) +
#   geom_histogram(aes(fill = type), bins = 50, color = "grey") +
#   facet_wrap(vars(type), ncol = 1, scale = "free_y") +
#   theme_cowplot()

save(rosace, file = "data/BRCA1/rosace/rosace.rda")

##### Run Rosace #####

rosace <- RunRosace(object = rosace,
                    name = "BRCA1",
                    type = "AssaySet",
                    savedir = "results/BRCA1/rosace/result",
                    pos.col = "position",
                    ctrl.col = "type",
                    ctrl.name = "synonymous",
                    thred = 7)
save(rosace, file = "results/BRCA1/rosace/result/rosace_eval.rda")

