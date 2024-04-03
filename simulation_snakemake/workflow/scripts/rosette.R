library(rosace)
library(cowplot)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

##### parse argument ##### 
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein of the experiment")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data

if (data == "OCT1") {
  key <- "1SM73"
} else if (data == "CARD11") {
  key <- "CARD11"
} else if (data == "MET") {
  key <- "MET1"
} 

##### loading directory #####
dir <- file.path("data", data)
sdir <- file.path("results", "rosette", data, "")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}

##### create rosette object #####
load(file.path(dir, "rosace.rda"))
rosace <- runSLR(rosace, name = paste(key, 1, sep = "_"), type = "Assay")
rosette <- CreateRosetteObject(object = rosace,
                               score.name = paste(key, 1, "SLR", sep = "_"),
                               pos.col = "position", mut.col = "mutation",
                               ctrl.col = "type", ctrl.name = "synonymous",
                               project.name = key)

### dispersion plot 
p1 <- PlotDisp(rosace, name = paste(key, 1, sep = "_"), 
         ctrl.col = "type", ctrl.name = "synonymous") +
  labs(title = data)
save_plot(p1, file = file.path(sdir, "disp.png"), base_width = 5, base_height = 4.5)

### hierarchical clustering                               
hclust <- HclustMutant(rosette, save.plot = TRUE, save.dir = sdir)
rosette <- GenMutLabel(rosette, hclust = hclust, Gmut = 4, save.plot = TRUE, save.dir = sdir)

### variant group clustering
rosette <- GenVarLabel(rosette, Gvar = 2)
# PlotScoreHist(rosette, var.group = FALSE, mut.group = FALSE)
# PlotScoreHist(rosette, var.group = TRUE, mut.group = TRUE)
p2 <- PlotScoreHist(rosette, var.group = TRUE, mut.group = FALSE)
p2 <- p2 + theme_cowplot() +
  scale_fill_manual(values = c("green", "blue", "red")) +
  labs(title = data)
save_plot(p2, file = file.path(sdir, "theta.png"), base_width = 5, base_height = 3)

### dirichlet 
rosette <- PMVCountDist(rosette, pos.missing = 0.2)

### save rosette object
save(rosette, file = file.path(sdir, "rosette.rda"))

##### real data effects #####
effects_real <- cbind(
  ExtractVarScore(rosace, name = paste(key, "1_SLR", sep = "_")),
  rosace@scores[[1]]@score[, -1]
)

##### generate test simulation #####
if (data == "OCT1") {
    cfg.clean <- CreateConfig(rosette,
                          n.sim = 1, save.sim = sdir, type.sim = "growth",
                          n.rep = 3, n.round = 3,
                          null.var.group = 'var1', wt.effect = -2,
                          seq.shrink = 1.2, seq.depth = 100,
                          lib.shrink = 1,
                          var.shrink = 1, pos.flag = FALSE)
} else if (data == "MET") {
    cfg.clean <- CreateConfig(rosette,
                          n.sim = 1, save.sim = sdir, type.sim = "growth",
                          n.rep = 3, n.round = 3,
                          null.var.group = 'var2', wt.effect = 2,
                          seq.shrink = 1.2, seq.depth = 100,
                          lib.shrink = 1,
                          var.shrink = 1, pos.flag = FALSE)
} else if (data == "CARD11") {
    cfg.clean <- CreateConfig(rosette,
                          n.sim = 1, save.sim = sdir, type.sim = "growth",
                          n.rep = 5, n.round = 1,
                          null.var.group = 'var2', wt.effect = 2,
                          seq.shrink = 1.2, seq.depth = 100,
                          lib.shrink = 1,
                          var.shrink = 1, pos.flag = FALSE)
}
runRosette(config = cfg.clean, save.tsv = FALSE, save.rosace = TRUE, save.enrich2 = FALSE)

##### process simulation result #####

### function
process_rosace <- function(rosace) {
  ##### preprocessing
  # rosace <- FilterData(rosace, key = key, na.rmax = 0.5)
  rosace <- ImputeData(rosace, key = key, impute.method = "knn", na.rmax = 0.5)
  rosace <- NormalizeData(rosace,
                          key = key,
                          normalization.method = "wt",
                          wt.var.names = rosace@var.data[[1]][str_detect(rosace@var.data[[1]], "ctrl$")],
                          wt.rm = FALSE)
  rosace <- IntegrateData(object = rosace, key = key)
  
  ##### run SLR
  rosace <- runSLR(rosace, name = key, type = "AssaySet")
  
  ##### create effects table
  effects <- cbind(
    ExtractVarScore(rosace, name = "simulation_SLR"),
    rosace@scores[[1]]@score[, -1]
  )
    
  return(effects)
}

### simulation data effects
key <- "simulation"
if (data %in% c("OCT1", "MET")) {
    rosace <- readRDS(file = file.path(sdir, "growth_rep3_rd3_clean/sim1/rosace/rosace.rds"))
} else if (data == "CARD11") {
    rosace <- readRDS(file = file.path(sdir, "growth_rep5_rd1_clean/sim1/rosace/rosace.rds"))
}
effects <- process_rosace(rosace)


##### plot simulation vs realdata #####
df_score <- rbind(
  effects %>% select(score = estimate, position = position) %>% mutate(type = "simulation"),
  effects_real %>% select(score = estimate, position = position) %>% mutate(type = "real data")
)

### score distribution plot
p1a <- ggplot(df_score, aes(x = score)) +
  geom_histogram(aes(y = ..density..), color = "grey", bins = 50) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(type), ncol = 2) +
  labs(x = "score estimated by naive method")
cowplot::save_plot(p1a, file = file.path(sdir, "simhist_row.png"), 
                   base_width = 12, base_height = 2.5)
p1b <- ggplot(df_score, aes(x = score)) +
  geom_histogram(aes(y = ..density..), color = "grey", bins = 50) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(type), ncol = 1) +
  labs(x = "score estimated by naive method")
cowplot::save_plot(p1b, file = file.path(sdir, "simhist_column.png"), 
                   base_width = 6, base_height = 5)

### banana curve plot 
df_score_pos <- df_score %>% group_by(position, type) %>%
  summarise(mean = mean(score), sd = sd(score))

p2 <- ggplot(df_score_pos, aes(mean, sd)) +
  geom_point(size = 0.7) +
  facet_wrap(vars(type), ncol = 2) +
  theme_cowplot() +
  labs(x = "position score mean", y = "position score standard error",
       title = data)
save_plot(p2, file = file.path(sdir, "simmeansd.png"), 
          base_width = 10, base_height = 4.5)
