#devtools::load_all("../Rosace/.")
library(tidyverse)
library(cowplot)
library(rosace)


### load rosace file
load("data/plot_rosette_OCT1/rosace_1SM73.RData")
sdir <- "data/plot_rosette_OCT1/"

rosace <- runSLR(rosace, name = "1SM73_2", type = "Assay")
# hist(rosace@scores[[1]]@score$estimate, breaks = 50)
df_score <- data.frame(score = rosace@scores[[1]]@score$estimate, type = "real")
# df_score_var2 <- 
#   cbind(ExtractVarScore(rosace, name = "1SM73_2_SLR"),
#         rosace@scores[[1]]@score)
# df_score_var2 <- df_score_var2 %>% group_by(position) %>% 
#   summarise(mean = mean(estimate), sd = sd(estimate))
# df_score_var2$type <- "real"


### create rosette
load("data/plot_rosette_OCT1/rosace_1SM73.RData")
rosace <- runSLR(rosace, name = "1SM73_2", type = "Assay")
rosette <- CreateRosetteObject(object = rosace,
                               score.name = "1SM73_2_SLR",
                               pos.col = "position", mut.col = "mutation",
                               ctrl.col = "type", ctrl.name = "synonymous",
                               project.name = "1SM73_2_SLR")

### plot dispersion
p1 <- PlotDisp(object = rosace, name = "1SM73_2",
               ctrl.col = "type", ctrl.name = "synonymous")
EstimateDisp(object = rosace, name = "1SM73_2",
         ctrl.col = "type", ctrl.name = "synonymous")
p1 <- p1 + theme_cowplot() + 
  geom_hline(yintercept = 6.9, color = "red", linetype = "dashed") +
  labs(y = "estimated eta", x = "number of synonymous variants randomly chosen")
save_plot(p1, file = file.path(sdir, "eta_real.png"), base_width = 5, base_height = 4)


# PlotDisp(object = rosace, name = "1SM73_2", t = 2,
#          ctrl.col = "type", ctrl.name = "synonymous")

# step 1: generate mutant group label
hclust <- HclustMutant(rosette, save.plot = FALSE)
rosette <- GenMutLabel(rosette, hclust = hclust, Gmut = 4, save.plot = TRUE,
                       save.dir = sdir) #p2
rm(hclust)

# step 2: generate variant group label
PlotScoreHist(rosette, var.group = FALSE, mut.group = FALSE)
rosette <- GenVarLabel(rosette, Gvar = 2)
PlotScoreHist(rosette, var.group = TRUE, mut.group = TRUE)
p3 <- PlotScoreHist(rosette, var.group = TRUE, mut.group = FALSE)
p3 <- p3 + theme_cowplot() +
  scale_fill_manual(values = c("green", "blue", "red"))
save_plot(p3, file = file.path(sdir, "theta.png"), base_width = 5, base_height = 2)


# step 3: fit dirichlet
rosette <- PMVCountDist(rosette, pos.missing = 0.2)
# save(rosette, file = "tests/testdata/rosette_1SM73.RData")
#save(rosette, file = "tests/testdata/rosette_MET.RData")

# run simulation
load("data/plot_rosette_OCT1/rosette_1SM73.RData")
cfg.clean <- CreateConfig(rosette,
                          n.sim = 1, save.sim = sdir, type.sim = "growth",
                          n.rep = 3, n.round = 3,
                          null.var.group = 'var1', wt.effect = -2,
                          seq.shrink = 1.2, seq.depth = 100,
                          lib.shrink = 1,
                          var.shrink = 1, pos.flag = FALSE)
runRosette(config = cfg.clean, save.tsv = TRUE, save.rosace = TRUE, save.enrich2 = TRUE)

##### load Rosace
rosace <- readRDS(file = file.path(sdir, "growth_rep3_rd3_clean/sim1/rosace/rosace.rds"))
# key <- "simulation"
alt <- "var2"
key <- "simulation"

##### preprocessing
##### run SLR
# rosace <- FilterData(rosace, key = key, na.rmax = 0.5)
rosace <- ImputeData(rosace, key = key, impute.method = "knn", na.rmax = 0.5)
rosace <- NormalizeData(rosace,
                        key = key,
                        normalization.method = "wt",
                        wt.var.names = rosace@var.data[[1]][str_detect(rosace@var.data[[1]], "ctrl$")],
                        wt.rm = FALSE)
rosace <- IntegrateData(object = rosace, key = key)
rosace <- runSLR(rosace, name = key, type = "AssaySet")

rosace::EstimateDisp(object = rosace, name = "simulation_3",
                     ctrl.col = "mutation", ctrl.name = "ctrl")
p1b <- PlotDisp(object = rosace, name = "simulation_3",
               ctrl.col = "mutation", ctrl.name = "ctrl")
p1b <- p1b + theme_cowplot() + 
  geom_hline(yintercept = 20, color = "red", linetype = "dashed") +
  labs(y = "estimated eta", x = "number of synonymous variants randomly chosen")

save_plot(p1b, file = file.path(sdir, "eta_sim_multi.png"), base_width = 5, base_height = 4)


##### create effects table
# load ground truth
effects <- rosace@misc$effects[, 1:3] # variants expected_effects
effects <- effects %>% mutate(expected_test = ifelse(expected_group == alt, TRUE, FALSE))
# load method results
for (score in rosace@scores) {
  df_method <- score@score[, 1:3] %>% filter(variants != "_wt")
  colnames(df_method) <- c("variants",
                           paste(score@method, "effects", sep = "_"),
                           paste(score@method, "tests", sep = "_"))
  effects <- left_join(effects, df_method)
}
# effects <- effects %>% mutate(SLR_tests = p.adjust(SLR_tests, method = "fdr"))
# effects <- effects %>% mutate(SLR_alt = (SLR_tests <= 0.1),
#                               SLR.1side_alt = SLR_alt & (SLR_effects > 0))
# hist(effects$SLR_effects, breaks = 50)
# comp_fdr(effects$expected_test, effects$SLR_alt)

df_score <- rbind(df_score,
                  data.frame(score = effects$SLR_effects, type = "simulation"))
p4 <- ggplot(df_score, aes(x = score)) +
  geom_histogram(aes(y = ..density..), color = "grey", bins = 50) +
  theme_cowplot() +
  facet_wrap(vars(type), ncol = 1) +
  labs(x = "score estimated by naive method")
save_plot(p4, file = file.path(sdir, "simhist.png"), 
          base_width = 6, base_height = 5)

# simulation
df_score_var <- 
  cbind(ExtractVarScore(rosace, name = "simulation_SLR"),
        rosace@scores[[1]]@score)
df_score_var <- df_score_var %>% group_by(position) %>% 
  summarise(mean = mean(estimate), sd = sd(estimate))
df_score_var$type <- "simulation"
# real
load("data/plot_rosette/rosace_1SM73.RData")
rosace <- runSLR(rosace, name = "1SM73_2", type = "Assay")
df_score_var2 <- 
  cbind(ExtractVarScore(rosace, name = "1SM73_2_SLR"),
        rosace@scores[[1]]@score)
df_score_var2 <- df_score_var2 %>% group_by(position) %>% 
  summarise(mean = mean(estimate), sd = sd(estimate), count = n()) %>% 
  filter(count > 18)
df_score_var2$type <- "real"
# plot
df_score_var_plot <- rbind(df_score_var, df_score_var2 %>% select(-count))
p5 <- ggplot(df_score_var_plot, aes(mean, sd)) +
  geom_point(size = 0.7) +
  facet_wrap(vars(type), ncol = 2) +
  theme_cowplot() +
  labs(x = "position score mean", y = "position score standard deviation")
save_plot(p5, file = file.path(sdir, "simmeansd.png"), 
          base_width = 7, base_height = 4)
  
  
