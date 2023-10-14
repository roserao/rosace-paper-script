library(rosace)
library(tidyverse)
library(cowplot)
library(readxl)

##### Load Data: Gene MSH2, Cell Line HAP1 #####
df <- read_tsv("data/MSH2/raw/GSE162130_MSH2_HAP1_raw_count.tsv")
key <- "MSH2"
type <- "growth"

##### Aggregate CODON variant into #####
df <- df %>% group_by(var_AA) %>%
  summarise(R1D_P0 = sum(R1D_P0), R1DB6_P2 = sum(R1DB6_P2), R1DB_P2 = sum(R1DB_P2),
            R2D_P0 = sum(R2D_P0), R2DB6_P2 = sum(R2DB6_P2), R2DB_P2 = sum(R2DB_P2),
            R3D_P0 = sum(R3D_P0), R3DB6_P2 = sum(R3DB6_P2), R3DB_P2 = sum(R3DB_P2))

##### Process variant data #####
var.data <- df[, 1]
label_type <- function(wt, mut) {
  if (mut == "wt") {
    return("wildtype")
  } else if (mut == "") {
    return("deletion")
  } else if (wt == mut) {
    return("synonymous")
  } else {
    return("missense")
  }
}
var.data <- var.data %>% 
  separate(var_AA, into = c("wildtype", "position", "mutation"), remove = FALSE) %>%
  rowwise()  %>% mutate(type = label_type(wildtype, mutation)) %>% ungroup()
df <- cbind(var.data, df[2:10])
rm(var.data)
df <- df %>%  # filter out control, weird behavior
  filter(position != "AA") %>%
  mutate(position = as.numeric(position),
         mutation = ifelse(mutation == "", "*", mutation)) 


##### Filter variants #####
# filter deletion
df <- df %>% filter(type != "deletion")

# time0 vs. no selection
df_nosel <- df %>% select(1:5, 6, 8, 9, 11, 12, 14)
df_nosel <- df_nosel %>% 
  filter(rowSums(df_nosel[6:11]) > 30, rowSums(df_nosel[6:11] > 0) > 3)

# time0 vs. selection
df_sel <- df %>% select(1:5, 6:7, 9:10, 12:13)
df_sel <- df_sel %>% 
  filter(rowSums(df_sel[6:11]) > 30, rowSums(df_sel[6:11] > 0) > 3)

##### Create rosace object #####
key <- "MSH2"
rosace <- create_rosace_MSH2(df_sel, key = key)
rosace_nosel <- create_rosace_MSH2(df_nosel, key = key)

save(rosace, file = "data/MSH2/rosace/rosace.rda")
save(rosace_nosel, file = "data/MSH2/rosace/rosace_nosel.rda")

create_rosace_MSH2 <- function(df, key) {
  type <- "growth"
  colnames(df)[1] <- "variants"
  
  assay1 <- CreateAssayObject(counts = as.matrix(df[, c(6:7)]),
                              var.names = df$variants,
                              key = key, rep = 1, type = type)
  assay2 <- CreateAssayObject(counts = as.matrix(df[, c(8:9)]),
                              var.names = df$variants,
                              key = key, rep = 2, type = type)
  assay3 <- CreateAssayObject(counts = as.matrix(df[, c(10:11)]),
                              var.names = df$variants,
                              key = key, rep = 3, type = type)
  rosace <- CreateRosaceObject(object = assay1)
  rosace <- AddAssayData(object = rosace, assay = assay2)
  rosace <- AddAssayData(object = rosace, assay = assay3)
  GetAssayName(rosace)
  
  rosace@var.data <- rosace@var.data %>% left_join(df[, 1:5])
  rosace <- NormalizeData(rosace, key = key,
                          normalization.method = "wt",
                          wt.var.names = (rosace@var.data %>% filter(type == "synonymous"))$variants,
                          wt.rm = FALSE)
  rosace <- IntegrateData(object = rosace, key = key)
 # rosace <- runSLR(rosace, name = key, type = "AssaySet")
  return(rosace)
}

##### Run Rosace #####

load("data/rosace_MSH2.rda")
rosace <- RunRosace(object = rosace,
                    name = "MSH2",
                    type = "AssaySet",
                    savedir = "results/MSH2/rosace", 
                    pos.col = "position", 
                    ctrl.col = "type",
                    ctrl.name = "synonymous")

