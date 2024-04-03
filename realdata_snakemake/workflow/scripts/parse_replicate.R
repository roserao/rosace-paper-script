library("rosace")
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

##### parse argument ##### 
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein of the experiment")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data

if (!(data %in% c("OCT1", "MET"))) {
    stop("This code only applies to OCT1 and MET dataset.")
}
if (data == "OCT1") {
  key <- "1SM73"
} else if (data == "MET") {
  key <- "MET1"
}

##### directory ##### 
### load data
load(file.path("data", data, "rosace.rda"))
rosace_raw <- rosace
### save directory
dir <- file.path("replicates", data)

##### 1 replicate  ##### 
### replicate 1
sdir <- file.path(dir, "rep1_1", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
tmp <- list(rosace::Assay2AssaySet(rosace_raw@assays[[1]]))
names(tmp) <- key
rosace@assay.sets <- list()
rosace@assay.sets <- tmp
save(rosace, file = file.path(sdir, "rosace.rda"))

### replicate 2
sdir <- file.path(dir, "rep1_2", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
tmp <- list(rosace::Assay2AssaySet(rosace_raw@assays[[2]]))
names(tmp) <- key
rosace@assay.sets <- list()
rosace@assay.sets <- tmp
save(rosace, file = file.path(sdir, "rosace.rda"))

### replicate 3
sdir <- file.path(dir, "rep1_3", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
tmp <- list(rosace::Assay2AssaySet(rosace_raw@assays[[3]]))
names(tmp) <- key
rosace@assay.sets <- list()
rosace@assay.sets <- tmp
save(rosace, file = file.path(sdir, "rosace.rda"))

##### 2 replicates  ##### 
### replicate 1 & 2
sdir <- file.path(dir, "rep2_1", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
tmp <- list(IntegrateData.Assay(rosace_raw@assays[[1]], rosace_raw@assays[[2]]))
names(tmp) <- key
rosace@assay.sets <- list()
rosace@assay.sets <- tmp
save(rosace, file = file.path(sdir, "rosace.rda"))

### replicate 2 & 3
sdir <- file.path(dir, "rep2_2", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
tmp <- list(IntegrateData.Assay(rosace_raw@assays[[2]], rosace_raw@assays[[3]]))
names(tmp) <- key
rosace@assay.sets <- list()
rosace@assay.sets <- tmp
save(rosace, file = file.path(sdir, "rosace.rda"))

### replicate 1 & 3
sdir <- file.path(dir, "rep2_3", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
tmp <- list(IntegrateData.Assay(rosace_raw@assays[[1]], rosace_raw@assays[[3]]))
names(tmp) <- key
rosace@assay.sets <- list()
rosace@assay.sets <- tmp
save(rosace, file = file.path(sdir, "rosace.rda"))

##### 3 replicates  ##### 
sdir <- file.path(dir, "rep3_1", "raw")
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
}
rosace <- rosace_raw
save(rosace, file = file.path(sdir, "rosace.rda"))



