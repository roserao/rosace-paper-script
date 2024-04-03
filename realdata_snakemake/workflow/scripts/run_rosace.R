library("rosace")
library("readr")
library("stringr")
library("dplyr")
library("tidyr")

##### parse argument ##### 
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "target protein of the experiment"),
    make_option(c("-m", "--mode"), type = "character", default = "P", 
        help = "position grouping or not"),
    make_option(c("-r", "--rep"), type = "character", default = "N", 
        help = "whether analyze parsed replicate data")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
mode <- opt$mode
rep <- opt$rep

if (data == "OCT1") {
  key <- "1SM73"
} else if (data == "CARD11") {
  key <- "CARD11"
} else if (data == "MSH2") {
  key <- "MSH2"
} else if (data == "BRCA1") {
  key <- "BRCA1"
} else if (data == "MET") {
  key <- "MET1"
} else if (data == "BRCA1-RING") {
  key <- "BRCA1R"
} else if (data == "Cohesin") {
  key <- "COH"
}

### load data
if (rep == "N") {
  load(file.path("data", data, "rosace.rda"))
  dir <- file.path("results", data)
} else {
  load(file.path("replicates", data, rep, "raw", "rosace.rda"))
  dir <- file.path("replicates", data, rep)
}


if (mode == "P") {

  ##### directory ##### 
  sdir <- file.path(dir, "rosace")
  if (!dir.exists(sdir)) {
      dir.create(sdir, recursive = TRUE)
  } 

  ##### run rosace #####
  cmdstanr::set_cmdstan_path("/u/home/r/roserao/.cmdstan/cmdstan-2.34.1")
  rosace <- RunRosace(object = rosace,
                      name = key,
                      type = "AssaySet",
                      savedir = sdir, 
                      pos.col = "position",
                      ctrl.col = "type",
                      ctrl.name = "synonymous",
                      mc.cores = 4,
                      debug = FALSE,
                      install = FALSE)
  save(rosace, file = file.path(sdir, "rosace_eval.rda")) 

} else if (mode == "N") {
    ##### directory ##### 
  sdir <- file.path(dir, "rosace_nopos")
  if (!dir.exists(sdir)) {
      dir.create(sdir, recursive = TRUE)
  } 

  ##### run rosace #####
  cmdstanr::set_cmdstan_path("/u/home/r/roserao/.cmdstan/cmdstan-2.34.1")
  rosace <- RunRosace(object = rosace,
                      name = key,
                      type = "AssaySet",
                      savedir = sdir, 
                      # pos.col = "position",
                      # ctrl.col = "type",
                      # ctrl.name = "synonymous",
                      mc.cores = 4,
                      debug = FALSE,
                      install = FALSE)
  save(rosace, file = file.path(sdir, "rosace_eval.rda")) 
}






    