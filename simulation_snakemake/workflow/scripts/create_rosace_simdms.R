library("rosace")
library("readr")
library("stringr")
library("dplyr")
library("tidyr")


for (data in c("binding", "growth")) {
  key <- "simulation"
  type <- "growth"
  dir <- file.path("results", "simdms", paste(data, "simulation", "depth", 200, sep = "_"))
  truth <- read_tsv(file.path(dir, "Tsv", "expected.tsv"))

  sdir <- file.path("results", "simdms", data)
  if (!dir.exists(sdir)) {
    dir.create(sdir)
  }

  ddir <- file.path(dir, "Data")
  count_file_vec <- dir(ddir, full.names = TRUE)
  count <- data.frame()
  # file <- count_file_vec[1]
  for (file in count_file_vec) {

    rep <- str_extract(file, "rep[0-9]+")
    time <- str_extract(file, "c_[0-9]+")
    mode <- str_extract(file, "[a-z]+_rep")
    mode <- substr(mode, 1, nchar(mode) - 4)

    tmp <- read_tsv(file)
    colnames(tmp) <- c("variants", "count")
    tmp$time <- time
    tmp$rep <- rep
    tmp$mode <- mode

    count <- rbind(count, tmp)
  }


  for (m in unique(count$mode)) {
    ssdir <- file.path(sdir, m)
    if (!dir.exists(ssdir)) {
      dir.create(ssdir)
    }

    ssdir_tsv <- file.path(ssdir, "tsv")
    ssdir_rosace <- file.path(ssdir, "rosace")
    if (!dir.exists(ssdir_tsv)) {
      dir.create(ssdir_tsv)
    }
    if (!dir.exists(ssdir_rosace)) {
      dir.create(ssdir_rosace)
    }

    count_sub <- count %>% filter(mode == m)
    count_sub <- count_sub %>% pivot_wider(names_from = c("rep", "time"), values_from = count) %>% select(-mode)
    write_tsv(count_sub, file.path(ssdir_tsv, "counts.tsv"))

    assay1 <- CreateAssayObject(counts = as.matrix(count_sub %>% select(2:7)),
                                var.names = count_sub$variants,
                                key = key, rep = 1, type = type,
                                na.rm = 1, min.count = 0)
    assay2 <- CreateAssayObject(counts = as.matrix(count_sub %>% select(8:13)),
                                var.names = count_sub$variants,
                                key = key, rep = 2, type = type,
                                na.rm = 1, min.count = 0)
    assay3 <- CreateAssayObject(counts = as.matrix(count_sub %>% select(14:19)),
                                var.names = count_sub$variants,
                                key = key, rep = 3, type = type,
                                na.rm = 1, min.count = 0)
    assay4 <- CreateAssayObject(counts = as.matrix(count_sub %>% select(20:25)),
                                var.names = count_sub$variants,
                                key = key, rep = 4, type = type,
                                na.rm = 1, min.count = 0)
    assay5 <- CreateAssayObject(counts = as.matrix(count_sub %>% select(26:31)),
                                var.names = count_sub$variants,
                                key = key, rep = 5, type = type,
                                na.rm = 1, min.count = 0)

    rosace <- CreateRosaceObject(object = assay1)
    rosace <- AddAssayData(object = rosace, assay = assay2)
    rosace <- AddAssayData(object = rosace, assay = assay3)
    rosace <- AddAssayData(object = rosace, assay = assay4)
    rosace <- AddAssayData(object = rosace, assay = assay5)

    rosace <- ImputeData(rosace, key = key, impute.method = "zero")
    rosace <- NormalizeData(rosace, key = key,
                            normalization.method = "wt",
                            wt.var.names = "_wt",
                            wt.rm = TRUE)
    rosace <- IntegrateData(object = rosace, key = key)
    rosace@var.data$truth <- truth$expected_score
    saveRDS(rosace, file = file.path(ssdir_rosace, "rosace.rds"))

  }
}









