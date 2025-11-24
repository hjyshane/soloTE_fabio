library(Seurat)
library(qs)
library(tidyverse)

set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

base_dir <-  "~/SoloTE/"

sample_list <- list(M2M_88 = c("M2M88", 8),
                    M1M_35 = c("M1M35", 4),
                    M0M_57 = c("M0M57", 7),
                    M0M_8 = c("M0M8", 3),
                    F2M_13 = c("F2M13", 6),
                    F1M_75 = c("F1M75", 2),
                    F0M_28 = c("F0M28", 1),
                    F0M_2 = c("F0M2", 5))

# csv_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "csv")
# plot_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "plot")
# 
# dir.create(csv_dir)
# dir.create(plot_dir)

# Process

for (i in seq_along(sample_list)) {
  # Directory
  qs_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "qsave")
  
  
  # Read object
  sobj <- qs::qread(file.path(qs_dir, paste0("02__filtered_obj_", sample_list[[i]][1], "_locustes.qs")))
 
  # process 
  Seurat::DefaultAssay(sobj) <- "RNA"
  
  sobj <- sobj %>% 
    Seurat::SCTransform(
      method = 'glmGamPoi', 
      assay = 'RNA', 
      vars.to.regress = c('percent_mt_rna'),
      verbose = T, conserve.memory = T) %>%
    Seurat::RunPCA(
      assay = 'SCT',
      verbose = T) %>%
    Seurat::RunUMAP(
      dims = 1:30, 
      reduction = 'pca', 
      assay = 'SCT', 
      verbose = T) 
  
  sobj <- Seurat::NormalizeData(sobj, assay = "RNA")
  
  # save
  qs::qsave(sobj, file.path(file.path(qs_dir, paste0("03__processed_obj_", sample_list[[i]][1], "_locustes.qs"))))
  
}
