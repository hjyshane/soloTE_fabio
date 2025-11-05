library(Seurat)
library(qs)

set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

base_dir <-  "~/SoloTE/"

sample_list <- list("M2M88", "M1M35", "M0M57", "M0M8", "F2M13", "F1M75", "F0M28", "F0M2")
csv_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "csv")
plot_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "plot")

dir.create(csv_dir)
dir.create(plot_dir)
# Process

for (i in sample_list) {
  # Directory
  input_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), paste0(sample_list[[i]][1], "_classtes_MATRIX"))
  qs_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "qsave")
  
  
  # Read object
  sobj <- qs::qread(file.path(qs_dir, paste0("02__filetered_obj_", sample_list[[i]][1], "_.qs")))
 
  # process 
  Seurat::DefaulteAssay(sobj) <- "RNA"
  
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
  qs::qsave(sobj_filtered, file.path(file.path(qs_dir, paste0("03__processed_obj_", sample_list[[i]], "_.qs"))))
  
}