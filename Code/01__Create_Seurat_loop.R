library(Seurat)
library(qs)

set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

base_dir <-  "~/SoloTE/"

folder_list <- list("M2M88", "M1M35", "M0M57", "M0M8", "F2M13", "F1M75", "F0M28", "F0M2")

for (i in folder_list) {
  # Directory
  input_dir <- file.path(base_dir, i, paste0(i, "_SoloTE_output"), paste0(i, "_classtes_MATRIX"))
  qs_dir <- file.path(base_dir, i, paste0(i, "_SoloTE_output"), "qsave")

  dir.create(qs_dir)

  
  # Read matrix
  em <- Seurat::ReadMtx(
    mtx = file.path(input_dir, "matrix.mtx"),
    features = file.path(input_dir, "features.tsv"),
    cells = file.path(input_dir, "barcodes.tsv")
  )
  
  # Turn into Seurat object
  sobj <- Seurat::CreateSeuratObject(
    counts = em,
    project = soloTE
    )
  
  
  
  # Add sample data
  sobj$SampleID <- i
  
  # Add mt percent
  sobj <- Seurat::PercentageFeatureSet(
    sobj,
    pattern = '^mt-',
    col.name = 'percent_mt_rna',
    assay = "RNA") 
  
  # save
  qs::qsave(sobj, file = file.path(qs_dir, paste0("01__seurate_obj_", i, "_.qs")))
}