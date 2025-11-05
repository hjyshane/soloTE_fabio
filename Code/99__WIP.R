library(Seurat)
library(ggplot2)
library(tidyverse)
# ----
set.seed(83)
# ----
base_dir <-  "~/SoloTE/F0M28/F0M28_SoloTE_output"
base_dir <- "/Volumes/hxy008-1/SoloTE/F0M28/F0M28_SoloTE_output"

directory <- list(
  mtx_dir = file.path(base_dir, "F0M28_classtes_MATRIX"),
  rds_dir = file.path(base_dir, "SeuratObject"),
  qsave_dir = file.path(base_dir, "qsave"),
  qc_dir    = file.path(base_dir, "QC"),
  csv_dir = file.path(base_dir, 'CSV')
  )

dir.create(directory[["qsave_dir"]])
dir.create(directory[["qc_dir"]])
dir.create(directory[["csv_dir"]])

# ----
expression_matrix <- ReadMtx(
  mtx = file.path(directory[["mtx_dir"]], "matrix.mtx"), 
  features = file.path(directory[["mtx_dir"]], "features.tsv"),
  cells = file.path(directory[["mtx_dir"]], "barcodes.tsv")
  )

sobj <- CreateSeuratObject(
  counts = expression_matrix, 
  project = "soloTE"
  )

# ----
sobj$orig.ident <- "FOM28"

# ----
saveRDS(sobj, file.path(directory[["rds_dir"]], "01__seurat_object.rds"))

# ----
solo <- readRDS("~/SoloTE/F0M28/F0M28_SoloTE_output/SeuratObject/01__seurat_object.rds")

# Calculate mt rna ratio
solo <- Seurat::PercentageFeatureSet(
  solo,
  pattern = '^mt-',
  col.name = 'percent_mt_rna',
  assay = "RNA") 

# Read final obj
final <- readRDS("~/Slingshot_Tracy/Fabio_finalobj.rds")
f0m28 <- subset(final, subset = sample == "F0M_28")

# get rid of _1 end 
rownames(f0m28@meta.data) <- sub("_1$", "", rownames(f0m28@meta.data))

# Filter based on f0m28 seurat object
keep <- intersect(colnames(solo), rownames(f0m28@meta.data))
solo_filtered <- subset(solo, cells = keep)

# add metadata
ord <- match(Cells(solo_filtered), rownames(f0m28@meta.data))
meta <- f0m28@meta.data[ord, , drop = FALSE]
solo_filtered@meta.data <- cbind(solo_filtered@meta.data, meta)

# save
saveRDS(solo_filtered, file.path(directory[["rds_dir"]], "02__filtered_obj.rds"))

# process
solo_process <- solo_filtered  %>% 
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

solo_process <- Seurat::NormalizeData(solo_process, assay = "RNA")

# save
saveRDS(solo_process, file.path(directory[["rds_dir"]], "02__processed_obj.rds"))


# qc
solo_process_qc <-  solo_process %>% subset(
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 200 & 
    nFeature_RNA < 7000 & 
    percent_mt_rna < 5)

# Basic report 
before <- ncol(solo_process)
after <- ncol(solo_process_qc)

hist(solo_process$percent_mt_rna)
hist(solo_process_qc$percent_mt_rna)


# save
saveRDS(solo_process_qc, file.path(directory[["rds_dir"]], "04__qc_filtered_obj.rds"))


