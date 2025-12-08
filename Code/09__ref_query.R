library(Seurat)
library(ggplot2)
library(patchwork)
library(qs)
library(tidyverse)
library(SeuratExtend)


# set directory
base_dir <-  "~/SoloTE"

directory <- list(
  qsave_dir = file.path(base_dir, "OUTPUTS", "qsave"),
  markers_dir = file.path(base_dir, "OUTPUTS", "markers"),
  plot_dir = file.path(base_dir, "OUTPUTS", "PLOTS")
)

ref_obj <- readRDS("/home/gdbedrosianlab/hxy008/projects/PTZ_ATAC_scRNA_072024/File/Annotation_file/Jesse_mousebrain_RDS/mouseatlas_cut_processed.rds")

fabio <-  qs::qread(file.path(directory[["qsave_dir"]], "00__fabio_final_JY.qs"))

Seurat::DefaultAssay(fabio) <- "SCT"
Seurat::DefaultAssay(ref_obj) <- "SCT"

# Find anchors between sobj and ref_obj
anchors <- Seurat::FindTransferAnchors(
reference = ref_obj,
query = fabio,
dims = 1:30,
normalization.method = 'SCT')

# Predict cell types
predictions <- Seurat::TransferData(
anchorset = anchors,
refdata = ref_obj$Description,   # or whichever column contains cell type labels
dims = 1:30
)

# Add metadata
fabio <- Seurat::AddMetaData(
object = fabio,
metadata = predictions)

# Ref cell type prediction will be stored in seurat_obj$predicted.id.
Idents(fabio) <- fabio$predicted.id

# Visualize
DimPlot(fabio, reduction = "umap", label = T, repel = T) + NoLegend()

long_names <- unique(fabio$predicted.id)

short_names <- c(
  "Gran_NB_DG",
  "NPC",
  "ExN_CA3",
  "Cajal_Retzius",
  "OPC",
  "Axonic_HP",
  "ExN_CA1",
  "ExN_Cortex",
  "CGE",
  "Trilaminar_HP",
  "InN_HP",
  "InN_CCK",
  "CGE",
  "Astrocytes",
  "Gran_Neu_DG",
  "DG_Glial",
  "In_HP_CTX",
  "MGE")

# Name vector for mapping
name_map <- setNames(short_names, long_names)

# Function to replace long names to short names
shortened <- function(name) {
  ifelse(name %in% names(name_map), name_map[name], name)
}

# Apply function
fabio$mouseref.short <- shortened(fabio$predicted.id)

# Check UMAP
Idents(fabio) <- "mouseref.short"
DimPlot(fabio, reduction = "umap", label = T, repel = T) + NoLegend() +
DimPlot(fabio, reduction = "umap", label = T, repel = T, group.by = "tracy_clusters") + NoLegend()

Idents(fabio) <- "predicted.id"
DimPlot(fabio, reduction = "umap", label = T, repel = T) + NoLegend()
