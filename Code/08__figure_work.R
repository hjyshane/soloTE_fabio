# Library load
library(Seurat)
library(ggplot2)
library(patchwork)
library(qs)
library(tidyverse)
library(SeuratExtend)

# set seed
set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

# set directory
base_dir <-  "~/SoloTE"

directory <- list(
  qsave_dir = file.path(base_dir, "OUTPUTS", "qsave"),
  markers_dir = file.path(base_dir, "OUTPUTS", "markers"),
  plot_dir = file.path(base_dir, "OUTPUTS", "PLOTS")
)

fianl <- qs::qread(file.path(directory[["qsave_dir"]], "03__fina_subfamily.qs"))

Idents(fianl) <- "tracy_clu"
FeaturePlot(fianl, features = c("SoloTE-L1Md-F3"), split.by = "group2", 
                                label = T, repel = T, cols = c("blue", "green"))            

VlnPlot(fianl, features = "SoloTE-L1Md-F2", group.by = "group2", split.by = "tracy_clusters")

# F2, F3
# F3 subset and deg and see
# UMAP split by zeroM/withM
# Dotplot with markers
# cell type proportion
# DEG vandiagram npc/nsc/dg upset with up/down expression
# check cell count between solo vs fabio
# GO term plot?? - hold
# FeaturePlot for L1 A, T, F3, F2 (dot or Vln)
# slingshot comparison between
# Cell cycle analysis btw nsc,npc,dg
# QC plot
# 