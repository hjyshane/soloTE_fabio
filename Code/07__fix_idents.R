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

# create directory
for (i in seq_along(directory)) {
  if (!dir.exists(directory[[i]])) {
    dir.create(directory[[i]])
    message(directory[[i]], " created")
  } else {message(directory[[i]], " exists.")}
}

# SoloTE list
# read
obj_sub <- qs::qread(file.path(directory[["qsave_dir"]], "02__processed_solote_subfamily.qs"))

# fix idents in cluster
Idents(obj_sub) <- "tracy_clusters"

levels(obj_sub)

# rename cluster
obj_sub <- RenameIdents(obj_sub, 
                        "MGE_IN" = "MGE_IN", 
                        "CGE_IN" = "CGE_IN", 
                        "CA1" = "CA1", 
                        "CA3" = "CA3", 
                        "DG" = "DG", 
                        "L2-4_EN" = "L2-4_EN",
                        "L5-6_EN" = "L5-6_EN", 
                        "KI67_Prog" = "NSC", 
                        "Pax6_Prog" = "NSC",
                        "Tbr2_Prog" = "NPC", 
                        "CR" = "CR", 
                        "Neurons" = "L2-4_EN",
                        "OPCs" = "OPCs")

Idents(obj_sub) -> obj_sub@meta.data$tracy_clusters

qs::qsave(obj_sub, file.path(directory[["qsave_dir"]], "03__fina_subfamily.qs"))
