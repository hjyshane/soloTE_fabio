library(Seurat)
library(ggplot2)
library(patchwork)
library(qs)
library(tidyverse)
library(sctransform)
library(SeuratExtend)

set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

base_dir <-  "~/SoloTE"

directory <- list(
  qsave_dir = file.path(base_dir, "OUTPUTS", "qsave")
)

for (i in seq_along(directory)) {
  dir.create(directory[[i]])
}

# sample list
sample_list <- list(F0M_28 = c("F0M28", 1),
                    F1M_75 = c("F1M75", 2),
                    M0M_8 = c("M0M8", 3),
                    M1M_35 = c("M1M35", 4),
                    F0M_2 = c("F0M2", 5),
                    F2M_13 = c("F2M13", 6),
                    M0M_57 = c("M0M57", 7),
                    M2M_88 = c("M2M88", 8))

# Read
obj <- qs::qread(file.path(directory[["qsave_dir"]], "01__merged_solote_subfamily.qs"))


# process
#run normalization, integration, and clustering
obj <- obj %>%
  SCTransform(conserve.memory = T) %>% 
  RunPCA(npcs = 30, verbose = F) %>% 
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4, cluster.name = "harmony_clusters") %>% 
  RunUMAP(dims = 1:30) 

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)

# join layers
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

# save
qs::qsave(obj, file.path(directory[["qsave_dir"]], "02__processed_solote_subfamily.qs"))

# read
obj <- qs::qread(file.path(directory[["qsave_dir"]], "02__processed_solote_subfamily.qs"))

# check
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "tracy_clusters"

# Examine some expected cell type markers
FeaturePlot(obj, features = c("Gad1", "Prox1", "Eomes", "Mki67", "Satb2", "Bcl11b"))
FeaturePlot(obj, features = c("Tbr1", "Tle4", "Pax6", "Vim", "Dcx", "Adarb2"))
FeaturePlot(obj, features = c("Reln", "Vip", "Sst", "P2ry12", "Gfap", "Aqp4"))
FeaturePlot(obj, features = c("Slc1a2", "Pdgfra", "Cspg4", "Olig2", "Grp", "Calb1"))
FeaturePlot(obj, features = c("Grik4", "Matn2", "Etv1", "Ndnf", "Cux2", "Foxp2"))

# make some nicer plots
VlnPlot(obj, features = c("Mki67", "Hopx", "Pax6", "Eomes", "Prox1", "Neurod1"), pt.size = 0)
VlnPlot(obj, features = c("Xist"), pt.size = 0, group.by = "sample") # female marker
VlnPlot(obj, features = c("Uty"), pt.size = 0, group.by = "sample") # male marker
VlnPlot(obj, features = c("Kdm5d"), pt.size = 0, group.by = "sample") # male marker

DotPlot(obj, features = c("Pdgfra", "Reln", "Prox1", "Tle4", "Bcl11b", "Adarb2", "Grik4", "Lhx6", "Satb2", "Eomes", "Pax6", "Mki67"))


# TE
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "tracy_clusters"
FeaturePlot(obj, features = c("SoloTE-LINE", "SoloTE-SINE"), split.by = "group2", label = T, repel = T, cols = c("blue", "green"))


