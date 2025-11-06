library(Seurat)
library(ggplot2)
library(patchwork)
library(qs)
library(tidyverse)
library(SeuratExtend)


# SoloTE list
# read
obj <- qs::qread(file.path(directory[["qsave_dir"]], "02__processed_solote_subfamily.qs"))

solote_features <- rownames(obj)[grep("^SoloTE-", rownames(obj))]
FeaturePlot(obj, features = c("SoloTE-LINE", "SoloTE-SINE"), split.by = "group2", label = T, repel = T, cols = c("blue", "green"))

# all markers
newmarkers <- FindAllMarkers(obj,only.pos=TRUE)
newmarkers_signif <- newmarkers[which(newmarkers$p_val_adj<=0.05),]
newmarkers_signif_te <- newmarkers_signif[grep("SoloTE",newmarkers_signif$gene),]
newmarkers_signif_te


# markers between group

celltype <- unique(obj$tracy_clusters)

Idents(obj) <- "group2"

markerlist <- list()

for (ct in celltype) {
  sub <- subset(obj, subset = tracy_clusters == ct)
  data <- FindMarkers(sub,
                      ident.1 = "zeroM",
                      ident.2 = "withM"
                      )
  
  data$cell <- ct
  
  markerlist[[ct]] <- data
}

# sig marker
sigmarkerlist <- list()
sigtelist <- list()

for (i in seq_along(markerlist)) {
  sigmarkerlist[[names(markerlist[i])]] <- markerlist[[i]][which(markerlist[[i]]$p_val_adj<=0.05),]
  sigtelist[[names(markerlist[i])]] <- sigmarkerlist[[i]][grep("SoloTE",sigmarkerlist$gene),]
  
}

