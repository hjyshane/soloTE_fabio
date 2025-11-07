library(Seurat)
library(ggplot2)
library(patchwork)
library(qs)
library(tidyverse)
library(SeuratExtend)


set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

base_dir <-  "~/SoloTE"

directory <- list(
  qsave_dir = file.path(base_dir, "OUTPUTS", "qsave"),
  markers_dir = file.path(base_dir, "OUTPUTS", "markers")
)

for (i in seq_along(directory)) {
  if (!dir.exists(directory[[i]])) {
    dir.create(directory[[i]])
  } else {message(directory[[i]], " exists.")}
}

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
  data$gene <- rownames(data)
  
  markerlist[[ct]] <- data
}

# combine
markers <- bind_rows(markerlist)
rownames(markers) <- NULL

# save
readr::write_csv(markers, file = file.path(directory[["markers_dir"]], "all_gene.csv"))

# sig marker
sigmarkerlist <- list()
sigtelist <- list()

for (i in seq_along(markerlist)) {
  sigmarkerlist[[names(markerlist[i])]] <- markerlist[[i]][which(markerlist[[i]]$p_val_adj<=0.05),]
}

# combine
for (i in sigmarkerlist){
  i$gene <- rownames(i)
  rownames(i) <- NULL
}

sigmarkers <- bind_rows(sigmarkerlist)
rownames(sigmarkers) <- NULL

sigtemarkers <- sigmarkers[grep("SoloTE",sigmarkerlist$gene),]

# check
sigmarkers
sigtemarkers

# save
readr::write_csv(sigmarkers, file = file.path(directory[["markers_dir"]], "sig_gene.csv"))
readr::write_csv(sigmtearkers, file = file.path(directory[["markers_dir"]], "sig_te_gene.csv"))
