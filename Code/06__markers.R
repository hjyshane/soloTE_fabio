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
obj_sub <- qs::qread(file.path(directory[["qsave_dir"]], "03__fina_subfamily.qs"))

# set idents
Idents(obj_sub) <- "tracy_clusters"

# test
# subfamily
solote_features_sub <- rownames(obj_sub)[grep("^SoloTE-", rownames(obj_sub))]
FeaturePlot(obj_sub, features = c("SoloTE-L1Md-A", "SoloTE-L1Md-T"), split.by = "group2", label = T, repel = T, cols = c("blue", "green")) %>%
  ggplot2::ggsave(file = file.path(directory[["plot_dir"]], "feature_L1Md.png"), width = 8, height = 6, dpi = 300, bg = "white")

# markers between group
celltype <- unique(obj_sub$tracy_clusters)

Idents(obj_sub) <- "group2"

markerlist <- list()

for (ct in celltype) {
  sub <- subset(obj_sub, subset = tracy_clusters == ct)
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
readr::write_csv(markers, file = file.path(directory[["markers_dir"]], "all_genes-subfamily.csv"))

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

sigtemarkers <- sigmarkers[grep("^SoloTE",sigmarkers$gene),]

# check
sigmarkers
sigtemarkers

# save
readr::write_csv(sigmarkers, file = file.path(directory[["markers_dir"]], "sig_gene-subfamily.csv"))
readr::write_csv(sigtemarkers, file = file.path(directory[["markers_dir"]], "sig_te_gene-subfamily.csv"))
