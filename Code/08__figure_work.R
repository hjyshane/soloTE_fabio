# Library load
library(Seurat)
library(ggplot2)
library(patchwork)
library(qs)
library(tidyverse)
library(SeuratExtend)
library(vioplot)


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

fabio <- readRDS("~/Slingshot_Tracy/Fabio_finalobj.rds")
solo <- qs::qread(file.path(directory[["qsave_dir"]], "03__fina_subfamily.qs"))

#### UMAP ----
umap <- SeuratExtend::DimPlot2(fabio,
                              label = T,
                              repel = T,
                              split.by = "group2",
                              theme = list(NoAxes(),
                                            NoLegend(),
                                            labs(title = NULL),
                                            theme_umap_arrows()))

#### Marker Dotplot ----
markers <- c("Pdgfra", "Reln", "Prox1", "Tle4", "Bcl11b", "Adarb2", "Grik4",
             "Lhx6", "Satb2", "Eomes", "Pax6", "Mki67")
markers  <- factor(markers, levels = markers)

Ident(fabio) <- "tracy_clusters"
dot <-
  SeuratExtend::DotPlot2(
    fabio,
    features = markers,
    group.by = "tracy_clusters",
    )

dot

#### Cell type proportion plot ----
# cell type table
Idents(fabio) <- "tracy_clusters"

table <- fabio@meta.data %>%
  as.data.frame() %>%
  dplyr::select(tracy_clusters, group2) %>%
  dplyr::group_by(tracy_clusters, group2) %>%
  dplyr::summarize(
    counts = n(),
    percentage = (counts / nrow(fabio@meta.data)) * 100)

# bar plot
bar <- ggplot(table, aes(x = group2, y = percentage, fill = tracy_clusters)) +
  geom_bar(stat = "identity",
          position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "right") +
  labs(x = "Cell Type", y = "Percentage of Cells (%)",
      title = "Cell Type Proportions")

# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "cellProportion.png"),
      plot = bar,
      width = 5, height = 8, dpi = 300)

#### DEG vandiagram ----
deg_dg <- readr::read_csv(file.path(directory[["markers_dir"]], "Tracy_deg", "DG_withM_DG_zeroM.csv"))
deg_npc <- readr::read_csv(file.path(directory[["markers_dir"]], "Tracy_deg", "NPC_withM_NPC_zeroM.csv"))
deg_nsc <- readr::read_csv(file.path(directory[["markers_dir"]], "Tracy_deg", "NSC_withM_NSC_zeroM.csv"))

library(VennDiagram)
library(ggvenn)
library(RColorBrewer)

deg_dg <- deg_dg %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(gene = ...1, direction = ifelse(avg_log2FC > 0, "up", "down"))

deg_npc <- deg_npc %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::mutate(gene = ...1, direction = ifelse(avg_log2FC > 0, "up", "down"))
deg_nsc <- deg_nsc %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::mutate(gene = ...1, direction = ifelse(avg_log2FC > 0, "up", "down"))

venn_list <- list(
  "DG with M" = deg_dg$gene,
  "NPC with M" = deg_npc$gene,
  "NSC with M" = deg_nsc$gene
)

venn_plot <- ggvenn(venn_list,
                    fill_color = brewer.pal(n = 3, name = "Set2"),
                    stroke_size = 0.5,
                    set_name_size = 4)

ggsave(filename = file.path(directory[["plot_dir"]], "venn_DEG.png"),
      plot = venn_plot,
      width = 6, height = 6, dpi = 300, bg = "white")

#### Upsetplot ----
library(ComplexUpset)

deg_list <- list(
  "DG with M" = deg_dg$gene,
  "NPC with M" = deg_npc$gene,
  "NSC with M" = deg_nsc$gene
)

direction_df <- bind_rows(
  deg_dg[, c("gene", "direction")],
  deg_npc[, c("gene", "direction")],
  deg_nsc[, c("gene", "direction")]
) %>%
  distinct(gene, .keep_all = TRUE)


# Create a data frame for the upset plot
all_genes <- unique(unlist(deg_list))
upset_data <- data.frame(gene = all_genes)
for (name in names(deg_list)) {
  upset_data[[name]] <- ifelse(upset_data$gene %in% deg_list[[name]], 1, 0)
}

upset_data <- left_join(upset_data, direction_df, by = "gene")

# upset plot
# color with up/down gene
upset_plot <- ComplexUpset::upset(
  upset_data,
  intersect = names(deg_list),
  name = "DEG Intersection",
  base_annotations = list(
    'Intersection size' = (
      ggplot(mapping = aes(x = intersection, fill = direction)) +
      geom_bar() +
      scale_fill_manual(values = c("up" = "red", "down" = "blue"))
    )
  )
)

upset_plot

ggsave(filename = file.path(directory[["plot_dir"]], "upset_DEG.png"),
      plot = upset_plot,
      width = 8, height = 6, dpi = 300, bg = "white")

#### FeaturePlot ----
Idents(solo) <- "tracy_clusters"

feature <- Seurat::FeaturePlot(solo,
                    features = c("SoloTE-L1Md-T", "SoloTE-L1Md-A", "SoloTE-L1Md-F2", "SoloTE-L1Md-F3"),
                    label = T,
                    repel = T,
                    split.by = "group2",
                    cols = c("blue", "yellow"))
feature

# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "feature_TE.png"),
      plot = feature,
      width = 16, height = 32, dpi = 300, bg = "white")



#### QC plots ----

qc <- Seurat::VlnPlot(
  solo,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  group.by = "group2",
  ncol = 2,
  pt.size = 0.1,
  alpha = 0
)

# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "QC_violin.png"),
      plot = qc,
      width = 8, height = 6, dpi = 300, bg = "white")


# qc feature
qc_scatter <- FeaturePlot(
  solo,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  split.by = "group2"
)

# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "QC_feature.png"),
      plot = qc_scatter,
      width = 12, height = 24, dpi = 300, bg = "white")


#### compare cell count solo vs fabio ----
cell_count_solo <- table(Idents(solo)) %>% as.data.frame()
cell_count_fabio <- table(Idents(fabio)) %>% as.data.frame()

cell_count_comparison <- left_join(cell_count_fabio, cell_count_solo, by = "Var1",
                                  suffix = c("_fabio", "_solo"))


#### cell cycle analysis ----
# Seurat cell cycle analysis
cc.genes.updated.2019 <- Seurat::cc.genes.updated.2019

Idents(fabio) <- "tracy_clusters"
fabio <- Seurat::CellCycleScoring(fabio,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes,
                                set.ident = TRUE)

# Feature plot of cell cycle scores
Idents(fabio) <- "tracy_clusters"
feature <- Seurat::FeaturePlot(fabio,
                    features = c("S.Score", "G2M.Score"),
                    split.by = "group2",
                    cols = c("blue", "yellow"),
                    label = T,
                    repel = T)


# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_feature.png"),
      plot = feature,
      width = 12, height = 6, dpi = 300, bg = "white")

# Violin plot of cell cycle scores
violin <- Seurat::VlnPlot(
  fabio,
  features = c("S.Score", "G2M.Score"),
  group.by = "tracy_clusters",
  split.by = "group2",
  ncol = 2,
  pt.size = 0.1,
  alpha = 0,
  split.plot = TRUE
)
# save
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_violin.png"),
      plot = violin,
      width = 8, height = 6, dpi = 300, bg = "white")

# F3 subset and deg and see
# GO term plot?? - hold
# slingshot comparison between
