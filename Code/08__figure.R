
# %%
# Library load
library(Seurat, lib.loc = .libPaths()[2])
library(ggplot2, lib.loc = .libPaths()[1])
library(patchwork, lib.loc = .libPaths()[1])
library(qs, lib.loc = .libPaths()[2])
library(tidyverse)
library(SeuratExtend)
library(vioplot)
library(VennDiagram, lib.loc = .libPaths()[1])
library(ggvenn, lib.loc = .libPaths()[1])
library(RColorBrewer, lib.loc = .libPaths()[2])
library(ComplexUpset, lib.loc = .libPaths()[1])
library(Nebulosa, lib.loc = .libPaths()[2])
library(clusterProfiler, lib.loc = .libPaths()[2])

# %%
# set seed
set.seed(83)

# %%
options(future.globals.maxSize = 128000 * 1024^2)

# %%
# set directory
base_dir <-  "~/SoloTE"

directory <- list(
  qsave_dir = file.path(base_dir, "OUTPUTS", "qsave"),
  markers_dir = file.path(base_dir, "OUTPUTS", "markers"),
  plot_dir = file.path(base_dir, "OUTPUTS", "PLOTS")
)

# %%
fabio <-  qs::qread(file.path(directory[["qsave_dir"]], "00__fabio_final_JY.qs"))
solo <- qs::qread(file.path(directory[["qsave_dir"]], "03__final_subfamily.qs"))
sub <-  readRDS("~/Slingshot_Tracy/Fabio_sub.rds")

# factor level
fabio$group2 <- factor(
  fabio$group2,
  levels = c("zeroM", "withM")
)

solo$group2 <- factor(
  solo$group2,
  levels = c("zeroM", "withM")
)

# save
qs::qsave(fabio, file.path(directory[["qsave_dir"]], "fabio_final_JY.qs"))
qs::qsave(solo, file.path(directory[["qsave_dir"]], "03__final_subfamily.qs"))

#### UMAP ----
Idents(fabio) <- "tracy_clusters"
umap <- SeuratExtend::DimPlot2(fabio,
                              label = T,
                              repel = T,
                              split.by = "group2",
                              theme = list(NoAxes(),
                                            NoLegend(),
                                            labs(title = NULL),
                                            theme_umap_arrows()))

umap_2 <- Seurat::DimPlot(fabio, group.by = "sample")

umap
umap_2
# save
ggsave(filename = file.path(directory[["plot_dir"]], "umap_by_group.png"),
      plot = umap,
      width = 6, height = 6, dpi = 300, bg = "white")

ggsave(filename = file.path(directory[["plot_dir"]], "umap_by_sample.png"),
      plot = umap_2,
      width = 6, height = 6, dpi = 300, bg = "white")

#### Marker Dotplot ----
markers <- c("Eomes","Mki67", "Pax6", "Prox1",  "Bcl11b", "Grik4",
            "Pdgfra", "Satb2", "Tle4", "Adarb2", "Lhx6", "Reln")
markers  <- factor(markers, levels = markers)


fabio$tracy_clusters <- factor(
  fabio$tracy_clusters,
  levels = c("NPC", "NSC", "DG", "CA1", "CA3", "OPCs",
            "L2-4_EN", "L5-6_EN", "CGE_IN", "MGE_IN", "CR"))

Idents(fabio) <- "tracy_clusters"

dot <- Seurat::DotPlot(
    fabio,
    scale.min = 0,
    col.min = 0,
    features = markers,
    cols = c("darkgrey", "indianred1"),
    dot.scale = 8,
    cluster.idents = F) +
    ggplot2::theme(
        axis.text.x =
            ggplot2::element_text(
                angle = 90,
                vjust = 0.8,
                hjust = 1,
                size = 10
            ),
        axis.text.y =
            ggplot2::element_text(
                size = 10
            ),
        axis.title.x =
            ggplot2::element_blank(),
        axis.title.y =
            ggplot2::element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
    ) +
    ggplot2::coord_flip()

dot


# save
ggsave(filename = file.path(directory[["plot_dir"]], "marker_dotplot.png"),
      plot = dot,
      width = 8, height = 6, dpi = 300, bg = "white")

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
          position = "fill") +
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
deg_dg_org <- read.csv(file.path(directory[["markers_dir"]], "Tracy_deg", "DG_withM_DG_zeroM.csv"))
deg_npc_org <- read.csv(file.path(directory[["markers_dir"]], "Tracy_deg", "NPC_withM_NPC_zeroM.csv"))
deg_nsc_org <- read.csv(file.path(directory[["markers_dir"]], "Tracy_deg", "NSC_withM_NSC_zeroM.csv"))

deg_dg <- deg_dg_org %>%
  dplyr::mutate(gene = X,
                direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

deg_npc <- deg_npc_org %>%
  dplyr::mutate(gene = X,
                direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

deg_nsc <- deg_nsc_org %>%
  dplyr::mutate(gene = X,
                direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

venn_list <- list(
  "DG" = deg_dg$gene,
  "NPC" = deg_npc$gene,
  "NSC" = deg_nsc$gene
)

venn_plot <- ggvenn(venn_list,
                    fill_color = brewer.pal(n = 3, name = "Set2"),
                    stroke_size = 0.5,
                    set_name_size = 4)

ggsave(filename = file.path(directory[["plot_dir"]], "venn_DEG.png"),
      plot = venn_plot,
      width = 6, height = 6, dpi = 300, bg = "white")

#### Upsetplot ----
te_deg <- read.csv(file.path(directory[["markers_dir"]], "sig_te_gene-subfamily_NSC-NPC-DG.csv")) %>%
dplyr::mutate(direction = ifelse(avg_log2FC > 0, "up", "down"))

dg <- te_deg %>%
dplyr::filter(cell == "DG")

nsc <- te_deg %>%
dplyr::filter(cell == "NSC")

npc <- te_deg %>%
dplyr::filter(cell == "NPC")
deg_list <- list(
  "DG" = dg$gene,
  "NPC" = npc$gene,
  "NSC" = nsc$gene
)

direction_df <- bind_rows(
  dg[, c("gene", "direction")],
  npc[, c("gene", "direction")],
  nsc[, c("gene", "direction")]
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

ggsave(filename = file.path(directory[["plot_dir"]], "upset_DEG_te.png"),
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

feature_2 <- Seurat::FeaturePlot(solo,
                    features = c("SoloTE-L1Md-T", "SoloTE-L1Md-A", "SoloTE-L1Md-F2", "SoloTE-L1Md-F3"),
                    label = T,
                    repel = T,
                    cols = c("blue", "yellow"))
feature_2

# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "feature_TE.png"),
      plot = feature,
      width = 16, height = 16, dpi = 300, bg = "white")


ggsave(filename = file.path(directory[["plot_dir"]], "feature_TE_2.png"),
      plot = feature_2,
      width = 16, height = 16, dpi = 300, bg = "white")

#### Vln plot TE ---
Idents(solo) <- "tracy_clusters"
subset_te <- subset(solo, subset = tracy_clusters %in% c("NSC", "NPC", "DG"))

vln_te <- Seurat::VlnPlot(subset_te,
                    features = c("SoloTE-L1Md-F2", "SoloTE-L1Md-F3"),
                    split.by = "group2",
                    alpha = 0) +
                    theme(legend.position = "right")
vln_te



ggsave(filename = file.path(directory[["plot_dir"]], "vln_TE_subset.png"),
      plot = vln_te,
      width = 8, height = 6, dpi = 300, bg = "white")

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
cc.genes <- Seurat::cc.genes


Idents(fabio) <- "tracy_clusters"
fabio <- Seurat::CellCycleScoring(fabio,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes,
                                set.ident = TRUE)
# save
qs::qsave(fabio, file.path(directory[["qsave_dir"]], "fabio_cellCycle.qs"))

fabio <- qs::qread(file.path(directory[["qsave_dir"]], "fabio_cellCycle.qs"))
fabio <- JoinLayers(fabio)
fabio <- Seurat::NormalizeData(fabio)

qs::qsave(fabio, file.path(directory[["qsave_dir"]], "fabio_final_JY.qs"))

# Feature plot of cell cycle scores
Idents(fabio) <- "tracy_clusters"

g2m <- Seurat::FeaturePlot(fabio,
                    features = "G2M.Score",
                    split.by = "group2",
                    cols = c("blue", "yellow"),
                    label = T,
                    repel = T) & theme(legend.position = "right")

feature <- Seurat::FeaturePlot(fabio,
                    features = c("S.Score", "G2M.Score"),
                    split.by = "group2",
                    cols = c("blue", "yellow"),
                    label = T,
                    repel = T) & theme(legend.position = "right")

feature_nub <- Nebulosa::plot_density(fabio,
                                    features = c("S.Score", "G2M.Score"),
                                    reduction = "umap",
                                    pal = "inferno")

# Save plots
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_g2m.png"),
      plot = g2m,
      width = 12, height = 6, dpi = 300, bg = "white")
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_feature.png"),
      plot = feature,
      width = 12, height = 6, dpi = 300, bg = "white")

ggsave(filename = file.path(directory[["plot_dir"]], "density_feature.png"),
      plot = feature_nub,
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
) + ggplot2::theme(legend.position = "right")  # 또는 "bottom", "top", "left")

violin_g2m <- Seurat::VlnPlot(
  fabio,
  features = c("G2M.Score"),
  group.by = "tracy_clusters",
  split.by = "group2",
  pt.size = 0.01,
  alpha = 0
) +
  theme(panel.background = element_blank()) +
  geom_boxplot() +
  ggplot2::theme(legend.position = "right")  + # 또는 "bottom", "top", "left")+
  labs(x = "Cell Type", y = "G2M Phase Score", title = "G2M Phase Score by Cell Type and Group")

violin_s <- Seurat::VlnPlot(
  fabio,
  features = c("S.Score"),
  group.by = "tracy_clusters",
  split.by = "group2",
  pt.size = 0.01,
  alpha = 0
) +
  theme(panel.background = element_blank()) +
  geom_boxplot() +
  ggplot2::theme(legend.position = "right")  + # 또는 "bottom", "top", "left")+
  labs(x = "Cell Type", y = "S Phase Score", title = "S Phase Score by Cell Type and Group")


# save
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_violin.png"),
      plot = violin,
      width = 4, height = 4, dpi = 300, bg = "white")
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_violin_g2m.png"),
      plot = violin_g2m,
      width = 12, height = 6, dpi = 300, bg = "white")
ggsave(filename = file.path(directory[["plot_dir"]], "cellCycle_violin_s.png"),
      plot = violin_s,
      width = 12, height = 6, dpi = 300, bg = "white")
#### DEG ----
fabio <- PrepSCTFindMarkers(fabio)

deg_dg <- Seurat::FindMarkers(
  fabio,
  ident.1 = "withM",
  group.by = "group2",
  subset.ident = "DG",
  ident.2 = "zeroM",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
readr::write_csv(deg_dg, file.path(directory[["markers_dir"]], "deg_DG_withM_zeroM.csv"))

deg_nsc <- Seurat::FindMarkers(
  fabio,
  ident.1 = "withM",
  group.by = "group2",
  subset.ident = "NSC",
  ident.2 = "zeroM",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
readr::write_csv(deg_nsc, file.path(directory[["markers_dir"]], "deg_nsc_withM_zeroM.csv"))

deg_npc <- Seurat::FindMarkers(
  fabio,
  ident.1 = "withM",
  group.by = "group2",
  subset.ident = "NPC",
  ident.2 = "zeroM",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
readr::write_csv(deg_npc, file.path(directory[["markers_dir"]], "deg_npc_withM_zeroM.csv"))

#### GO term ----
# rerun go term analysis with new deg list
sig_dg <- deg_dg %>%
  dplyr::mutate(gene = rownames(deg_dg),
                direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

sig_npc <- deg_npc %>%
  dplyr::mutate(gene = rownames(deg_npc),
                direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

sig_nsc <- deg_nsc %>%
  dplyr::mutate(gene = rownames(deg_nsc),
                direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

rownames(deg_dg_org) <- deg_dg_org$X
rownames(deg_dg) <- deg_dg$gene

rownames(deg_nsc_org) <- deg_nsc_org$X
rownames(deg_nsc) <- deg_nsc$gene

rownames(deg_npc_org) <- deg_npc_org$X
rownames(deg_npc) <- deg_npc$gene

go_dg <- clusterProfiler::enrichGO(
  gene = rownames(deg_dg), #significant DEGs, you define what's significant based on logfc and adj pvalue
  OrgDb = 'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  universe = rownames(deg_dg_org), #the background set of genes, ie all genes expressed in your cluster. You can get this by setting logfc.threshold = 0 in your FindMarkers function so all the genes expressed are identified not just the significantly different ones.
  pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_nsc <- clusterProfiler::enrichGO(
  gene = rownames(deg_nsc), #significant DEGs, you define what's significant based on logfc and adj pvalue
  OrgDb = 'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  universe = rownames(deg_nsc_org), #the background set of genes, ie all genes expressed in your cluster. You can get this by setting logfc.threshold = 0 in your FindMarkers function so all the genes expressed are identified not just the significantly different ones.
  pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_npc <- clusterProfiler::enrichGO(
  gene = rownames(deg_npc), #significant DEGs, you define what's significant based on logfc and adj pvalue
  OrgDb = 'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  universe = rownames(deg_npc_org), #the background set of genes, ie all genes expressed in your cluster. You can get this by setting logfc.threshold = 0 in your FindMarkers function so all the genes expressed are identified not just the significantly different ones.
  pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# save
qs::qsave(go_dg, file.path(directory[["qsave_dir"]], "go_dg.qs"))
qs::qsave(go_nsc, file.path(directory[["qsave_dir"]], "go_nsc.qs"))
qs::qsave(go_npc, file.path(directory[["qsave_dir"]], "go_npc.qs"))

# go_df
dg_results <- go_dg@result
npc_results <- go_npc@result
nsc_results <- go_nsc@result

# save
readr::write_csv(dg_results, file.path(directory[["markers_dir"]], "go_dg_results.csv"))
readr::write_csv(npc_results, file.path(directory[["markers_dir"]], "go_npc_results.csv"))
readr::write_csv(nsc_results, file.path(directory[["markers_dir"]], "go_nsc_results.csv"))

# filter
dg_results <- dg_results %>%
  dplyr::filter(p.adjust < 0.05)

nsc_results <- nsc_results %>%
  dplyr::filter(p.adjust < 0.05)

npc_results <- npc_results %>%
  dplyr::filter(p.adjust < 0.05)

# save
readr::write_csv(dg_results, file.path(directory[["markers_dir"]], "go_dg_results_filtered.csv"))
readr::write_csv(nsc_results, file.path(directory[["markers_dir"]], "go_nsc_results_filtered.csv"))
readr::write_csv(npc_results, file.path(directory[["markers_dir"]], "go_npc_results_filtered.csv"))

# check clusterprofiler docs for visualization options
library(enrichplot)
barplot(go_dg, showCategory=20)
barplot(go_nsc, showCategory=20)
barplot(go_npc, showCategory=20)

# enrichplot
dg_enrich <- clusterProfiler::dotplot(go_dg, showCategory=10) +
  ggtitle("DG GO Enrichment")

nsc_enrich <- clusterProfiler::dotplot(go_nsc, showCategory=10) +
  ggtitle("NSC GO Enrichment")

npc_enrich <- clusterProfiler::dotplot(go_npc, showCategory=10) +
  ggtitle("NPC GO Enrichment")

# centplot
dg_cnet <- enrichplot::cnetplot(go_dg, showCategory = 5, foldChange = deg_dg$avg_log2FC)
nsc_cnet <- enrichplot::cnetplot(go_nsc, showCategory = 5, foldChange = deg_nsc$avg_log2FC)
npc_cnet <- enrichplot::cnetplot(go_npc, showCategory = 5, foldChange = deg_npc$avg_log2FC)

# heatplot
dg_heat <- enrichplot::heatplot(go_dg, showCategory = 5, foldChange = deg_dg$avg_log2FC)
nsc_heat <- enrichplot::heatplot(go_nsc, showCategory = 5, foldChange = deg_nsc$avg_log2)
npc_heat <- enrichplot::heatplot(go_npc, showCategory = 5, foldChange = deg_npc$avg_log2FC)

# save
ggsave(dg_enrich, filename = file.path(directory[["plot_dir"]], "go_dg_enrichment.png"), width = 8, height = 6, dpi = 300, bg = "white")
ggsave(nsc_enrich, filename = file.path(directory[["plot_dir"]], "go_nsc_enrichment.png"), width = 8, height = 6, dpi = 300, bg = "white")
ggsave(npc_enrich, filename = file.path(directory[["plot_dir"]], "go_npc_enrichment.png"), width = 8, height = 6, dpi = 300, bg = "white")

ggsave(dg_cnet, filename = file.path(directory[["plot_dir"]], "go_dg_enrichment_cnet.png"), width = 8, height = 6, dpi = 300, bg = "white")
ggsave(nsc_cnet, filename = file.path(directory[["plot_dir"]], "go_nsc_enrichment_cnet.png"), width = 8, height = 6, dpi = 300, bg = "white")
ggsave(npc_cnet, filename = file.path(directory[["plot_dir"]], "go_npc_enrichment_cnet.png"), width = 8, height = 6, dpi = 300, bg = "white")

#### feature plot solo like heatmap ----
features <- c("SoloTE-L1Md-T", "SoloTE-L1Md-A", "SoloTE-L1Md-F2", "SoloTE-L1Md-F3")

plot_list <- lapply(features, function(f) {
  p_list <- lapply(unique(solo$group2), function(g) {
    subset_obj <- subset(solo, subset = group2 == g)
    Nebulosa::plot_density(subset_obj, features = f, reduction = "umap", pal = "inferno") +
      ggtitle(paste(f, "-", g))
  })
  wrap_plots(p_list, ncol = length(unique(solo$group2)))
})

density_feature <- wrap_plots(plot_list, ncol = 1)

# save
ggsave(density_feature, filename = file.path(directory[["plot_dir"]], "density_feature_v2.png"), width = 12, height = 12, dpi = 300, bg = "white")

#### other plots to show TE difference ----
features <- c("SoloTE-L1Md-T", "SoloTE-L1Md-A", "SoloTE-L1Md-F2", "SoloTE-L1Md-F3")


#### plots for Hopx and Rbfox3 ----
# show only NSC, NPC, DG cluster
Idents(fabio) <- "tracy_clusters"
DefaultAssay(fabio) <- "RNA"
fabio <- JoinLayers(fabio)
fabio <- Seurat::NormalizeData(fabio)


vln <- Seurat::VlnPlot(fabio_sub,
                      features = c("Hopx", "Rbfox3"),
                      group.by = "tracy_clusters",
                      split.by = "group2",
                      ncol = 2,
                      pt.size = 0.1,
                      alpha = 0
                      ) & ggplot2::theme(axis.text.x = element_text(angle = 0, hjust = 1),
                      axis.title.x = element_blank())

# save
ggsave(vln, filename = file.path(directory[['plot_dir']], "vln_hopx_rbfox3.png"), width = 10, height = 8, dpi = 300, bg = "white")

#### Cell cycle plot more ----
# more cell cycle plots to show difference between withM and zeroM in NSC and NPC, DG

#### vlnplot ----
Idents(fabio) <- "tracy_clusters"
dg <- subset(fabio, subset = tracy_clusters %in% c("DG"))
nsc <- subset(fabio, subset = tracy_clusters %in% c("NSC"))
npc <- subset(fabio, subset = tracy_clusters %in% c("NPC"))

fabio_sub <- subset(fabio, subset = tracy_clusters %in% c("NSC", "NPC", "DG"))

dg$tracy_clusters <- droplevels(dg$tracy_clusters)
nsc$tracy_clusters <- droplevels(nsc$tracy_clusters)
npc$tracy_clusters <- droplevels(npc$tracy_clusters)
fabio_sub$tracy_clusters <- droplevels(fabio_sub$tracy_clusters)


nsc_vln <- Seurat::VlnPlot(nsc,
                      features = "Hopx",
                      group.by = "group2",
                      pt.size = 0.1,
                      alpha = 0
) & ggplot2::theme(axis.text.x = element_text(angle = 0, hjust = 1),
                  axis.title.x = element_blank(),
                  legend.position = "right")

nsc$group2 <- factor(nsc$group2,  levels = c("zeroM","withM"))

nsc_vln_2 <- Seurat::VlnPlot(nsc,
                      features = c("Pax6","Hes5"),
                      group.by = "group2",
                      pt.size = 0.1,
                      alpha = 0
) + ggplot2::theme(legend.position = "right") & ggplot2::theme(axis.text.x = element_text(angle = 0, hjust = 1),
                  axis.title.x = element_blank())


dg_vln <- Seurat::VlnPlot(dg,
                      features = "Rbfox3",
                      group.by = "group2",
                      pt.size = 0.1,
                      alpha = 0
) + ggplot2::theme(axis.text.x = element_text(angle = 0, hjust = 1),
                  axis.title.x = element_blank(),
                  legend.position = "right")


vln2 <- Seurat::VlnPlot(fabio_sub,
                      features = c("S.Score", "G2M.Score"),
                      group.by = "tracy_clusters",
                      split.by = "group2",
                      split.plot = TRUE,
                      ncol = 2,
                      pt.size = 0.1,
                      alpha = 0
) & ggplot2::theme(axis.text.x = element_text(angle = 0, hjust = 1),
                  axis.title.x = element_blank(),
                  legend.position = "right")

vln2

# save
ggsave(vln2, filename = file.path(directory[['plot_dir']], "vln_cellCycle_NSC_NPC_DG.png"), width = 6, height = 6, dpi = 300, bg = "white")
ggsave(dg_vln, filename = file.path(directory[['plot_dir']], "vln_cellCycle_DG_Rbfox3.png"), width = 6, height = 6, dpi = 300, bg = "white")
ggsave(nsc_vln, filename = file.path(directory[['plot_dir']], "vln_cellCycle_NSC_Hopx.png"), width = 6, height = 6, dpi = 300, bg = "white")
ggsave(nsc_vln_2, filename = file.path(directory[['plot_dir']], "vln_cellCycle_NSC_Pax6_Hse5.png"), width = 6, height = 6, dpi = 300, bg = "white")

# more plots for cell cycle in NSC, NPC, DG
fabio_sub <- subset(fabio, subset = tracy_clusters %in% c("NSC", "NPC", "DG"))

fabio_sub$tracy_clusters <- droplevels(fabio_sub$tracy_clusters)

Idents(fabio_sub) <- "group2"


feature_cc <- Seurat::FeaturePlot(fabio_sub,
                    features = c("S.Score", "G2M.Score"),
                    split.by = "tracy_clusters",
                    cols = c("blue", "yellow"),
                    label = T,
                    repel = T)

nebulosa_cc <- Nebulosa::plot_density(fabio_sub,
                                    features = c("S.Score", "G2M.Score"),
                                    reduction = "umap",
                                    pal = "inferno")

#### Heatmap with DEG by sample ----
Idents(sub) <- "tracy_clusters"
deg_dg <- read.csv(file.path(directory[["markers_dir"]], "Tracy_deg", "DG_withM_DG_zeroM.csv"))
deg_npc <- read.csv(file.path(directory[["markers_dir"]], "Tracy_deg", "NPC_withM_NPC_zeroM.csv"))
deg_nsc <- read.csv(file.path(directory[["markers_dir"]], "Tracy_deg", "NSC_withM_NSC_zeroM.csv"))

sig_dg <- deg_dg %>%
  dplyr::mutate(gene = X,
                direction = ifelse(avg_log2FC > 0, "up", "down"), X = NULL) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

sig_npc <- deg_npc %>%
  dplyr::mutate(gene = X,
                direction = ifelse(avg_log2FC > 0, "up", "down"), X = NULL) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

sig_nsc <-deg_nsc %>%
  dplyr::mutate(gene = X,
                direction = ifelse(avg_log2FC > 0, "up", "down"), X = NULL) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)


create_deg_heatmap <- function(seurat_obj, deg_df, cell_type, top_n = 50) {

  # 특정 세포 타입만 subset
  cells_subset <- seurat_obj@meta.data %>%
    dplyr::filter(tracy_clusters == cell_type) %>%
    rownames()

  # DEG 유전자 선택 (top N개 또는 전체)
  if(nrow(deg_df) > top_n) {
    selected_genes <- deg_df %>%
      dplyr::arrange(desc(abs(avg_log2FC))) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(gene)
  } else {
    selected_genes <- deg_df %>% dplyr::pull(gene)
  }

  # Expression data 추출
  exp_data <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', layer = "data") %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% selected_genes) %>%
    dplyr::select(all_of(cells_subset))

  # Metadata 준비
  meta_subset <- seurat_obj@meta.data[cells_subset, ] %>%
    dplyr::select(group2, sample) %>%
    dplyr::arrange(group2, sample)

  # 샘플별 평균 expression 계산
  sample_list <- unique(meta_subset$sample)
  avg_exp_by_sample <- sapply(sample_list, function(smp) {
    cells <- rownames(meta_subset)[meta_subset$sample == smp]
    if(length(cells) > 0) {
      rowMeans(exp_data[, cells, drop = FALSE])
    } else {
      rep(NA, nrow(exp_data))
    }
  })

  # 컬럼 순서 정렬 (group2 기준)
  sample_group <- meta_subset %>%
    dplyr::distinct(sample, group2) %>%
    dplyr::arrange(group2, sample)

  avg_exp_by_sample <- avg_exp_by_sample[, sample_group$sample]

  # Z-score normalization
  heatmap_scaled <- t(scale(t(avg_exp_by_sample)))

  # NA 제거 (모든 값이 같은 유전자)
  heatmap_scaled <- heatmap_scaled[complete.cases(heatmap_scaled), ]

  # Column annotation
  col_ha <- HeatmapAnnotation(
    Group = sample_group$group2,
    col = list(Group = c("withM" = "#E64B35", "zeroM" = "#4DBBD5")),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Group = list(title = "Group")
    )
  )

  # Row annotation (Up/Down regulation)
  gene_direction <- deg_df %>%
    dplyr::filter(gene %in% rownames(heatmap_scaled)) %>%
    dplyr::arrange(match(gene, rownames(heatmap_scaled)))

  # row_ha <- rowAnnotation(
  #   Direction = ifelse(gene_direction$avg_log2FC > 0, "Up", "Down"),
  #   col = list(Direction = c("Up" = "#E64B35", "Down" = "#4DBBD5")),
  #   annotation_legend_param = list(
  #     Direction = list(title = "Regulation")
  #   )
  # )

  # Color scale
  col_fun <- colorRamp2(c(-2, 0, 2), c("#3B4992", "white", "#EE0000"))

  # Heatmap 생성
  ht <- Heatmap(
    heatmap_scaled,
    name = "Z-score",
    col = col_fun,

    # Row settings
    cluster_rows = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    row_names_side = "left",
    # left_annotation = row_ha,

    # Column settings
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    top_annotation = col_ha,

    # Legend
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = "Z-score",
      title_position = "leftcenter-rot",
      legend_height = unit(4, "cm")
    ),

    # Title
    column_title = paste0(cell_type, " DEGs"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),

    # Size
    width = unit(10, "cm"),
    height = unit(14, "cm")
  )

  return(ht)
}

# 각 세포 타입별 Heatmap 생성
ht_dg <- create_deg_heatmap(sub, sig_dg, "DG", top_n = 1000)
ht_npc <- create_deg_heatmap(sub, sig_npc, "NPC", top_n = 1000)
ht_nsc <- create_deg_heatmap(sub, sig_nsc, "NSC", top_n = 1000)

# 개별 출력
png(file.path(directory[["plot_dir"]], "DG_DEG_heatmap_all.png"), width = 12, height = 16, units = "in", res = 300, bg = "white")
draw(ht_dg)
dev.off()

png(file.path(directory[["plot_dir"]], "NPC_DEG_heatmap_all.png"), width = 12, height = 16, units = "in", res = 300, bg = "white")
draw(ht_npc)
dev.off()

png(file.path(directory[["plot_dir"]], "NSC_DEG_heatmap_all.png"), width = 12, height = 16, units = "in", res = 300, bg = "white")
draw(ht_nsc)
dev.off()


#### Subset SoloTE-L1Md-F3 cells ----
# SoloTE-L1Md-F3 subset and deg and see
# cells expressing SoloTE-L1Md-F3
Idents(fabio) <- "tracy_clusters"
cells_L1MdF3 <- WhichCells(solo, expression = `SoloTE-L1Md-F3` > 0)

fabio_L1MdF3 <- subset(fabio, cells = cells_L1MdF3)
fabio_notL1MdF3 <- subset(fabio, cells = setdiff(Cells(fabio), cells_L1MdF3))


#### cell ratio for gene expression ----
# calculate the ratio of cells expressing a gene in each group
calculate_cell_ratio <- function(seurat_obj, gene) {
  # Get the expression matrix for the gene
  gene_expr <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')[gene, ]
  # Get the metadata for the group
  group_meta <- seurat_obj$group2
  # Calculate the ratio of cells expressing the gene in each group
  cell_ratio <- table(group_meta[gene_expr > 0]) / table(group_meta)
  return(cell_ratio)
}

# Calculate the cell ratio for each gene
genes_of_interest <- c("Mki67", "Hopx", "Eomes", "Sox2", "Gfap")

cell_ratios <- lapply(genes_of_interest, function(gene) {
  calculate_cell_ratio(fabio, gene)
})

# Co-express or not co-express
co_express <- function(seurat_obj, gene1, gene2) {
  # Get the expression matrix for the genes
  gene1_expr <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')[gene1, ]
  gene2_expr <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')[gene2, ]
  # Get the metadata for the group
  group_meta <- seurat_obj$group2
  # Calculate the ratio of cells co-expressing the genes in each group
  co_expr <- table(group_meta[gene1_expr > 0 & gene2_expr > 0]) / table(group_meta)
  return(co_expr)
}

# Calculate the co-expression ratio for each pair of genes
co_expression_ratios <- list()
gene_pairs <- list(c("Hopx", "Mki67"), c("Eomes", "Mki67"), c("Sox2", "Mki67"), c("Mki67", "Hopx"), c("Mki67", "Eomes"), c("Mki67", "Sox2"), c("Hopx", "Mki67"), c("Eomes", "Hopx"), c("Sox2", "Hopx"), c("Eomes", "Sox2"), c("Sox2", "Eomes"))
for (pair in gene_pairs) {
  co_expression_ratios[[paste(pair, collapse = "_")]] <- co_express(fabio, pair[1], pair[2])
}

# not co-express
not_co_express <- function(seurat_obj, gene1, gene2) {
  # Get the expression matrix for the genes
  gene1_expr <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')[gene1, ]
  gene2_expr <- Seurat::GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')[gene2, ]
  # Get the metadata for the group
  group_meta <- seurat_obj$group2
  # Calculate the ratio of cells not co-expressing the genes in each group
  not_co_expr <- table(group_meta[gene1_expr > 0 & gene2_expr == 0]) / table(group_meta)
  return(not_co_expr)
}

# Calculate the not co-expression ratio for each pair of genes
not_co_expression_ratios <- list()
for (pair in gene_pairs) {
  not_co_expression_ratios[[paste(pair, collapse = "_")]] <- not_co_express(fabio, pair[1], pair[2])
}

#### gene expression ----
fabio$group <- factor(fabio$group, levels = c("zeroM", "oneM", "twoM"))

meta <- fabio$group

hopx <- Seurat::GetAssayData(fabio, assay  = "RNA", slot = "data")["Hopx", ]
mki67 <- Seurat::GetAssayData(fabio, assay  = "RNA", slot = "data")["Mki67", ]
tbr2 <-  Seurat::GetAssayData(fabio, assay  = "RNA", slot = "data")["Eomes", ]
neurod1 <- Seurat::GetAssayData(fabio, assay  = "RNA", slot = "data")["Neurod1", ]

hkratio <- table(meta[hopx > 0.3 & mki67 > 0.3]) / table(meta[hopx>0.3])
tkratio <- table(meta[tbr2 > 0.3 & mki67 > 0.3]) / table(meta[tbr2>0.3])

nd1ratio <- table(meta[neurod1 > 0.3]) / table(meta)
hopxratio <- table(meta[hopx > 0.3]) / table(meta)
tbr2ratio <- table(meta[tbr2 > 0.3]) / table(meta)
mki67ratio <- table(meta[mki67 > 0.3]) / table(meta)


plot_data <- data.frame(
  Group = names(hopxratio),
  Hopx_Mki67 = as.numeric(hkratio),
  Tbr2_Mki67 = as.numeric(tkratio),
  Neurod1 = as.numeric(nd1ratio),
  Hopx = as.numeric(hopxratio),
  Tbr2 = as.numeric(tbr2ratio),
  Mki67 = as.numeric(mki67ratio)
)


plot_data_long <- plot_data %>%
  pivot_longer(cols = -Group, 
               names_to = "Marker", 
               values_to = "Ratio")


ggplot(plot_data_long, aes(x = factor(Group, levels =  c("zeroM", "oneM", "twoM")), y = Ratio, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Marker, scales = "free_y") +
  theme_bw() +
  labs(title = "Marker Expression Ratios by Group",
       x = "Group",
       y = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

DefaultAssay(sub) <- "RNA"
Idents(sub) <- "tracy_clusters"
sub <- NormalizeData(sub)
Seurat::VlnPlot(sub, features = c("Hopx", "Mki67", "Eomes"), layer = "data", )

