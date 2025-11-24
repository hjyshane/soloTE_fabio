library(Seurat)
library(qs)

set.seed(83)

options(future.globals.maxSize = 128000 * 1024^2)

base_dir <-  "~/SoloTE/"

sample_list <- list(M2M_88 = c("M2M88", 8),
                    M1M_35 = c("M1M35", 4),
                    M0M_57 = c("M0M57", 7),
                    M0M_8 = c("M0M8", 3),
                    F2M_13 = c("F2M13", 6),
                    F1M_75 = c("F1M75", 2),
                    F0M_28 = c("F0M28", 1),
                    F0M_2 = c("F0M2", 5))


# grap cells
# Read final obj
final <- readRDS("~/Slingshot_Tracy/Fabio_finalobj.rds")

for (i in seq_along(sample_list)) {
  # Directory
  qs_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "qsave")

    # Read object
  sobj <- qs::qread(file.path(qs_dir, paste0("01__seurate_obj_", sample_list[[i]][1], "_locustes.qs")))
  
  # Find from final object
  fobj <- subset(final, subset = sample == names(sample_list[i]))
  
  # Clean barcode
  rownames(fobj@meta.data) <- sub(paste0("_", sample_list[[i]][2]), "", rownames(fobj@meta.data))
  
  # Get intersect
  keep <- intersect(colnames(sobj), rownames(fobj@meta.data))
  sobj_filtered <- subset(sobj, cells = keep)
  
  # add metadata
  ord <- match(Seurat::Cells(sobj_filtered), rownames(fobj@meta.data))
  meta <- fobj@meta.data[ord, , drop = FALSE]
  sobj_filtered@meta.data <- cbind(sobj_filtered@meta.data, meta)
  
  # save
  qs::qsave(sobj_filtered, file.path(file.path(qs_dir, paste0("02__filtered_obj_", sample_list[[i]][1], "_locustes.qs"))))

}
