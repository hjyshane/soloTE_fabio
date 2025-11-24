library(Seurat)
library(qs)

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

# list of object
obj_list <- list()

for (i in seq_along(sample_list)) {
  # Directory
  qs_dir <- file.path(base_dir, sample_list[[i]][1], paste0(sample_list[[i]][1], "_SoloTE_output"), "qsave")

  # Read object
  obj_list[[names(sample_list[i])]] <- qs::qread(file.path(qs_dir, paste0("03__processed_obj_", sample_list[[i]][1], "_locustes.qs")))

  # rename cells to prevent dup cell names
  obj_list[[names(sample_list[i])]] <- RenameCells(obj_list[[names(sample_list[i])]], add.cell.id = i)

}

# prepare merge list
merge_list <- obj_list[2:8]

# merge
mobj <- merge(x = obj_list[[1]], y = merge_list)

# check cell num
print(table(mobj@meta.data$sample))

# save
qs::qsave(mobj, file.path(directory[["qsave_dir"]], "01__merged_solote_locustes.qs"))
