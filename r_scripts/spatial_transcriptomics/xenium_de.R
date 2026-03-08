library(Seurat)
library(SeuratObject)
library(future)
library(ggplot2)
library(dplyr)
library(magrittr)
library(harmony)
library(spacexr)
library(presto)


combined_xenium_obj <- readRDS("xenium_objects_with_rctd.rds")

combined_xenium_obj <- JoinLayers(combined_xenium_obj)

#DEG BASED ON CELL TYPE
hr_samples <- c("2-3","175","214","510","304")
lr_samples <- c("1-4","2-1","2-2","1-3","249","686","868")
combined_xenium_obj$group <- ifelse(combined_xenium_obj$sample %in% hr_samples, "HR",
                                    ifelse(combined_xenium_obj$sample %in% lr_samples, "LR", "WT"))

# Normalize
combined_xenium_obj <- NormalizeData(combined_xenium_obj)
combined_xenium_obj <- FindVariableFeatures(combined_xenium_obj)

#DE
Idents(combined_xenium_obj) <- "group"
markers_groups <- FindAllMarkers(combined_xenium_obj, only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.5)
write.table(markers_groups, file = "subtype_markers.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

deg_results <- FindMarkers(combined_xenium_obj, ident.1 = "HR", ident.2 = "LR")
write.table(deg_results, file = "hr_vs_lr_sample_DEG_xenium.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)

weight_matrix <- as.matrix(combined_xenium_obj@meta.data[, c("rctd_Immune", "rctd_NormalBrain", "rctd_Tumour")])
max_indices <- max.col(weight_matrix)
combined_xenium_obj$celltype <- c("Immune", "NormalBrain", "Tumour")[max_indices]


table(combined_xenium_obj$group)
Idents(combined_xenium_obj) <- "celltype"
celltypes_to_analyze <- c("Immune", "Tumour", "NormalBrain")

for(ct in celltypes_to_analyze) {
  deg_results <- FindMarkers(
    combined_xenium_obj,
    ident.1 = "HR",
    ident.2 = "LR",
    group.by = "group",
    subset.ident = ct,
    min.pct = 0.05,
    logfc.threshold = 0.5
  )

  filename <- paste0(tolower(ct), "_hr_vs_lr_DEG_xenium.csv")
  write.table(deg_results,
              file = filename,
              quote = FALSE, sep = ",",
              row.names = TRUE, col.names = NA)

  print(paste("Completed DEG for", ct, "- found", nrow(deg_results), "genes"))
}
