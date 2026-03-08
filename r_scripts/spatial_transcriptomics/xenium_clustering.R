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
hr_samples <- c("2-3","3-1","3-2","4-1","3-4")
lr_samples <- c("1-4","2-1","2-2","1-3","3-3","4-2","4-3")
combined_xenium_obj$group <- ifelse(combined_xenium_obj$sample %in% hr_samples, "HR", ifelse(combined_xenium_obj$sample %in% lr_samples, "LR", "WT"))

astro_samples <- c("2-1","3-2","4-1","3-4","3-3","4-2","4-3")
oligo_samples <- c("1-4","2-3","2-2","1-3","3-1")
combined_xenium_obj$subtype <- ifelse(combined_xenium_obj$sample %in% astro_samples, "astro", ifelse(combined_xenium_obj$sample %in% oligo_samples, "oligo", "WT"))

weight_matrix <- as.matrix(combined_xenium_obj@meta.data[, c("rctd_Immune", "rctd_NormalBrain", "rctd_Tumour")])
max_indices <- max.col(weight_matrix)
combined_xenium_obj$celltype <- c("Immune", "NormalBrain", "Tumour")[max_indices]

#subset to exclude NormalBrain
combined_xenium_obj <- subset(combined_xenium_obj, subset = !(celltype == "NormalBrain"))


Idents(combined_xenium_obj) <- "group"
combined_xenium_obj <- NormalizeData(combined_xenium_obj)
combined_xenium_obj <- FindVariableFeatures(combined_xenium_obj)

ncol(combined_xenium_obj@assays$Xenium)

combined_xenium_obj <- SketchData(
  object = combined_xenium_obj,
  ncells = 40000,
  method = "LeverageScore",
  sketched.assay = "sketch")

ncol(combined_xenium_obj@assays$sketch)

DefaultAssay(combined_xenium_obj) <- "sketch"
combined_xenium_obj <- FindVariableFeatures(combined_xenium_obj, verbose = FALSE)
combined_xenium_obj <- ScaleData(combined_xenium_obj, verbose = FALSE)
combined_xenium_obj <- RunPCA(combined_xenium_obj, verbose = FALSE)

elbow_plot <- ElbowPlot(combined_xenium_obj,ndims = 50)
ggsave(filename = "elbow_plot_rctd_nonormbrain.png", plot = elbow_plot, width = 6, height = 4, dpi = 300)
DimHeatmap(combined_xenium_obj, dims = 1:9, cells = 500, balanced = TRUE,ncol=3)


combined_xenium_obj <- RunUMAP(combined_xenium_obj, dims = 1:30,return.model = TRUE, verbose = FALSE)
combined_xenium_obj <- FindNeighbors(combined_xenium_obj, dims = 1:30,k.param = 30)
combined_xenium_obj <- FindClusters(combined_xenium_obj,algorithm=4, resolution = 0.4)

saveRDS(combined_xenium_obj, "xenium_sketch_clustered_nonormbrain.rds")
combined_xenium_obj <- readRDS("xenium_sketch_clustered_nonormbrain.rds")


umap_plot <- DimPlot(combined_xenium_obj, group.by = "celltype") + ggtitle("UMAP Clustering")
ggsave("umap_clustering_celltype_rctd_nornormbrain.png", plot = umap_plot, width = 8, height = 6, dpi = 300)

umap_plot <- DimPlot(combined_xenium_obj, group.by = "seurat_clusters") + ggtitle("UMAP Clustering")
ggsave("umap_clustering_rctd_nornormbrain.png", plot = umap_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "group") + ggtitle("UMAP: HR vs LR")
ggsave("clustering_group_rctd_nornormbrain.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "sample") + ggtitle("UMAP: Samples")
ggsave("clustering_sample_rctd_nornormbrain.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "subtype") + ggtitle("UMAP: Subtype")
ggsave("clustering_subtype_rctd_nornormbrain.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)


#project to full dataset
library(future)
options(future.globals.maxSize = 3000 * 1024^2)
#to do: projection to full dataset, cell niche
combined_xenium_obj <- ProjectData(
  object = combined_xenium_obj,
  assay = "Xenium",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  full.reduction = "full.pca",
  dims = 1:20,
  refdata = list(seurat_clusters = "seurat_clusters")  # This maps cluster labels
)

DefaultAssay(combined_xenium_obj) <- "sketch"
combined_xenium_obj <- ProjectUMAP(
  query = combined_xenium_obj,
  query.reduction = "full.pca",
  reference = combined_xenium_obj,
  reference.reduction = "pca",
  reduction.model = "umap",
  reduction.name = "full.umap"
)

# verify all cells have embeddings
nrow(combined_xenium_obj@reductions$full.umap@cell.embeddings)


# Verify all cells have cluster labels
table(is.na(combined_xenium_obj$seurat_clusters))
cat("Total cells:", ncol(combined_xenium_obj), "\n")
cat("Cells with cluster labels:", sum(!is.na(combined_xenium_obj$seurat_clusters)), "\n")

saveRDS(combined_xenium_obj, "xenium_full_projected_nonormbrain.rds")
combined_xenium_obj <- readRDS("xenium_full_projected_nonormbrain.rds")

DefaultAssay(combined_xenium_obj) <- "Xenium"
names(combined_xenium_obj@reductions)
nrow(combined_xenium_obj@reductions$full.umap@cell.embeddings)

#DE - cluster markers
combined_xenium_obj[["Xenium"]] <- JoinLayers(combined_xenium_obj[["Xenium"]])
all_markers<- presto::wilcoxauc(combined_xenium_obj, assay = "data",seurat_assay = "Xenium")
write.csv(all_markers, file = "cluster_markers_presto.csv", row.names = TRUE)

#example feature plot
ft_plot <- FeaturePlot(combined_xenium_obj, features = c("APOE"),reduction = "umap")
ggsave("feature_plot_APOE.png", plot = ft_plot, width = 8, height = 6, dpi = 300)

#get counts
cluster_counts <- table(combined_xenium_obj$seurat_clusters, combined_xenium_obj$group)
write.csv(cluster_counts, file = "cluster_counts_group_rctd_baysor_nonormbrain.csv", row.names = TRUE)
chisq_test <- chisq.test(cluster_counts)
print(chisq_test)

cluster_counts <- table(combined_xenium_obj$seurat_clusters, combined_xenium_obj$subtype)
write.csv(cluster_counts, file = "cluster_counts_subtype_rctd_baysor_nonormbrain.csv", row.names = TRUE)
chisq_test <- chisq.test(cluster_counts)
print(chisq_test)

cluster_counts <- table(combined_xenium_obj$seurat_clusters, combined_xenium_obj$sample)
write.csv(cluster_counts, file = "cluster_counts_case_rctd_baysor_nonormbrain.csv", row.names = TRUE)
chisq_test <- chisq.test(cluster_counts)
print(chisq_test)

#check if enriched in certain groups HRvsLR IDHmutant
library(speckle)
#subset to only IDHmutant samples
subset_xenium_obj <- subset(combined_xenium_obj, subset = group %in% c("HR","LR"))
prop <- propeller(clusters = subset_xenium_obj$seurat_clusters, sample = subset_xenium_obj$sample,
                  group = subset_xenium_obj$group, robust=TRUE,trend = TRUE)

write.csv(prop, "enriched_clusters_by_group_hr_lr.csv")
