library(Seurat)
library(future)
library(ggplot2)
library(dplyr)
library(magrittr)
options(future.globals.maxSize = 2000 * 1024^2)

#csv from xenium explorer to exclude a section for another case
csv_file <- '2-3_cells_include.csv'
cell_ids_to_split <- read.csv(csv_file)$Cell_ID

#attention + hazard scores
csv_cell_paths <- list("xenium_case_split/cells_1-3.csv.gz","xenium_case_split/cells_1-4.csv.gz","xenium_case_split/cells_2-1.csv.gz","xenium_case_split/cells_2-2.csv.gz","xenium_case_split/cells_2-3.csv.gz")
csv_metadata_list <- lapply(csv_cell_paths, function(csv_path) {
  metadata <- read.csv(gzfile(csv_path))
  rownames(metadata) <- metadata$cell_id  # Set cell IDs as rownames
  median_attn <- median(metadata$attn, na.rm = TRUE)
  metadata$attn_flag <- ifelse(metadata$attn > median_attn, 1, 0)
  return(metadata)
})
gene_list_file <- "tumour_genes.txt"
tumour_markers <- readLines(gene_list_file)
xenium_paths <- list(
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027509__Region_1__20240807__201856/output-XETG00098__0027509__Region_1__20240807__201856',
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027509__Region_2__20240807__201856/output-XETG00098__0027509__Region_2__20240807__201856',
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027509__Region_3__20240807__201856/output-XETG00098__0027509__Region_3__20240807__201856',
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027509__Region_4__20240807__201856/output-XETG00098__0027509__Region_4__20240807__201856',
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027512__Region_1__20240807__201856/output-XETG00098__0027512__Region_1__20240807__201856',
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027512__Region_2__20240807__201856/output-XETG00098__0027512__Region_2__20240807__201856',
  '/Users/katherinerich/Downloads/LGG-FINALIZED/data/xenium_output/output-XETG00098__0027512__Region_3__20240807__201856/output-XETG00098__0027512__Region_3__20240807__201856'
)

xenium_objects <- lapply(1:length(xenium_paths), function(i) {
  xenium_data <- LoadXenium(xenium_paths[[i]], fov = "fov", assay = "Xenium")
  xenium_data <- subset(xenium_data, subset = nCount_Xenium > 5)
  expression_matrix <- GetAssayData(xenium_data, slot = "counts")[tumour_markers, ] > 0
  cells_expressing_any_gene <- colSums(expression_matrix) > 0
  xenium_data <- subset(xenium_data, cells = colnames(xenium_data)[cells_expressing_any_gene])
  if(i!=1&i!=2){
  metadata <- csv_metadata_list[[i-2]]
  xenium_data <- subset(xenium_data, cells = metadata$cell_id)
  xenium_data <- AddMetaData(xenium_data, metadata = metadata,col.name ='haz')
  xenium_data <- AddMetaData(xenium_data, metadata = metadata,col.name ='attn_flag')
  xenium_data <- AddMetaData(xenium_data, metadata = metadata,col.name ='attn')
  }
  
  if(i==7){
    xenium_data <- subset(xenium_data, cells = cell_ids_to_split)
  }
  return(xenium_data)
})

sample_labels <- c("1-1","1-2","1-3", "1-4", "2-1", "2-2", "2-3")
for (i in seq_along(xenium_objects)) {
  xenium_objects[[i]]$sample <- sample_labels[i]
}

#orig_xenium_obj <- merge(xenium_objects[[1]], y = xenium_objects[-1], add.cell.ids = sample_labels, project = "combined_xenium")
#combined_xenium_obj <- subset(orig_xenium_obj, subset = sample %in% c("1-1","1-2") | attn_flag == 1)
combined_xenium_obj <- merge(xenium_objects[[1]], y = xenium_objects[-1], add.cell.ids = sample_labels, project = "combined_xenium")
wt_samples <- c("1-1","1-2")
odg_samples <- c("1-3","1-4","2-2","2-3")
astro_samples <- c("2-1")
combined_xenium_obj$subtype <- ifelse(combined_xenium_obj$sample %in% odg_samples, "ODG", ifelse(combined_xenium_obj$sample %in% astro_samples, "Astro", "WT"))


combined_xenium_obj$group <- ifelse(combined_xenium_obj$haz > -3.065747 & combined_xenium_obj$subtype != "WT" , "HR", ifelse(combined_xenium_obj$subtype == "WT", "WT", "LR"))
#combined_xenium_obj$group <- ifelse(combined_xenium_obj$haz > (-3.120119+0.5) & combined_xenium_obj$subtype != "WT" , "HR", ifelse(combined_xenium_obj$haz < (-3.120119-0.5) & combined_xenium_obj$subtype != "WT", "LR",ifelse(combined_xenium_obj$subtype == "WT", "WT", "MR") ))

Idents(combined_xenium_obj) <- "group"

combined_xenium_obj@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = ..count..),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")


combined_xenium_obj@meta.data %>% 
  ggplot(aes(x=group, fill=group)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = ..count..),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Risk Group")

metadata <- combined_xenium_obj@meta.data
metadata %>% 
  ggplot(aes(color=sample, x=nCount_Xenium, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 

metadata <- combined_xenium_obj@meta.data
metadata %>% 
  ggplot(aes(color=group, x=nCount_Xenium, fill= group)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() 

combined_xenium_obj <- NormalizeData(combined_xenium_obj)
combined_xenium_obj <- FindVariableFeatures(combined_xenium_obj)

#DEG
Idents(combined_xenium_obj) <- "group"
combined_xenium_obj <- JoinLayers(combined_xenium_obj)
subset_xenium_obj <- subset(combined_xenium_obj, subset = attn_flag==1)
deg_results <- FindMarkers(subset_xenium_obj, ident.1 = "HR", ident.2 = "LR")
write.table(deg_results, file = "hr_vs_lr_roi_DEG_xenium_median.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)

hr_samples <- c("2-3")
lr_samples <- c("1-4","2-1","2-2","1-3")
combined_xenium_obj$sample_risk <- ifelse(combined_xenium_obj$sample %in% hr_samples, "HR", ifelse(combined_xenium_obj$sample %in% lr_samples, "LR", "WT"))
Idents(combined_xenium_obj) <- "sample_risk"
subset_xenium_obj <- subset(combined_xenium_obj, subset = attn_flag==1)
deg_results <- FindMarkers(subset_xenium_obj, ident.1 = "HR", ident.2 = "LR")
write.table(deg_results, file = "hr_vs_lr_sample_DEG_xenium_median.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)


combined_xenium_obj <- SketchData(
  object = combined_xenium_obj,
  ncells = 30000,
  method = "LeverageScore",
  sketched.assay = "sketch")

DefaultAssay(combined_xenium_obj) <- "sketch"


combined_xenium_obj <- FindVariableFeatures(combined_xenium_obj,assay="sketch", verbose = F)
top10 <- head(VariableFeatures(combined_xenium_obj), 10)
plot1 <- VariableFeaturePlot(combined_xenium_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(filename = "plot2-2.png", plot = plot2, width = 6, height = 4, dpi = 300)


combined_xenium_obj <- ScaleData(combined_xenium_obj,assay="sketch", verbose = F)
combined_xenium_obj <- RunPCA(combined_xenium_obj,npcs=30,assay="sketch", verbose = F)
elbow_plot <- ElbowPlot(combined_xenium_obj,ndims = 30)
ggsave(filename = "elbow_plot.png", plot = elbow_plot, width = 6, height = 4, dpi = 300)
DimHeatmap(combined_xenium_obj, dims = 10:18, cells = 500, balanced = TRUE,ncol=3)

combined_xenium_obj <- FindNeighbors(combined_xenium_obj,k.param=30,assay="sketch", dims = 1:30)
combined_xenium_obj <- FindClusters(combined_xenium_obj,algorithm=4, resolution = 0.5)
combined_xenium_obj <- RunUMAP(combined_xenium_obj,n.neighbors = 30L,assay="sketch",dims = 1:30,min.dist=0.1)

umap_plot <- DimPlot(combined_xenium_obj, group.by = "seurat_clusters") + ggtitle("UMAP Clustering")
ggsave("umap_clustering.png", plot = umap_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "group") + ggtitle("UMAP: HR vs LR")
ggsave("clustering_group.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "sample") + ggtitle("UMAP: Samples")
ggsave("clustering_sample.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "subtype") + ggtitle("UMAP: Subtype")
ggsave("clustering_subtype.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)
cluster_plot <- DimPlot(combined_xenium_obj, group.by = "sample_risk") + ggtitle("UMAP: Sample Risk")
ggsave("clustering_sample_risk.png", plot = cluster_plot, width = 8, height = 6, dpi = 300)


cluster_counts <- table(combined_xenium_obj$seurat_clusters, combined_xenium_obj$sample)
write.csv(cluster_counts, file = "cluster_counts_case.csv", row.names = TRUE)
chisq_test <- chisq.test(cluster_counts)
print(chisq_test)


cluster_counts <- table(combined_xenium_obj$seurat_clusters,combined_xenium_obj$sample, combined_xenium_obj$group)
write.csv(cluster_counts, file = "cluster_counts_group_sample.csv", row.names = TRUE)
cluster_counts <- table(combined_xenium_obj$seurat_clusters, combined_xenium_obj$group)
write.csv(cluster_counts, file = "cluster_counts.csv", row.names = TRUE)
#combined_xenium_obj <- JoinLayers(combined_xenium_obj)
combined_xenium_obj[["sketch"]] <- JoinLayers(combined_xenium_obj[["sketch"]])

cluster_proportions <- prop.table(cluster_counts, margin = 2)
print(cluster_proportions)

write.csv(cluster_proportions, file = "cluster_proportions.csv", row.names = TRUE)
Idents(combined_xenium_obj) <- "group"
markers_clusters <- FindAllMarkers(combined_xenium_obj,assay="sketch", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.5)
write.table(markers_clusters, file = "group_markers.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

Idents(combined_xenium_obj) <- "seurat_clusters"
markers_clusters <- FindAllMarkers(combined_xenium_obj,assay="sketch", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.5)
write.table(markers_clusters, file = "cluster_markers.csv", sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)


DefaultAssay(combined_xenium_obj) <- "Xenium"
combined_xenium_obj <- JoinLayers(combined_xenium_obj)

Idents(combined_xenium_obj) <- "group"
deg_results <- FindMarkers(combined_xenium_obj, ident.1 = "HR", ident.2 = "LR")
write.table(deg_results, file = "hr_vs_lr_roi_DEG_xenium.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)

Idents(combined_xenium_obj) <- "sample_risk"
deg_results <- FindMarkers(combined_xenium_obj, ident.1 = "HR", ident.2 = "LR")
write.table(deg_results, file = "hr_vs_lr_sample_DEG_xenium.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)


Idents(combined_xenium_obj) <- "group"
gene_list_file <- "genes_roi.txt"

top_genes <- readLines(gene_list_file)
print("TOP GENES FROM FILE")
print(top_genes)
print(dim(combined_xenium_obj))
for(assay in SeuratObject::Assays(combined_xenium_obj)) {
  if("SCTAssay" %in% class(combined_xenium_obj[[assay]])) {
    slot(object = combined_xenium_obj[[assay]],
         name = "SCTModel.list")[[1]]@feature.attributes %<>%
      `[`(intersect(rownames(.), top_genes), )
  }
}

DefaultAssay(combined_xenium_obj) <- "sketch"

filtered_xenium_obj <- subset(x=combined_xenium_obj,downsample = 10000,features=top_genes)
print(dim(filtered_xenium_obj))

heatmap_plot <- DoHeatmap(filtered_xenium_obj, features = top_genes, group.by = "sample")
ggsave("heatmap_cases.png", plot = heatmap_plot, width = 8, height = 15, dpi = 300)

heatmap_plot <- DoHeatmap(filtered_xenium_obj, features = top_genes, group.by = "group")
ggsave("heatmap.png", plot = heatmap_plot, width = 8, height = 15, dpi = 300)

heatmap_plot <- DoHeatmap(filtered_xenium_obj, features = top_genes, group.by = "seurat_clusters")
ggsave("heatmap_clusters.png", plot = heatmap_plot, width = 25, height = 15, dpi = 300)

heatmap_plot <- DoHeatmap(filtered_xenium_obj, features = top_genes, group.by = "sample_risk")
ggsave("heatmap_sample_risk.png", plot = heatmap_plot, width = 10, height = 15, dpi = 300)

heatmap_plot <- DoHeatmap(filtered_xenium_obj, features = top_genes, group.by = "subtype")
ggsave("heatmap_subtype.png", plot = heatmap_plot, width = 10, height = 15, dpi = 300)

DefaultAssay(combined_xenium_obj) <- "sketch"



Idents(combined_xenium_obj) <- "seurat_clusters"
sample_to_plot <- "2-3"  # Change this to your desired sample
subset_xenium_obj <- subset(combined_xenium_obj, subset = sample == sample_to_plot)
p1=ImageDimPlot(subset_xenium_obj, size = 0.75)
ggsave("2-3.png", plot = p1, width = 8, height = 8, dpi = 300)

head(subset_xenium_obj [[]])

# Extract the row names (index) and the 'seurat_clusters' column
selected_data <- subset_xenium_obj@meta.data %>% mutate(Index = rownames(subset_xenium_obj@meta.data)) %>%  # Add row names as a new column
  select(Index, seurat_clusters)         # Select the Index and seurat_clusters columns

# Remove the first 4 characters from the 'Index' column
selected_data <- selected_data %>%
  mutate(Index = substr(Index, 5, nchar(Index)))
# Save the modified data frame to a compressed CSV file
write.csv(selected_data, file = gzfile("2-3.csv.gz"), row.names = FALSE)

tibble(cellID = colnames(combined_xenium_obj), clusterID = Idents(combined_xenium_obj)) %>%
  write_csv(file = sprintf("cluster_mappings.csv", Sys.Date()))

#Idents(combined_xenium_obj) <- "attn_flag"
#PLOTTINGCLUSTER1
sample_to_plot <- "2-3"  # Change this to your desired sample
subset_xenium_obj <- subset(combined_xenium_obj, subset = sample == sample_to_plot & seurat_clusters%in% c(7,10))
#subset_xenium_obj <- subset(combined_xenium_obj, subset = sample == sample_to_plot)
#cluster_colors <- rep("gray80", length(unique(subset_xenium_obj$seurat_clusters)))
#cluster_colors["12"] <- "red"
#p1=ImageDimPlot(subset_xenium_obj, size = 0.75,cols = cluster_colors)
p1=ImageDimPlot(subset_xenium_obj, size = 0.75)
ggsave("2-3_7-10.png", plot = p1, width = 8, height = 6, dpi = 300)

Idents(combined_xenium_obj) <- "group"
sample_to_plot <- "2-2"  # Change this to your desired sample
subset_xenium_obj <- subset(combined_xenium_obj, subset = sample == sample_to_plot)
group_colors <- c("HR" = "red", "LR" = "blue", "Low-Attention" = "gray")
p1=ImageDimPlot(subset_xenium_obj, size = 0.75,cols=group_colors)
ggsave("2-2_group.png", plot = p1, width = 8, height = 6, dpi = 300)

subset_xenium_obj <- subset(combined_xenium_obj, subset = sample %in% c("1-3","1-4","2-1","2-2","2-3") & seurat_clusters%in% c(1))

xenium_obj <- subset(subset_xenium_obj, subset = attn_flag==1)
print(mean(xenium_obj$haz))

subset_xenium_obj@meta.data %>%
  filter(!is.na(haz) & !is.na(attn)) %>%
  ggplot(aes(x = haz, fill = factor(attn_flag))) + 
  geom_density(alpha = 0.6) + 
  scale_fill_manual(values = c("blue", "red"), labels = c("Low Attention", "High Attention")) +
  geom_vline(xintercept = -3.06, linetype = "dashed", color = "black", size = 0.5) +  # Add vertical line
  theme_classic() +
  labs(fill = "Attention Level")

normalized_attn <- subset_xenium_obj$attn / sum(subset_xenium_obj$attn)

# Multiply normalized attn by the corresponding haz values
weighted_haz <- normalized_attn * subset_xenium_obj$haz

# Calculate the sum of the weighted haz values
sum_weighted_haz <- sum(weighted_haz)

# Print the result
print(sum_weighted_haz)


