library(Seurat)
library(SeuratObject)
library(future)
library(ggplot2)
library(dplyr)
library(magrittr)
library(harmony)
library(spacexr)
library(presto)

combined_xenium_obj <- readRDS("xenium_full_projected_nonormbrain.rds")
DefaultAssay(combined_xenium_obj) <- "Xenium"

#match cells with hazard to xenium object
sample_labels <- c("1-1","1-2","1-3", "1-4", "2-1", "2-2", "2-3", "3-1","3-2","3-3","3-4","4-1","4-2","4-3")
csv_cell_paths <- list(
  "xenium_haz_labels/cells_1-1.csv.gz",
  "xenium_haz_labels/cells_1-2.csv.gz",
  "xenium_haz_labels/cells_1-3.csv.gz",
  "xenium_haz_labels/cells_1-4.csv.gz",
  "xenium_haz_labels/cells_2-1.csv.gz",
  "xenium_haz_labels/cells_2-2.csv.gz",
  "xenium_haz_labels/cells_2-3.csv.gz",
  "xenium_haz_labels/cells_3-1.csv.gz",
  "xenium_haz_labels/cells_3-2.csv.gz",
  "xenium_haz_labels/cells_3-3.csv.gz",
  "xenium_haz_labels/cells_3-4.csv.gz",
  "xenium_haz_labels/cells_4-1.csv.gz",
  "xenium_haz_labels/cells_4-2.csv.gz",
  "xenium_haz_labels/cells_4-3-left.csv.gz",
  "xenium_haz_labels/cells_4-3-right.csv.gz"
)

csv_metadata_list <- lapply(1:length(csv_cell_paths), function(i) {
  csv_path <- csv_cell_paths[[i]]
  if (i <= 13) {
    sample_name <- sample_labels[i]
  } else {
    sample_name <- "868"  # Both 868 files
  }

  metadata <- read.csv(gzfile(csv_path))
  rownames(metadata) <- metadata$cell_id
  metadata$attn <- metadata$attn
  metadata$sample_label <- sample_name
  return(metadata)
})

combined_metadata <- do.call(rbind, csv_metadata_list)

combined_xenium_obj$haz <- NA
combined_xenium_obj$attn <- NA

csv_to_sample_map <- list(
  "1-1" = 1,
  "1-2" = 2,
  "1-3" = 3,
  "1-4" = 4,
  "2-1" = 5,
  "2-2" = 6,
  "2-3" = 7,
  "3-1" = 8,
  "3-2" = 9,
  "3-3" = 10,
  "3-4" = 11,
  "4-1" = 12,
  "4-2" = 13,
  "4-3" = c(14, 15)
)

haz_named <- rep(NA_real_, ncol(combined_xenium_obj))
names(haz_named) <- colnames(combined_xenium_obj)
attn_named <- haz_named

for (sample_name in sample_labels) {
  seurat_cells <- colnames(combined_xenium_obj)[combined_xenium_obj$sample == sample_name]
  original_cell_ids <- gsub(paste0("^", sample_name, "_"), "", seurat_cells)

  if (sample_name == "4-3") {
    sample_metadata <- rbind(csv_metadata_list[[14]], csv_metadata_list[[15]])
  } else {
    sample_metadata <- csv_metadata_list[[csv_to_sample_map[[sample_name]]]]
  }

  cat("Processing", sample_name, "\n")
  cat("  Seurat cells:", length(seurat_cells), "\n")
  cat("  Metadata cells:", nrow(sample_metadata), "\n")

  matching_cells <- intersect(original_cell_ids, rownames(sample_metadata))
  seurat_matched <- seurat_cells[match(matching_cells, original_cell_ids)]

  haz_named[seurat_matched] <- as.numeric(sample_metadata[matching_cells, "haz"])
  attn_named[seurat_matched] <- as.numeric(sample_metadata[matching_cells, "attn"])

  cat("  Matched:", length(matching_cells), "\n")
}

# Use AddMetaData for proper cell ID alignment
combined_xenium_obj <- AddMetaData(combined_xenium_obj, metadata = haz_named, col.name = "haz")
combined_xenium_obj <- AddMetaData(combined_xenium_obj, metadata = attn_named, col.name = "attn")


cat("\n=== MATCHING SUMMARY ===\n")
cat("Total cells:", ncol(combined_xenium_obj), "\n")
cat("Cells with haz scores:", sum(!is.na(combined_xenium_obj$haz)), "\n")
cat("Cells with attn scores:", sum(!is.na(combined_xenium_obj$attn)), "\n")
cat("Percentage matched:", round(100 * sum(!is.na(combined_xenium_obj$haz)) / ncol(combined_xenium_obj), 2), "%\n")


#plot hazard
sample_labels_mut <- c("1-3", "1-4", "2-1", "2-2", "2-3", "3-1","3-2","3-3","3-4","4-1","4-2","4-3")
for(s in sample_labels_mut){
  subset_xenium_obj <- subset(combined_xenium_obj, subset = sample == s)

  haz_values <- as.numeric(as.character(subset_xenium_obj$haz))
  coords <- GetTissueCoordinates(subset_xenium_obj)
  coords$haz <- haz_values

  coords_na    <- coords[is.na(coords$haz), ]
  coords_valid <- coords[!is.na(coords$haz), ]
  coords_valid <- coords_valid[order(coords_valid$haz), ]  # plot low values first
  x_range <- diff(range(coords$x, na.rm = TRUE))
  y_range <- diff(range(coords$y, na.rm = TRUE))
  aspect_ratio <- y_range / x_range
  plot_width <- 8
  plot_height <- max(3, min(12, plot_width * aspect_ratio))  # clamp between 3 and 12


  p1 <- ggplot(coords, aes(x = x, y = y, color = haz)) +
    geom_point(size = 0.5) +
    scale_color_gradientn(
      colors = c("white", "red"),
      limits = c(-4, 0),
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white")
    )
  ggsave(paste0(s,"_haz_heatmap.png"), plot = p1, width = plot_width, height = plot_height, dpi = 300)
}


#get patch-level markers
med_haz <- median(combined_xenium_obj$haz,na.rm = TRUE)
combined_xenium_obj$risk_group <- ifelse(
  combined_xenium_obj$haz > med_haz,
  "HighRisk", "LowRisk"
)
table(combined_xenium_obj$sample, combined_xenium_obj$risk_group)


Idents(combined_xenium_obj) <- "risk_group"
combined_xenium_obj[["Xenium"]] <- JoinLayers(combined_xenium_obj[["Xenium"]])
markers <- FindMarkers(
  combined_xenium_obj,
  ident.1 = "HighRisk",
  ident.2 = "LowRisk",
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0.1  # lower threshold to be less stringent
)
write.csv(markers, "markers_patch_medhaz.csv", row.names = TRUE)
