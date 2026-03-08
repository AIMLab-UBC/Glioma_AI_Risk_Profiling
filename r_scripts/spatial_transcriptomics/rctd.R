library(Seurat)
library(SeuratObject)
library(future)
library(ggplot2)
library(dplyr)
library(magrittr)
library(harmony)
library(spacexr)
library(presto)

#load in single cell reference see https://bioc.r-universe.dev/spacexr/doc/rctd-tutorial.html#reference
ref.sc <- readRDS("rctd_reference_scp503.rds")

xenium_paths <- list(
  '/xenium_data_baysor/1-1',
  '/xenium_data_baysor/1-2',
  '/xenium_data_baysor/1-3',
  '/xenium_data_baysor/1-4',
  '/xenium_data_baysor/2-1',
  '/xenium_data_baysor/2-2',
  '/xenium_data_baysor/2-3',
  '/xenium_data_baysor/3-1',
  '/xenium_data_baysor/3-2',
  '/xenium_data_baysor/3-3',
  '/xenium_data_baysor/3-4',
  '/xenium_data_baysor/4-1',
  '/xenium_data_baysor/4-2',
  '/xenium_data_baysor/4-3'
)

xenium_objects <- lapply(1:length(xenium_paths), function(i) {

  xenium_data <- LoadXenium(xenium_paths[[i]], fov = "fov", assay = "Xenium")
  #threshold set at 25 to remove low-quality cells
  xenium_data <- subset(xenium_data, subset = nCount_Xenium > 25)

  return(xenium_data)
})

sample_labels <- c("1-1","1-2","1-3", "1-4", "2-1", "2-2", "2-3", "3-1","3-2","3-3","3-4","4-1","4-2","4-3")
for (i in seq_along(xenium_objects)) {
  xenium_objects[[i]]$sample <- sample_labels[i]
}

createPuck <- function(idx) {
  xenium_obj <- xenium_objects[[idx]]

  # Get coordinates
  coor.df <- GetTissueCoordinates(xenium_obj, which = "centroids")
  rownames(coor.df) <- coor.df$cell
  coor.df <- coor.df[, c("x", "y")]

  # Get counts
  st.df <- GetAssayData(xenium_obj, assay = "Xenium", layer = "counts")

  # Get nUMI
  numi.df <- xenium_obj$nCount_Xenium

  # Create puck
  puck <- SpatialRNA(coor.df, st.df, numi.df)
  return(puck)
}

idx.list <- 1:length(xenium_objects)
puck.list <- lapply(idx.list, createPuck)
names(puck.list) <- sample_labels

#processing can take a while, saving here to avoid losing progress
saveRDS(puck.list, "puck_list.rds")
saveRDS(xenium_objects, "xenium_objects.rds")


dir.create("rctd_results_baysor", showWarnings = FALSE)

runRCTD <- function(idx) {
  name <- sample_labels[[idx]]
  if(idx < 12){
    name <- sample_labels[[idx]]
    rctd.obj <- readRDS(paste0("rctd_results_baysor/",name,"_rctd_results.rds"))
    return(rctd.obj)
  }
  rctd.obj <- create.RCTD(puck.list[[idx]], ref.sc,UMI_min=10,UMI_min_sigma=50, max_cores = 2)
  rctd.obj <- run.RCTD(rctd.obj, doublet_mode = 'doublet')
  results.df <- rctd.obj@results$results_df
  weights.df <- rctd.obj@results$weights %>% as.data.frame()

  # Add cell type weights
  results.df$Immune <- weights.df$Immune
  results.df$NormalBrain <- weights.df$NormalBrain
  results.df$Tumour <- weights.df$Tumour

  # Save to CSV
  write.csv(
    results.df,
    file.path("rctd_results_baysor", paste0(name, ".rctd.csv"))
  )
  saveRDS(rctd.obj, paste0(name,"_rctd_results.rds"))
  return(rctd.obj)
}

# Run RCTD
rctd.results <- lapply(idx.list, runRCTD)
names(rctd.results) <- sample_labels

# Save RCTD results object
saveRDS(rctd.results, "rctd_results.rds")

#add RCTD labels
for (i in 1:length(xenium_objects)) {
  sample_name <- sample_labels[i]
  results.df <- rctd.results[[sample_name]]@results$results_df
  weights.df <- rctd.results[[sample_name]]@results$weights %>% as.data.frame()

  # Match cells
  cell_names <- colnames(xenium_objects[[i]])

  # Add weights
  xenium_objects[[i]]$rctd_Immune <- weights.df[cell_names, "Immune"]
  xenium_objects[[i]]$rctd_NormalBrain <- weights.df[cell_names, "NormalBrain"]
  xenium_objects[[i]]$rctd_Tumour <- weights.df[cell_names, "Tumour"]
  #xenium_objects[[i]]$celltype <- weights.df[cell_names, "first_type"]
}

# Save labeled xenium objects
combined_xenium_obj <- merge(xenium_objects[[1]], y = xenium_objects[-1], add.cell.ids = sample_labels, project = "combined_xenium")
saveRDS(combined_xenium_obj, "xenium_objects_with_rctd.rds")
