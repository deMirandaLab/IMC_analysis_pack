################
# dependencies #
################

library(stringr)
library(dplyr)
library(SpatialExperiment)
library(imcRtools)
library(pheatmap)

load(paste0(results_folder, "/Cell_matrix.RData"))
load(paste0(results_folder, "/Metadata.RData"))




##########
# Script #
##########

all_cells_subset <- all_cells_combined %>% filter(!grepl(phenotypes_to_exclude, all_cells_combined$pheno_merge ))  
#all_cells_subset <- all_cells_combined %>% filter(grepl("Tcells|macro|mono|Bcell|plasma|granulo|dendritic|ILC|vessels|fibro", all_cells_combined$pheno_merge )) 
counts <- as.matrix(all_cells_subset[,  ..markers])

# Create the DataFrame for spatial coordinates:
meta_data <- data.frame(
  x = all_cells_subset$Location_Center_X,
  y = all_cells_subset$Location_Center_Y,
  Cell_id = all_cells_subset$Cell_id,
  ROI = all_cells_subset$ROI
)

meta_data <-inner_join(meta_data, metadata, by = c("ROI"))

# Make sure the dimensions of assay and colData align correctly
# by transposing the 'y' matrix if needed:
if (ncol(counts) != ncol(meta_data)) {
  counts <- t(counts)
}

# Create the SpatialExperiment object:
spe <- SpatialExperiment(
  assay = counts, 
  colData = meta_data, 
  spatialCoordsNames = c("x", "y")
)


pheno <- data.frame(all_cells_subset$pheno_merge)
colnames(pheno) <- c("pheno")
colData(spe) <- cbind(colData(spe), pheno)


#spatial analysis
spe <- buildSpatialGraph(spe, img_id = "ROI", type = "knn", k= knn_n, coords = c("x", "y"))
spe <- buildSpatialGraph(spe, img_id = "ROI", type = "expansion", threshold = 10, coords = c("x", "y")) 

save(spe, file = paste0(results_folder, "/spatial_dataframe_allcells_expansion_10.Rdata"))



png(paste0(results_folder, "/expansionplot_all_ROIs.png"), 
    width = 20, height = 20, res = 300, units = "in")

plotSpatial(spe[,spe$ROI == "S08_T2_1"],
            coords = c("x", "y"),
            node_color_by = "pheno",
            node_size_fix = 1.5,
            img_id = "ROI", 
            draw_edges = TRUE, 
            colPairName = "expansion_interaction_graph", 
            nodes_first = FALSE,
            edge_color_fix = "white") + #scale_color_manual(values = color_mapping) +
  theme(
    panel.background = element_rect(fill = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

dev.off()

png(paste0(results_folder, "/knnplot_all_ROIs.png"), 
    width = 20, height = 20, res = 300, units = "in")

plotSpatial(spe,
            coords = c("x", "y"),
            node_color_by = "pheno",
            node_size_fix = 1.5,
            img_id = "ROI", 
            draw_edges = TRUE, 
            colPairName = "knn_interaction_graph", 
            nodes_first = FALSE,
            edge_color_fix = "white") + #scale_color_manual(values = color_mapping) +
  theme(
    panel.background = element_rect(fill = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

dev.off()


if (analysis == 'knn') {
  counts <- countInteractions(spe,
                              group_by = "ROI",
                              label = "pheno",
                              colPairName = "knn_interaction_graph",
                              method = c("classic"),
                              patch_size = NULL
  ) %>% as.data.frame()
  
} else {
  counts <- countInteractions(spe,
                            group_by = "ROI",
                            label = "pheno",
                            colPairName = "expansion_interaction_graph",
                            method = c("classic"),
                            patch_size = NULL
  ) %>% as.data.frame()
}

names(counts)[names(counts) == "group_by"] <- "ROI"
counts <- inner_join(counts, metadata, by = c("ROI"))


spe_filtered <- spe[, grepl("pre", colData(spe)$group)]

#test significance
if (analysis == 'knn') {
  test <- testInteractions(spe_filtered,
                           group_by = "ROI",
                           label = "pheno",
                           colPairName = "knn_interaction_graph",
                           method = c("classic"),
                           patch_size = NULL,
                           iter = 500,
                           p_threshold = 0.05,
                           #return_samples = FALSE
                           #tolerance = sqrt(.Machine$double.eps),
                           #BPPARAM = SerialParam()
  )
} else {
  test <- testInteractions(spe_filtered,
                           group_by = "ROI",
                           label = "pheno",
                           colPairName = "expansion_interaction_graph",
                           method = c("classic"),
                           patch_size = NULL,
                           iter = 500,
                           p_threshold = 0.05,
                           #return_samples = FALSE
                           #tolerance = sqrt(.Machine$double.eps),
                           #BPPARAM = SerialParam()
  )
}

test_df <- as.data.frame(test)

names(test_df)[names(test_df) == "group_by"] <- "ROI"
test_df <- inner_join(test_df, metadata, by = c("ROI"))

if (analysis == 'knn') {
  save(test_df, file = paste0(results_folder, "/Spatial_neighbour_analysis_knn_", knn_n, ".RData"))
  
} else {
  save(test_df, file = paste0(results_folder, "/Spatial_neighbour_analysis_expansion_", expansion_dist, "all.RData"))
}
