################
# dependencies #
################

library(ggplot2)
library(dplyr)
library(viridis)


if (analysis == 'knn') {
  load(paste0(results_folder, "/Spatial_neighbour_analysis_knn_", knn_n, ".RData"))
  
} else {
  load(paste0(results_folder, "/Spatial_neighbour_analysis_expansion_", expansion_dist, ".RData"))
}




##########
# Script #
##########


test_df <- test_df %>%
  group_by(from_label, to_label, type, response_simple) %>%
  summarise(
    mean_ct = mean(ct, na.rm = TRUE),
    p_mean = mean(p_gt, na.rm = TRUE)
  )

pdf(file = paste0(results_folder, "/interaction_heatmap_R.pdf"), width = 16, height = 8)

ggplot(test_df, aes(x=from_label, y=to_label,  fill=mean_ct)) +
  #scale_y_discrete(limits=pheno_order) + scale_x_discrete(limits=pheno_order) +
  geom_tile() + scale_fill_viridis(option= "viridis", na.value= "yellow", limits=c(0,2)) +
  ylab("Phenotype in neighbourhood") + xlab("Phenotype of interest") +
  guides(fill = guide_colourbar(title= "relative n interactions/cell")) +
  geom_point(aes(colour=p_mean)) +
  facet_wrap(facets = vars(response_simple), scales = "free", ncol = 7) + 
  scale_colour_gradient("P-value", low = "white", high = "black", breaks=c(0,0.05, 0.1), limits=c(0,0.1), na.value="transparent") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()


