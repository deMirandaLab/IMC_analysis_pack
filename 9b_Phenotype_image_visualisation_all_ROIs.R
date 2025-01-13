################
# dependencies #
################

library(ggplot2)
library(dplyr)

load(paste0(results_folder, "/Cell_matrix.Rdata"))


##########
# Script #
##########


pdf(file = paste0(results_folder, "/phenotype_image_distribution.pdf"), width = 15, height = 200)

ggplot(all_cells_combined, aes(x = Location_Center_X, y = 1000- Location_Center_Y, color = major_groups)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = palette) +
  facet_wrap(~ ROI, ncol = 3) +
  #scale_color_gradient(low = "black", high = "white") +
  theme(
    panel.background = element_rect(fill = "black"),
    #legend.position = "none",
    panel.grid = element_blank()
  )

dev.off()


