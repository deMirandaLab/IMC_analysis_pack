################
# dependencies #
################

library(ggplot2)
library(dplyr)

load(paste0(results_folder, "/Cell_matrix.Rdata"))


######################
# Variables          #
# adapt for own data #
######################

image_to_plot = "ROI_1"




##########
# Script #
##########

all_cells_sample <- all_cells_combined %>% filter(grepl(image_to_plot, ROI))
all_cells_sample <- all_cells_sample %>% mutate(row_number = row_number())


ggplot(all_cells_sample, aes(x = Location_Center_X, y = 1000-Location_Center_Y, color = major_groups)) +
  geom_point() +
  scale_color_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "black"),
    #legend.position = "none",
    panel.grid = element_blank()
  )

