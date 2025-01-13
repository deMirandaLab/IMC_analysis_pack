library(ggplot2)
library(dplyr)


# Loop through each ROI
for (roi in 1:length(unique(all_cells_combined$ROI))) { 
  # Filter data for the current ROI
  name = rois[roi]
  
  all_cells_sample <- all_cells_combined %>% 
    filter(grepl(name, ROI)) %>%
    mutate(row_number = row_number())
  
  # Create the plot
  plot <- ggplot(all_cells_sample, aes(x = Location_Center_X, y = 1000 - Location_Center_Y, color = major_groups)) +
    geom_point() +
    scale_color_manual(values = palette) +
    theme(
      panel.background = element_rect(fill = "black"),
      #legend.position = "none",
      panel.grid = element_blank()
    )
  
  # Save the plot as a PNG file
  ggsave(
    filename = paste0(results_folder, "/pheno_images/plot_", name, ".png"), # Save file named with ROI
    plot = plot,
    width = 8, height = 6, dpi = 300, # Adjust dimensions and resolution as needed
    bg = "transparent" # Optional: set background transparency
  )
}
