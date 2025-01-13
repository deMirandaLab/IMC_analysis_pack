#### Dependencies ####

library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(data.table)
library(tibble)
library(dplyr)

results_folder <- paste0(getwd(), "/Results" )  # Full path to the subfolder
if (!file.exists(results_folder)) {
  dir.create(results_folder)
}

#required data files
load(paste0(results_folder, "/intensity_counts_sample.RData"))



checkpoint_adapt <- checkpoints_counts_sample

#checkpoint_adapt <- checkpoint_adapt %>% filter(grepl("ibroblasts|essels", pheno))

checkpoint_adapt <- checkpoint_adapt %>% filter(grepl("TGF", marker))

checkpoint_comb <- checkpoint_adapt %>% 
  group_by(patient.x, group.x, marker, response.x) %>% 
  summarise(totalpheno = sum(totalpheno, na.rm = TRUE), 
            meanpositive = mean(meanpositive, na.rm = TRUE),
            mean_stroma = mean(stroma_area, na.rm = TRUE))

checkpoint_comb$counts_mm_stroma <- checkpoint_comb$meanpositive / (checkpoint_comb$mean_stroma/1000000)

checkpoint_comb$perc_positive <- (checkpoint_comb$meanpositive/checkpoint_comb$totalpheno)*100

#write.csv2(checkpoint_comb, paste0(results_folder, "/CD4_checkpoint_counts.csv"), row.names = FALSE)

#pdf(file = paste0(results_folder, "/Checkpoint_cells_mm2_", marker2,".pdf"), width = 15, height = 3)

group_colors_trans <- c(alpha("green", 0.5), alpha("green", 0.5), alpha("red", 0.5), alpha("red", 0.5))

checkpoint_comb$group_tp  <- paste(checkpoint_comb$response.x, checkpoint_comb$group.x, sep = " ")
checkpoint_comb$group_tp <- factor(checkpoint_comb$group_tp , levels = c("R", "NR", "control"))

ncol <- 4 
num_facets <- length(unique(checkpoint_comb$marker))
nrow <- ceiling(num_facets / ncol)
total_height <- nrow * 2.5

# Step 5: Create the PDF with dynamic height
pdf(file = paste0(results_folder, "/checkpoint_perc_tgf.pdf"), width = 8, height = total_height)



ggplot(checkpoint_comb, aes(x=group_tp, y= perc_positive, color=group.x)) + 
  geom_line(aes(group = patient.x), color = "grey")+ 
  geom_point(color = "black") + scale_color_manual(values=c("CR" = "black", "NR" = "black")) +
  facet_wrap(facets = vars(marker), scales="free_y", ncol=ncol) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=45, size=6, hjust=1), 
        strip.text=element_text(size=5, color="black"),
        strip.background = element_blank(),
        legend.position = "none") + 
  stat_summary(
    fun.data = function(x) {
      median <- median(x)
      lower_quartile <- quantile(x, 0.25)
      upper_quartile <- quantile(x, 0.75)
      data.frame(y = median, ymin = lower_quartile, ymax = upper_quartile)
    },
    geom = "crossbar",
    position = position_dodge(width = 0.9),
    width = 0.5,
    #color = "black",
    #fill = (alpha("white", 0.25)),
    aes(fill = group_tp),
    size = 0.5
  ) +
  scale_fill_manual(values = group_colors_trans) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))



dev.off()

  
