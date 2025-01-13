################
# dependencies #
################

library(stringr)
library(ggplot2)
library(dplyr)

load(paste0(results_folder, "/catcounts_sample.RData"))



##########
# Script #
##########

phenotypes_to_plot <- pheno_counts_sample
phenotypes_to_plot <- phenotypes_to_plot %>% filter(!grepl(phenotypes_to_exclude, pheno)) 

phenotypes_to_plot <- phenotypes_to_plot %>% filter(!grepl("epithel|fibro|vessel|undef", pheno)) 

group_colors_trans <- c(alpha("green", 0.5), alpha("green", 0.5), alpha("red", 0.5), alpha("red", 0.5))

phenotypes_to_plot$group_tp  <- paste(phenotypes_to_plot$response, phenotypes_to_plot$group, sep = " ")
phenotypes_to_plot$group_tp <- factor(phenotypes_to_plot$group_tp , levels = c("R", "NR", "control"))

ncol <- 4 
num_facets <- length(unique(phenotypes_to_plot$pheno))
nrow <- ceiling(num_facets / ncol)
total_height <- nrow * 2.5

# Step 5: Create the PDF with dynamic height
pdf(file = paste0(results_folder, "/cat_abundance_graphs_stroma_cor.pdf"), width = 8, height = total_height)



ggplot(phenotypes_to_plot, aes(x=group_tp, y= counts_mm_stroma, color=group)) + 
  geom_line(aes(group = patient), color = "grey")+ 
  geom_point(color = "black") + scale_color_manual(values=c("CR" = "black", "NR" = "black")) +
  facet_wrap(facets = vars(pheno), scales="free_y", ncol=ncol) + 
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




