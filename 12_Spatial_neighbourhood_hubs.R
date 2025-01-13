################
# dependencies #
################

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)
library(SpatialExperiment)
library(imcRtools)
library(ggplot2)
library(tidyr)


load(paste0(results_folder,"/spatial_dataframe_allcells_expansion_20.Rdata"))




##########
# Script #
##########

if (analysis == 'knn') {
  spe <- aggregateNeighbors(spe, 
                            colPairName = "knn_interaction_graph", 
                            aggregate_by = "metadata", 
                            count_by = "pheno")
  
} else {
  spe <- aggregateNeighbors(spe, 
                            colPairName = "expansion_interaction_graph", 
                            aggregate_by = "metadata", 
                            count_by = "pheno")
}



set.seed(220705)

pretreatment_idx <- colData(spe)$group == "pre"

# Subset the aggregatedNeighbors for pretreatment samples
pretreatment_agg <- spe$aggregatedNeighbors[pretreatment_idx, ]


cn_1 <- kmeans(pretreatment_agg, centers = 30)
spe$cn_celltypes <- NA
spe$cn_celltypes[pretreatment_idx] <- as.factor(cn_1$cluster)

#spe$cn_celltypes <- as.factor(cn_1$cluster)
spe$cn_celltypes <- as.character(spe$cn_celltypes)


plotSpatial(spe[,spe$ROI == "ROI_1"], 
            coords = c("x", "y"),
            node_color_by = "cn_celltypes", 
            img_id = "ROI", 
            node_size_fix = 1.5) +
  #scale_color_manual(values = color_mapping) + 
  theme(panel.background = element_rect(fill = "black"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
  )












          #####################################
          #                                   #
          #       CLEAN UP NEXT PART          #
          #                                   #
          #####################################

for_plot <- prop.table(table(spe$cn_celltypes, spe$'pheno'), 
                       margin = NULL)

for_plot <- t(for_plot)
for_plot <- scale(for_plot)


#for_plot[for_plot < 0] <- 0  
#for_plot[for_plot > 4] <- 4 

pdf(file = paste0(results_folder, "/IMC_expression_hubs_30.pdf"), width = 10, height = 8)

Heatmap(for_plot, 
        col = viridis(100, option = "viridis"),
        row_title = "Cellular neighbourhoods",
        column_title = "phenotype clusters",
        column_title_side = "bottom",
        row_title_side = "left",
        show_column_names = TRUE,
        show_row_names = TRUE,
        cluster_columns = TRUE,
        cluster_rows = TRUE)


dev.off()




sample_plot <- table(spe$cn_celltypes, spe$ROI) %>% as.data.frame
names(sample_plot)[names(sample_plot) == "Var2"] <- "ROI"
names(sample_plot)[names(sample_plot) == "Var1"] <- "cn_celltypes"

sample_plot <- merge(sample_plot, metadata, by = "ROI")
sample_plot <- sample_plot %>% group_by(patient, cn_celltypes, group, response) %>% summarise(Freq = mean(Freq, na.rm=TRUE))
total <- sample_plot %>% group_by(patient, group, response) %>% summarise(totalcells = sum(Freq, na.rm=TRUE))
sample_plot <- merge(sample_plot, total, by = c("patient", "group"))
sample_plot$percentage <- (sample_plot$Freq / sample_plot$totalcells) * 100

sample_plot <- sample_plot %>% filter(!grepl("31|41", patient)) 
sample_plot <- sample_plot %>% filter(patient != "8")

sample_plot <- sample_plot %>% filter(grepl("pre", group))


plot <- sample_plot

plot$response.x <- factor(plot$response.x  , levels = c("R", "NR"))
group_colors_trans <- c(alpha("green", 0.5), alpha("red", 0.5), alpha("red", 0.5), alpha("red", 0.5))

ncol <- 4 
num_facets <- length(unique(plot$cn_celltypes))
nrow <- ceiling(num_facets / ncol)
total_height <- nrow * 2.5

# Step 5: Create the PDF with dynamic height
pdf(file = paste0(results_folder, "/cn_pre_immune_graphs_perc_30.pdf"), width = 8, height = total_height)


ggplot(plot, aes(x=response.x, y= percentage, color=response.x)) + 
  #geom_line(aes(group = patient), color = "grey")+ 
  geom_point(color = "black") + scale_color_manual(values=c("CR" = "black", "NR" = "black")) +
  facet_wrap(facets = vars(cn_celltypes), scales="free_y", ncol=ncol) + 
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
    aes(fill = response.x),
    size = 0.5
  ) +
  scale_fill_manual(values = group_colors_trans) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))



dev.off()






counts_samples <- sample_plot[,-5] %>%
  pivot_wider(names_from = cn_celltypes, values_from = percentage, values_fill = 0) %>% as.data.frame()

group = counts_samples$group
patient = counts_samples$patient
response = counts_samples$response.x

col_ha = columnAnnotation(group = group, response = response,
                          simple_anno_size = unit(3, "mm"), border = TRUE, 
                          col = list(group = c(group_colours),
                                     response = c(group_colours))
)




counts_samples$patient <- as.character(counts_samples$patient)
rownames(counts_samples) <- counts_samples$patient
counts_samples <- as.matrix(counts_samples[,-c(1:5)])
counts_samples <- t(counts_samples)
counts_samples <- scale(counts_samples)
#counts_samples[counts_samples < -2] <- -2  
#counts_samples[counts_samples > 2] <- 2 
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-4, 4))


pdf(file = paste0(results_folder, "/Hubs_30_perc_heatmap_patient.pdf"), width = 7, height = 5)


Heatmap(counts_samples, 
        col = col_fun,
        row_title = "Cellular neighbourhoods",
        column_title = "FOV",
        show_column_names = TRUE,
        column_title_side = "bottom",
        row_title_side = "left",
        show_row_names = TRUE,
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        #right_annotation = row_ha,
        top_annotation = col_ha)

dev.off()

