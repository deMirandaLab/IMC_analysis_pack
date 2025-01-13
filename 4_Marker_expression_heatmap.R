################
# dependencies #
################

#library(stringr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)


load(paste0(results_folder,"/Cell_matrix_pre_post.RData"))


##########
# Script #
##########

cluster_expression <- all_cells_combined %>% group_by(pheno_merge) %>% summarise(across(everything(), mean), .groups = 'drop') %>%
  as.data.frame()

cluster_expression <- cluster_expression %>% filter(!grepl(phenotypes_to_exclude, pheno_merge))

rownames(cluster_expression) <- cluster_expression$pheno_merge
cluster_expression = subset(cluster_expression, select = c(markers+1))
cluster_expression_matrix <- as.matrix(cluster_expression)
#cluster_expression_matrix <- scale(cluster_expression_matrix)

col_fun = colorRamp2(c(0.05, 0.25, 0.5), c("darkblue", "yellow", "red"))
col_fun(seq(-6, 6))


pdf(file = paste0(results_folder, "/category_expression_heatmap.pdf"), width = 12, height = 8)
if (exists("pheno_order")) {
  Heatmap(cluster_expression_matrix,
          col = col_fun,
          name = "relative\nintensity",
          column_title = "Marker expression per phenotype cluster",
          width = ncol(cluster_expression_matrix) * unit(4, "mm"),
          column_names_gp = gpar(fontsize = 10),
          row_names_gp = gpar(fontsize = 10),
          row_order = pheno_order
  )
} else {
  Heatmap(cluster_expression_matrix,
          col = col_fun,
          name = "relative\nintensity",
          column_title = "Marker expression per phenotype cluster",
          width = ncol(cluster_expression_matrix) * unit(3, "mm"),
          column_names_gp = gpar(fontsize = 10),
          row_names_gp = gpar(fontsize = 10),
          cluster_rows = TRUE
  )
}

dev.off()
