################
# dependencies #
################

library(stringr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)
library(tidyr)

load(paste0(results_folder, "/Cell_matrix.RData"))
load(paste0(results_folder, "/phenocounts_sample_stroma_cor.RData"))



##########
# Script #
##########

pheno_abundance = subset(pheno_counts_sample, select = -c(totalsample, meansample, mean_stroma))
pheno_abundance <- pheno_abundance %>% filter(!grepl(phenotypes_to_exclude, pheno)) 
pheno_abundance <- pheno_abundance %>% filter(!grepl("patientname", patient)) 

pheno_abundance <- pheno_abundance %>% filter(grepl('post', group)) #select for a specific group
#pheno_abundance <- pheno_abundance %>% filter(!grepl('epit|ibroblast|essel|Myeloid_other', pheno)) #select for specific phenotypes
pheno_abundance <- pheno_abundance %>% filter(grepl('granu|dendritic|mono|macro', pheno)) 

pheno_abundance <- pheno_abundance %>%
  pivot_wider(names_from = pheno, values_from = counts_mm_stroma) %>%
  as.data.frame()

tumour <- pheno_abundance$sample
group <- pheno_abundance$group
response <- pheno_abundance$response
rownames(pheno_abundance) <- pheno_abundance$sample
variable = ncol(pheno_counts_sample) - 2
pheno_abundance = subset(pheno_abundance, select = -c(1:4))

pheno_abundance <- as.matrix(pheno_abundance)
pheno_abundance <- scale(pheno_abundance)
pheno_abundance <- t(pheno_abundance)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-4, 4))

col_ha = columnAnnotation(group = group, response = response,
                          simple_anno_size = unit(3, "mm"), border = TRUE, 
                          col = list(group = c(group_colours),
                                     response = c(group_colours))
)

row_ha = rowAnnotation(Log = anno_boxplot(pheno_abundance, height = unit(4, "cm"),   axis_param = list(labels_rot=270)))


pdf(file = paste0(results_folder, "/pheno_abundance_heatmap_myeloid_post_stroma_norm.pdf"), width = 10, height =10)

if (exists("pheno_order")) {
  Heatmap(pheno_abundance, col = col_fun,
          name = "scaled\ncounts",
          column_title = "Phenotype distribution",
          width = ncol(pheno_abundance) * unit(3, "mm"),
          column_names_gp = gpar(fontsize = 6),
          column_names_rot = 45,
          column_names_centered = FALSE,
          row_names_gp = gpar(fontsize = 6),
          row_names_side = "left",
          show_column_names = TRUE,
          border_gp = gpar(col = "black"),
          cluster_rows = FALSE,
          #top_annotation = col_ha,
          right_annotation = row_ha,
          row_order = pheno_order
  )
} else {
  Heatmap(pheno_abundance, col = col_fun,
          name = "scaled\ncounts",
          column_title = "Phenotype distribution",
          width = ncol(pheno_abundance) * unit(3, "mm"),
          height = nrow(pheno_abundance) * unit(3, "mm"),
          column_names_gp = gpar(fontsize = 6),
          column_names_rot = 45,
          column_names_centered = FALSE,
          row_names_gp = gpar(fontsize = 6),
          row_names_side = "right",
          show_column_names = TRUE,
          border_gp = gpar(col = "black"),
          cluster_rows = TRUE,
          top_annotation = col_ha,
          #right_annotation = row_ha
  )
}

dev.off()

