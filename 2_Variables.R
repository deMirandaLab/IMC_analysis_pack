rm(list = ls())

############################################################
#                                                          #
#     This script contains all changeable variables used   # 
#     in further scripts. Make sure to change and load     #
#     before running the scripts                           #
#                                                          #
############################################################

library(scales)

results_folder <- paste0(getwd(), "/Results_pre_post" )  

if (!file.exists(results_folder)) {
  dir.create(results_folder)
}

load(paste0(results_folder, "/Metadata_sample.RData"))
load(paste0(results_folder, "/Cell_matrix_pre_post.Rdata"))




##############################
# Plotting variables         #
##############################

samples_to_exlude = c("sample_name")

markers = c(5, 11:46, 50:52)  # Change to intensity columns in all_cells_combined (or cell_matrix) with markers
marker_names = colnames(all_cells_combined)[markers]

phenotypes_to_exclude = c('rtefact|negative|nclustered|ndefined') # Change to phenotypes to exclude from the analyses

pheno_order = c("Cancercells_Ki67-",            #if a specific order for plotting is desired in the heatmaps or graphs, make here a list of that order
                "Cancercells_Ki67+",
                "Tcells_CD4+", 
                "Tcells_CD8+",
                "Tcells_FOXP3+",
                "ILC",
                "Bcells", "PlasmaBcells_CD138-", "PlasmaBcells_CD138+",
                "HLADR+cells", "Monocytes_other", "Monocytes_HLADR+",
                "Monocytes_CD163+","monocytes_CD163+CD204+", "Monocytes_CD163+CD204+HLADR+", 
                "Macrophages_HLADR+",  "Macrophages_CD204+", "Macrophages_CD204+HLADR+",
                "Macrophages_CD163+CD204+", "Macrophages_CD163+CD204_HLADR+",  
                "Macrophages_sp", "Granulocytes", 
                "Fibroblasts","Vessels")


group_colours <- c("NR" = "red",
                   "R" = "green",
                   "pre" = "black",
                   "post" = "white")

group_colours_trans <- c("R" = alpha("green", 0.5),   #in case you want your colours somewhat transparant, change here
                         "NR" = alpha("red", 0.5))

group <- metadata_sample$group           #annotate which column contains your group




##############################
#       Major lineages       #
##############################

all_cells_combined$major_groups <- all_cells_combined$pheno_merge

all_cells_combined <- all_cells_combined %>%        #if not all groups are represented here, add these
  mutate(major_groups = case_when(
    grepl('plasmabcells', major_groups) ~ 'PlasmaBcells',
    grepl('Bcells', major_groups) ~ 'Bcells',
    grepl('CD8+|CD4+', major_groups) ~ 'Tcells',
    grepl('ILC', major_groups) ~ 'ILC',
    grepl('granulocytes', major_groups) ~ 'Granulocytes',
    grepl('acrophages', major_groups) ~ 'Macrophages',
    grepl('onocytes', major_groups) ~ 'Monocytes',
    grepl('dendritic', major_groups) ~ 'Dendritic_cells',
    grepl('epithelial', major_groups) ~ 'Epithelial',
    grepl('ibroblasts|essels', major_groups) ~ 'Stromal cells',
    TRUE ~ 'Undefined'
  ))


palette <- c("Bcells" = 'yellow',
             "Epithelial" = '#666666',
             "Granulocytes" = 'darkgreen', 
             "ILC" = 'darkred',
             "Macrophages" = 'blue',
             "Monocytes" = 'lightblue',
             "Tcells" = 'red', 
             "Stromal cells" = '#a6761d', 
             "Dendritic cells" = 'cyan',
             "PlasmaBcells" = 'purple',
             "Undefined" = 'black') #add some more colours if more groups were added



##############################
# Spatial analysis variables #
##############################

analysis = "expansion"  #choose knn or expansion

expansion_dist =20      #choose distance for expansion analysis 
knn_n = 5               #choose number of neighbours for knn


