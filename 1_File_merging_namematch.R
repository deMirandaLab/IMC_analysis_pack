rm(list = ls())

################
# dependencies #
################

library(stringr)
library(data.table)
library(dplyr)


results_folder <- paste0(getwd(), "/Results" )  # Full path to the subfolder

# Create the subfolder if it doesn't exist
if (!file.exists(results_folder)) {
  dir.create(results_folder)
}


######################
# Variables          #
# adapt for own data #
######################

category_folder <- "[add_path]"
phenotype_folder <- "[add_path]"
cell_intensity_folder <- "[add_path]"
metadata_path <- "[add_path]"
coords_path <- "[add_path]"





##########
# Script #
##########

category_files <- list.files(category_folder, pattern = ".csv", full.names = TRUE)    
phenotype_files <- list.files(phenotype_folder, pattern = ".csv", full.names = TRUE)
cell_intensity_files <- list.files(cell_intensity_folder, pattern = ".csv", full.names = TRUE)
metadata <- read.csv((metadata_path), header = TRUE, sep = ";")
metadata <- metadata %>% rename(ROI = id)
#coords_files <- list.files(coords_folder, pattern = ".csv", full.names = TRUE) 
coords <- read.csv((coords_path), header = TRUE, sep = ",")
coords <- coords %>% 
  mutate(ROI = str_extract(PathName_nucleus, "[^\\\\]+$") %>% str_remove("\\.ome$"))
coords <- coords %>% select(Location_Center_X, Location_Center_Y, ROI, ObjectNumber)


all_cells <- data.frame()

for (i in 1:length(phenotype_files)) { 
  
  ROI <- metadata$ROI[i]
  
  
  phenotype <- fread(phenotype_files[grep(ROI, phenotype_files)], header = FALSE, sep = ",")
  category <- fread(category_files[grep(ROI, category_files)], header = FALSE, sep = ",")   
  metrics<- fread(cell_intensity_files[grep(ROI, cell_intensity_files)], header = TRUE)                        
  metrics <- metrics %>% mutate(ObjectNumber = row_number())
  metrics$phenotype <- phenotype$V1
  metrics$category <- category$V1
  metrics <- cbind(metrics, ROI)
  
  #print(paste0(ROI," ", nrow(phenotype)," ", nrow(metrics), " ", nrow(coords)))
  
  all_cells <- rbind(all_cells, metrics)
  
  
}


all_cells_combined  <- inner_join(all_cells, coords, by = c("ROI", "ObjectNumber"))
all_cells_combined <-inner_join(all_cells_combined, metadata, by = c("ROI"))
save(all_cells_combined, file = paste0(results_folder, "/Cell_matrix.Rdata"))

metadata_sample <- metadata %>% select(-ROI) %>% distinct()
save(metadata, file = paste0(results_folder, "/Metadata.Rdata"))
save(metadata_sample, file = paste0(results_folder, "/Metadata_sample.Rdata"))
print("created & saved cell matrix")


