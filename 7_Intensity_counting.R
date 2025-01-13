################
# dependencies #
################

library(stringr)
library(dplyr)

load(paste0(results_folder, "/Cell_matrix.Rdata"))
load(paste0(results_folder, "/phenocounts_sample.RData"))


##########
# Script #
##########

phenotype_list <- all_cells_combined %>% count(pheno_merge) 
pheno_names <- phenotype_list$pheno_merge
pheno_names_rep <- str_replace_all(pheno_names,"\\+","_plus_") #note: as the plus isn't recognized, the name is changed to cotain a 'plus' instead of +...
all_cells_combined$Image_id <- as.character(all_cells_combined$Image_id)
rm(checkpoints_counts_combined)
checkpoints_counts_combined <- data.frame()
ROIs <- unique(all_cells_combined$ROI)

#loop through all files and count each phenotype, these are combined into one 'phenocountscomb' file

for (u in 1:length(marker_names)) {
  rm(checkpoint_counts)
  checkpoint_counts <- data.frame()
  marker <- marker_names[u]
  
  for (j in 1:length(ROIs)) { 
    

    ROI_name <- ROIs[j]
    
    
    rm(all_cells_subset)
    all_cells_subset <- all_cells_combined %>% filter(grepl(ROI_name, all_cells_combined$ROI))  
    all_cells_subset$pheno_merge <- str_replace_all(all_cells_subset$pheno_merge,"\\+","_plus_")
    all_cells_subset <- all_cells_subset %>% filter(.data[[marker]] >= 0.2)  #select the threshold, standard use 0.2
    
    for (t in 1:length(pheno_names_rep)) { 
      
      pheno <- pheno_names_rep[t]
      
      count <- sum(all_cells_subset$pheno_merge == pheno)
      count_pheno <- data.frame(count, pheno, ROI_name, marker)
      checkpoint_counts<- rbind(checkpoint_counts, count_pheno)  
      
    }
    
    #print (paste0("counted ", ROI_name))
  }
  checkpoint_counts <- checkpoint_counts %>% filter(!grepl(samples_to_exlude, ROI_name))
  checkpoint_counts$totalphenotype <- pheno_counts$count
  checkpoints_counts_combined <- rbind(checkpoints_counts_combined, checkpoint_counts)
  print(paste0("looped ", marker))
  
}

print("finished generating checkpoint files")


#change the 'plus' back to '+' and transform the counts to counts per sample instead of per ROI
checkpoints_counts_combined$pheno <- str_replace_all(checkpoints_counts_combined$pheno,"_plus_","+")
checkpoints_counts_combined <- checkpoints_counts_combined %>% rename(ROI = ROI_name)
checkpoints_counts_combined <-inner_join(checkpoints_counts_combined, metadata, by = c("ROI"))

checkpoints_counts_sample <- checkpoints_counts_combined %>% group_by(pheno, sample, marker) %>% summarise(meanpositive = mean(count, na.rm=TRUE))
tempmean <- (checkpoints_counts_combined %>% group_by(pheno, sample, marker) %>% summarise(meansample = mean(totalphenotype, na.rm=TRUE)))
checkpoints_counts_sample$totalpheno <- tempmean$meansample
checkpoints_counts_sample <-inner_join(checkpoints_counts_sample, metadata_sample, by = c("sample"))


#Two tables are generated and saved in the workindirectory containing the phenotype counts per sample and per image
write.csv2(checkpoints_counts_combined, paste0(results_folder, "/intensity_counts_image.csv"), row.names = FALSE)
write.csv2(checkpoints_counts_sample, paste0(results_folder, "/intensity_counts_sample.csv"), row.names = FALSE)
save(checkpoints_counts_sample, file = paste0(results_folder, "/intensity_counts_sample.RData"))
save(checkpoints_counts_combined, file = paste0(results_folder, "/intensity_counts_image.RData"))
print("saved checkpoint counts")
