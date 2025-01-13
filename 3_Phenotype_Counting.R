################
# dependencies #
################

library(stringr)
library(dplyr)

load(paste0(results_folder, "/Cell_matrix_pre_post.RData"))
load(paste0(results_folder, "/Metadata_pre_post.Rdata"))
load(paste0(results_folder, "/Metadata_sample_pre_post.Rdata"))


##########
# Script #
##########

phenotype_list <- all_cells_combined %>% count(phenotype) 
pheno_names <- phenotype_list$phenotype
pheno_names_rep <- str_replace_all(pheno_names,"\\+","_plus_") #note: as the plus isn't recognized, the name is changed to cotain a 'plus' instead of +...
rm(pheno_counts)
pheno_counts <- data.frame()
ROIs <- unique(all_cells_combined$ROI)

#loop through all files and count each phenotype, these are combined into one 'phenocountscomb' file
for (j in 1:length(unique(all_cells_combined$ROI))) { 
  
  ROI_name <- ROIs[j]
  rm(all_cells_subset)
  all_cells_subset <- all_cells_combined %>% filter(grepl(ROI_name, all_cells_combined$ROI))  
  all_cells_subset$phenotype <- str_replace_all(all_cells_subset$phenotype,"\\+","_plus_")
  
  for (t in 1:length(pheno_names_rep)) { 
    
    pheno <- pheno_names_rep[t]
    
    count <- sum(all_cells_subset$phenotype == pheno)
    count_pheno <- data.frame(count, pheno, ROI_name)
    pheno_counts<- rbind(pheno_counts, count_pheno)  
    
  }
  
  print (paste0("counted"," ", ROI_name))
}

print("looped")

#change the 'plus' back to '+' and transform the counts to counts per sample instead of per ROI
pheno_counts$pheno <- str_replace_all(pheno_counts$pheno,"_plus_","+")
pheno_counts <- pheno_counts %>% rename(ROI = ROI_name)
pheno_counts <-inner_join(pheno_counts, metadata, by = c("ROI"))
pheno_counts <- pheno_counts %>% filter(!grepl(samples_to_exlude, ROI)) 

pheno_counts_sample <- pheno_counts %>% 
  group_by(pheno, patient, group, response, pre_anno) %>% 
  summarise(totalsample = sum(count, na.rm = TRUE), 
            meansample = mean(count, na.rm = TRUE),
            mean_stroma = mean(stroma_area, na.rm = TRUE))



pheno_counts_sample$counts_mm_stroma <- pheno_counts_sample$meansample / (pheno_counts_sample$mean_stroma/1000000)
#pheno_counts_sample <-inner_join(pheno_counts_sample, metadata_sample, by = c("sample"))
#pheno_counts_sample$meansample_corrected <- pheno_counts_sample$meansample/pheno_counts_sample$area


#Two tables are generated and saved in the workindirectory containing the phenotype counts per sample and per image
write.csv2(pheno_counts, paste0(results_folder, "/phenocounts_image.csv"), row.names = FALSE)
write.csv2(pheno_counts_sample, paste0(results_folder, "/phenocounts_sample_stroma_cor.csv"), row.names = FALSE)
save(pheno_counts_sample, file = paste0(results_folder, "/phenocounts_sample_stroma_cor.RData"))
save(pheno_counts, file = paste0(results_folder, "/phenocounts_image.RData"))
print("saved phenotype counts")
