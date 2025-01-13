# IMC_analysis_pack
Software &amp; Workflow for the analysis of imaging mass cytometry data

<br/>

This repository contains the standard pipelines and scripts used by the cancer immunogenomics group for the analysis of imaging mass cytometry data. 
Detailed descriptions of the steps to follow are available in the [IMC analysis Guide](IMC_analysis_guidelines_CIG_v10.pdf)

In short, IMC data is exported from the MCD-viewer following the steps detailed in the guide. Next, data can be normalised using [PENGUIN](https://github.com/deMirandaLab/PENGUIN) and cell segmentation is performed using cellprofiler and the dedicated [cellprofiler pipeline](cp_pipeline_2024.cpproj). 
After data visualisation and cell phenotyping, using ImaCytE and Cytosplore, further analyses can be performed in R-studio.

<br/>

## Scripts
The repository consists of the following downstream visualisation/analysis scripts:
- [1_File_merging_namematch](1_File_merging_namematch.R)	To combine phenotypes, intensities, cell coordinates and metadata in a large data frame for subsequent analysis
- [2_Variables](2_Variables.R) Contains all variables required for downstream scripts
- [3_Phenotype_Counting](3_Phenotype_Counting.R)	To count the abundance of each phenotype per image and sample
- [4_Marker_expression_heatmap](4_Marker_expression_heatmap.R)	To create a heatmap of marker expression per phenotype 
- [5_Phenotype_abundance_heatmap](5_Phenotype_abundance_heatmap.R)	To create a heatmap with the abundance of each phenotype per sample
- [6_Phenotype_abundance_graphs](6_Phenotype_abundance_grap)	To create graphs for each phenotype and their abundance between groups of samples based on provided metadata
- [7_Intensity_counting](7_Intensity_counting.R)	To count the abundance of each phenotype positive for markers of interest 
- [8_Markerpositivity_graphs](8_markerpositivity_graphs.R)	To create graphs of the counts generated in script 6
- [9a_Phenotype_image_visualisation](9a_Phenotype_image_visualisation.R)	To visualise the phenotypes onto a selected image
- [9b_Phenotype_image_visualisation_all_ROIs](9c_Phenotype_image_visualisation_all_ROIs_PNG.R)	To visualise the phenotypes on all images
- [10_Spatial_neighbour_detection](10_Spatial_neighbour_detection.R)	To identify cell-cell neighbourhoods
- [11_Spatial_neighbourhood_plots](11_Spatial_neighbourhood_plots.R)	To plot cell-cell neighbourhood abundances
- [12_Spatial_neighbourhood_hubs](12_Spatial_neighbourhood_hubs.R)	To identify multicellular neighbourhoods and plot those


<br/>

## Overview of the pipeline 
![image](https://github.com/user-attachments/assets/e73e6824-e772-458b-8074-551e499bfe10)
