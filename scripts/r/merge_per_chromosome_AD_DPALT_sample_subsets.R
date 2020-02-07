library(tidyverse)
library(stringr)
library(stringi)
SUBSET <- Sys.getenv("AD_SUBSET")
subset_files <- list.files("final_AD_DPALT_files_per_chromosome_sample_subsets", pattern = SUBSET)
AD_subset_files <- subset_files[str_detect(subset_files, "PDGSC_AD_matrix")]
DPALT_subset_files <- subset_files[str_detect(subset_files, "PDGSC_DPALT_matrix")]
AD_NEWFILENAME <- str_c("PDGSC_AD_matrix_chrALL_", SUBSET)
DPALT_NEWFILENAME <- str_c("PDGSC_DPALT_matrix_chrALL_", SUBSET)
AD_full_matrix <- NULL
for (SUBSET in AD_subset_files){
    ad_subset <- read_delim(str_c("final_AD_DPALT_files_per_chromosome_sample_subsets/", SUBSET), delim="\t", col_types = cols(.default="c"))
    AD_full_matrix <- bind_rows(AD_full_matrix, ad_subset)
}
write.table(AD_full_matrix, str_c("final_AD_DPALT_files_sample_subsets/", AD_NEWFILENAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(AD_full_matrix)
rm(ad_subset)
gc()
DPALT_full_matrix <- NULL
for (SUBSET in DPALT_subset_files){
    dpalt_subset <- read_delim(str_c("final_AD_DPALT_files_per_chromosome_sample_subsets/", SUBSET), delim="\t", col_types = cols(.default="c"))
    DPALT_full_matrix <- bind_rows(DPALT_full_matrix, dpalt_subset)
}
write.table(DPALT_full_matrix, str_c("final_AD_DPALT_files_sample_subsets/", DPALT_NEWFILENAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(DPALT_full_matrix)
rm(dpalt_subset)
gc()