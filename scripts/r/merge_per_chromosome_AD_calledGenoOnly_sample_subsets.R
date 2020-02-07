library(tidyverse)
library(stringr)
library(stringi)
SUBSET <- Sys.getenv("AD_SUBSET")
subset_files <- list.files("AD_calledGenoOnly_per_chromosome_sample_subsets", pattern = SUBSET)
AD_subset_files <- subset_files[str_detect(subset_files, "PDGSC_AD_matrix")]
AD_NEWFILENAME <- str_c("PDGSC_AD_matrix_calledGenoOnly_chrALL_", SUBSET)
AD_full_matrix <- NULL
for (subset in AD_subset_files){
    ad_subset <- read_delim(str_c("AD_calledGenoOnly_per_chromosome_sample_subsets/", subset), delim="\t", col_types = cols(.default="c"))
    print (str_c("File ", subset, " read in!"))
    AD_full_matrix <- bind_rows(AD_full_matrix, ad_subset)
    print (str_c("File ", subset, " merged!"))
}
write.table(AD_full_matrix, str_c("final_AD_calledGenoOnly_files_sample_subsets/", AD_NEWFILENAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(AD_full_matrix)
rm(ad_subset)
gc()
