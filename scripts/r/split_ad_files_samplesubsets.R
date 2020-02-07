library(tidyverse)
library(stringr)
library(stringi)
DPALT_FILE <- Sys.getenv("DPALT_FILE")
AD_FILE <- Sys.getenv("AD_FILE")
AD_NEWFILENAME <- AD_FILE %>% str_replace(".tab", "_samplesubset")
DPALT_NEWFILENAME <- DPALT_FILE %>% str_replace(".tab", "_samplesubset")
DATASETS <- c(AD_FILE, DPALT_FILE)
for (DATASET in DATASETS){
    if(DATASET == AD_FILE) OUTPUT_FILENAME <- AD_NEWFILENAME
    if(DATASET == DPALT_FILE) OUTPUT_FILENAME <- DPALT_NEWFILENAME
    ad_matrix <- read_delim(str_c("final_AD_DPALT_files/", DATASET), delim="\t", col_types = cols(.default="c"))
    sample_ids <- names(ad_matrix)
    n_ids <- length(sample_ids)
    for (id_start in seq(2, n_ids, 750)){
        id_end <- id_start+749
        ad_subset <- ad_matrix[,c(1,id_start:id_end)]
        write.table(ad_subset, str_c("final_AD_DPALT_files_per_chromosome_sample_subsets/", OUTPUT_FILENAME, "_", id_start, ".tab"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        rm(ad_subset)
        gc()
    }
rm(ad_matrix)
gc()
}