library(tidyverse)
library(stringr)
library(stringi)
AD_FILE <- Sys.getenv("AD_FILE")
OUTPUT_FILENAME <- AD_FILE %>% str_replace(".tab", "_samplesubset")
ad_matrix <- read_delim(str_c("final_depths_called_geno_only_nsites/", AD_FILE), delim="\t", col_types = cols(.default="c"))
sample_ids <- names(ad_matrix)
n_ids <- length(sample_ids)
for (id_start in seq(2, n_ids, 750)){
    id_end <- id_start+749
    ad_subset <- ad_matrix[,c(1,id_start:id_end)]
    write.table(ad_subset, str_c("AD_calledGenoOnly_per_chromosome_sample_subsets/", OUTPUT_FILENAME, "_", id_start, ".tab"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    rm(ad_subset)
    gc()
}
rm(ad_matrix)
gc()
