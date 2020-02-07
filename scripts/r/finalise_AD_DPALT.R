GT_FILE <- Sys.getenv("GT_FILE")
AD_FILE <- Sys.getenv("AD_FILE")
library(tidyverse)
library(stringr)
library(stringi)
AD_NEWFILENAME <- str_c("PDGSC_", AD_FILE) %>% str_replace(".csv", ".tab")
DPALT_NEWFILENAME <- str_c("PDGSC_", GT_FILE) %>% str_replace("genotypes", "DPALT_matrix")
ad <- read_delim(str_c("per_chromosome_AD_DPALT_files/", AD_FILE), delim=",", col_types = cols(.default="c"))
sample_ids <- names(ad)[-1]
ad_sorted <- ad %>% separate(ID, c("CHROM", "POS"), remove=FALSE) %>% mutate(CHROM = NULL) %>% arrange(as.numeric(POS))
rm(ad)
gc()
gt <- read_delim(str_c("allsamples_genotypes/", GT_FILE), delim="\t", col_types = cols(.default="c"))
gt_sorted <- gt %>% mutate(`#CHROM` = str_c(`#CHROM`, POS, sep=":")) %>% arrange(as.numeric(POS))
rm(gt)
gc()
names(gt_sorted)[1] <- "ID"
stopifnot(all(ad_sorted$ID==gt_sorted$ID))
gt_sorted_noids <- gt_sorted[,-c(1,2)]
rm(gt_sorted)
gc()
find_alt <- function(X){
    !stri_detect_regex(X, "^(0\\/0|\\.\\/\\.)$")
}
gt_sorted_noids_logical <- as.data.frame(t(apply(gt_sorted_noids, 1, find_alt)))
rm(gt_sorted_noids)
gc()
ad_sorted_noids <- ad_sorted[,-c(1,2)]
ad_sorted_noids_numeric <- as.data.frame(matrix(NA, dim(ad_sorted_noids)[1],dim(ad_sorted_noids)[2]))
ad_sorted_noids_numeric[] <- lapply(ad_sorted_noids, function(x) as.numeric(as.character(x)))
rm(ad_sorted_noids)
gc()
names(ad_sorted_noids_numeric) <- sample_ids
names(gt_sorted_noids_logical) <- sample_ids
dpalt <- ad_sorted_noids_numeric * gt_sorted_noids_logical
rm(ad_sorted_noids_numeric)
rm(gt_sorted_noids_logical)
gc()
dpalt[dpalt==0] <- NA
dpalt_ids <- bind_cols(ad_sorted[,1], dpalt)
rm(dpalt)
gc()
write.table(dpalt_ids, str_c("final_AD_DPALT_files/", DPALT_NEWFILENAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(dpalt_ids)
gc()
ad_final <- ad_sorted[,-2]
rm(ad_sorted)
gc()
write.table(ad_final, str_c("final_AD_DPALT_files/", AD_NEWFILENAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(ad_final)
gc()