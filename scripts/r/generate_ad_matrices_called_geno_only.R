GT_FILE <- Sys.getenv("GT_FILE")
AD_FILE <- Sys.getenv("AD_FILE")
library(tidyverse)
library(stringr)
library(stringi)
AD_NEWFILENAME <- AD_FILE %>% str_replace("matrix_", "matrix_calledGenotypesOnly_")
ad <- read_delim(str_c("final_AD_DPALT_files/", AD_FILE), delim="\t", col_types = cols(.default="c"))
print("AD data loaded in!")
sample_ids <- names(ad)[-1]
ad <- ad %>% separate(ID, c("CHROM", "POS"), remove=FALSE) %>% mutate(CHROM = NULL) %>% arrange(as.numeric(POS))
print("AD data sorted!")
gt <- read_delim(str_c("allsamples_genotypes/", GT_FILE), delim="\t", col_types = cols(.default="c"))
print("GT data loaded in!")
gt <- gt %>% mutate(`#CHROM` = str_c(`#CHROM`, POS, sep=":")) %>% arrange(as.numeric(POS))
print("GT data sorted!")
names(gt)[1] <- "ID"
stopifnot(all(ad$ID==gt$ID))
gt <- gt[,-c(1,2)]
find_called <- function(X){
    !stri_detect_regex(X, "^(\\.\\/\\.)$")
}
gt <- as.data.frame(t(apply(gt, 1, find_called)))
print("Called genotypes found!")
ad_variantids <- ad[,1]
ad <- ad[,-c(1,2)]
ad_sorted_noids_numeric <- as.data.frame(matrix(NA, dim(ad)[1],dim(ad)[2]))
ad_sorted_noids_numeric[] <- lapply(ad, function(x) as.integer(as.character(x)))
print("AD data converted to integer format!")
rm(ad)
gc()
names(ad_sorted_noids_numeric) <- sample_ids
names(gt) <- sample_ids
ad_sorted_noids_numeric <- ad_sorted_noids_numeric * gt
print("AD and GT data has been multiplied!")
rm(gt)
gc()
ad_sorted_noids_numeric[ad_sorted_noids_numeric==0] <- NA
print("Uncalled genotypes set as NA!")
ad_sorted_noids_numeric <- bind_cols(ad_variantids, ad_sorted_noids_numeric)
write.table(ad_sorted_noids_numeric, str_c("final_depths_called_geno_only_nsites/", AD_NEWFILENAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(ad_sorted_noids_numeric)
gc()
