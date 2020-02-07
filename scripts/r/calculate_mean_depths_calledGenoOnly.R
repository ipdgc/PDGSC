library(tidyverse)
library(stringr)
library(stringi)
AD_SUBSET <- Sys.getenv("AD_SUBSET")
POSTQC_VARIANTLIST <- Sys.getenv("POSTQC_VARIANTLIST")
AD_NEWFILENAME <- AD_SUBSET %>% str_replace("matrix_calledGenoOnly_chrALL_", "means_nsites_calledGenoOnly")
# postqc_variantlist <- read_delim(POSTQC_VARIANTLIST, delim="\t", col_types = cols(.default="c"), col_names=FALSE)
names(postqc_variantlist)[1] <- "variant"
ad <- read_delim( AD_SUBSET, delim="\t", col_types = cols(.default="c"))
print("Depth data read in!")
# ad <- subset(ad, ID%in%postqc_variantlist$variant)
# print("Depth data subsetted for to only include post-QC variants!")
ad <- ad %>%
    mutate_all(funs(type.convert)) %>%
    mutate_if(is.factor, as.character)
print("Depth data converted to numeric!")
mean_ads <- data.frame(sapply(ad, mean, na.rm=TRUE)[-1])
print("Mean depths calculated!")
n_sites <- data.frame(sapply(ad, function(x) sum(!is.na(x)))[-1])
print("N sites calculated!")
names(mean_ads) <- "mean_depth"
names(n_sites) <- "n_sites"
mean_ads$iid <- rownames(mean_ads)
n_sites$iid <- rownames(n_sites)
combined <- inner_join(mean_ads, n_sites, by = "iid")
write.table(combined, AD_NEWFILENAME_GROUP1, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
