library(tidyverse)
library(stringr)
DIFFDEPTH_FILE <- Sys.getenv("DIFFDEPTH_FILE")
DIFFMISS_FILE <- Sys.getenv("DIFFMISS_FILE")
HWE_FILE <- Sys.getenv("HWE_FILE")
FREQ_COUNTS_FILE <- Sys.getenv("FREQ_COUNTS_FILE")
MISSRATE_FILE <- Sys.getenv("MISSRATE_FILE")
MISSRATE_CASES_FILE <- Sys.getenv("MISSRATE_CASES_FILE")
MISSRATE_CONTROLS_FILE <- Sys.getenv("MISSRATE_CONTROLS_FILE")
MULTIALLELIC_FILE <- Sys.getenv("MULTIALLELIC_FILE")
STAR_ALLELES_FILE <- Sys.getenv("STAR_ALLELES_FILE")
FREQ_CASES_FILE <- Sys.getenv("FREQ_CASES_FILE")
FREQ_CONTROLS_FILE <- Sys.getenv("FREQ_CONTROLS_FILE")
ADSP_INTERSECTION_FILE <- Sys.getenv("ADSP_INTERSECTION_FILE")
OUTPUT_FILENAME_QC_MASTERLIST <- Sys.getenv("OUTPUT_FILENAME_QC_MASTERLIST")
OUTPUT_FILENAME_COVARIATES <- Sys.getenv("OUTPUT_FILENAME_COVARIATES")
OUTPUT_FILENAME_PCA <- Sys.getenv("OUTPUT_FILENAME_PCA")
FILENAME_STEM_PERSUBSTUDY_DEPTHS <- Sys.getenv("FILENAME_STEM_PERSUBSTUDY_DEPTHS")
PLINK_PERSUBSTUDY_STEM <- Sys.getenv("PLINK_PERSUBSTUDY_STEM")
OUTPUT_FILENAME_QC_MASTERLIST_DEPTHS <- Sys.getenv("OUTPUT_FILENAME_QC_MASTERLIST_DEPTHS")
OUTPUT_FILENAME_QC_MASTERLIST_SUBSTUDIES <- Sys.getenv("OUTPUT_FILENAME_QC_MASTERLIST_SUBSTUDIES")
samplesizes <- read_delim("../clustering_results_keepmono_50depth_50anc_2to1/samplesizes.txt", delim="\t")

list_of_substudies <- read_delim("../clustering_results_keepmono_50depth_50anc_2to1/list_of_substudies.txt", delim="\t", col_types = cols(.default="c"), col_names=FALSE)$X1
list_of_substudies <- list_of_substudies[list_of_substudies != "group2"]
list_of_substudies <- c("Baylor", "cambridge", "Luebeck", "majamaa", "oslo", "oxford", "tuebingen", "UCL1", "UCL2")

diffdepth <- read_delim(DIFFDEPTH_FILE, delim = "\t", col_types = cols(.default="c"))
diffmiss <- read_delim(DIFFMISS_FILE, delim = "\t", col_types = cols(.default="c"))
hwe <- read_delim(HWE_FILE, delim = "\t", col_types = cols(.default="c"))
freq_counts <- read_delim(FREQ_COUNTS_FILE, delim = "\t", col_types = cols(.default="c"))
missrate <- read_delim(MISSRATE_FILE, delim = "\t", col_types = cols(.default="c"))
missrate_cases <- read_delim(MISSRATE_CASES_FILE, delim = "\t", col_types = cols(.default="c"))
missrate_controls <- read_delim(MISSRATE_CONTROLS_FILE, delim = "\t", col_types = cols(.default="c"))
multiallelic <- read_delim(str_c("../", MULTIALLELIC_FILE), delim = "\t", col_types = cols(.default="c"), col_names=FALSE)
star_alleles <- read_delim(str_c("../", STAR_ALLELES_FILE), delim = "\t", col_types = cols(.default="c"), col_names=FALSE)
freq_cases <- read_delim(FREQ_CASES_FILE, delim = "\t", col_types = cols(.default="c"), col_names=TRUE)
freq_controls <- read_delim(FREQ_CONTROLS_FILE, delim = "\t", col_types = cols(.default="c"), col_names=TRUE)
adsp_intersection <- read_delim(ADSP_INTERSECTION_FILE, delim = "\t", col_types = cols(.default="c"), col_names=TRUE)

names(diffdepth)[names(diffdepth)=="variant_id"] <- "SNP"
diffmiss <- diffmiss[,c("SNP", "P")]
names(diffmiss)[names(diffmiss)=="P"] <- "p_diffmiss"
hwe_all <- subset(hwe, TEST == "ALL")
hwe_cases <- subset(hwe, TEST == "AFF")
hwe_controls <- subset(hwe, TEST == "UNAFF")
hwe_all <- hwe_all[,c("SNP", "P")]
hwe_cases <- hwe_cases[,c("SNP", "P")]
hwe_controls <- hwe_controls[,c("SNP", "P")]
names(hwe_all)[names(hwe_all)=="P"] <- "p_hwe_all"
names(hwe_cases)[names(hwe_cases)=="P"] <- "p_hwe_cases"
names(hwe_controls)[names(hwe_controls)=="P"] <- "p_hwe_controls"
freq_counts <- freq_counts[,c("SNP", "C1", "C2")]
missrate$callrate_total <- 1-as.numeric(missrate$F_MISS)
missrate <- missrate[,c("SNP", "callrate_total")]
missrate_cases$callrate_cases <- 1-as.numeric(missrate_cases$F_MISS)
missrate_cases <- missrate_cases[,c("SNP", "callrate_cases")]
missrate_controls$callrate_controls <- 1-as.numeric(missrate_controls$F_MISS)
missrate_controls <- missrate_controls[,c("SNP", "callrate_controls")]
names(multiallelic) <-"SNP"
multiallelic$multiallelic <- TRUE
names(star_alleles) <-"SNP"
star_alleles$star_alleles <- TRUE
freq_cases <- freq_cases[,c("SNP", "MAF")]
names(freq_cases)[names(freq_cases)=="MAF"] <- "maf_cases"
freq_controls <- freq_controls[,c("SNP", "MAF")]
names(freq_controls)[names(freq_controls)=="MAF"] <- "maf_controls"
adsp_intersection <- adsp_intersection[,c("SNP")]
adsp_intersection$adsp_intersection <- TRUE

qc_masterlist <- hwe_all %>% left_join(hwe_cases, by = "SNP") %>%
  left_join(hwe_controls, by = "SNP") %>%
#  left_join(diffdepth, by = "SNP") %>%
  left_join(diffmiss, by = "SNP") %>%
  left_join(freq_counts, by = "SNP") %>%
  left_join(missrate, by = "SNP") %>%
  left_join(missrate_cases, by = "SNP") %>%
  left_join(missrate_controls, by = "SNP") %>%
  left_join(multiallelic, by = "SNP") %>%
  left_join(star_alleles, by = "SNP") %>%
  left_join(freq_cases, by = "SNP") %>%
  left_join(freq_controls, by = "SNP") %>%
  left_join(adsp_intersection, by = "SNP")
  
qc_masterlist$star_alleles[is.na(qc_masterlist$star_alleles)] <- FALSE
qc_masterlist$multiallelic[is.na(qc_masterlist$multiallelic)] <- FALSE
qc_masterlist$adsp_intersection[is.na(qc_masterlist$adsp_intersection)] <- FALSE

qc_masterlist_depth <- qc_masterlist

for(i in seq_along(list_of_substudies)){
  substudy <- list_of_substudies[i]
  depth_path <- str_c(FILENAME_STEM_PERSUBSTUDY_DEPTHS, substudy, ".tab")
  temp_depth <- read_delim(depth_path, delim="\t", col_types=cols(.default="c"))
  names(temp_depth) <- str_c(substudy, "_", names(temp_depth))
  names(temp_depth)[1] <- "SNP"
  qc_masterlist_depth <- left_join(qc_masterlist_depth, temp_depth, by = "SNP")
}

write.table(qc_masterlist_depth, OUTPUT_FILENAME_QC_MASTERLIST_DEPTHS, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

qc_masterlist_substudies <- qc_masterlist_depth

for(i in seq_along(list_of_substudies)){
  substudy <- list_of_substudies[i]
  temp_diffmiss <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_diffmiss_merged.missing.tab"), delim = "\t", col_types = cols(.default="c"))
  temp_hwe <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_diffmiss_merged.hwe.tab"), delim = "\t", col_types = cols(.default="c"))
  temp_missrate <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_merged.lmiss.tab"), delim = "\t", col_types = cols(.default="c"))
  temp_missrate_cases <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_merged_cases.lmiss.tab"), delim = "\t", col_types = cols(.default="c"))
  temp_missrate_controls <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_merged_controls.lmiss.tab"), delim = "\t", col_types = cols(.default="c"))
  temp_freq_cases <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_merged_cases.frq.tab"), delim = "\t", col_types = cols(.default="c"), col_names=TRUE)
  temp_freq_controls <- read_delim(str_c(PLINK_PERSUBSTUDY_STEM, substudy, "_merged_controls.frq.tab"), delim = "\t", col_types = cols(.default="c"), col_names=TRUE)
  
  temp_diffmiss <- temp_diffmiss[,c("SNP", "P")]
  names(temp_diffmiss)[names(temp_diffmiss)=="P"] <- "p_diffmiss"
  temp_hwe_all <- subset(temp_hwe, TEST == "ALL")
  temp_hwe_cases <- subset(temp_hwe, TEST == "AFF")
  temp_hwe_controls <- subset(temp_hwe, TEST == "UNAFF")
  temp_hwe_all <- temp_hwe_all[,c("SNP", "P")]
  temp_hwe_cases <- temp_hwe_cases[,c("SNP", "P")]
  temp_hwe_controls <- temp_hwe_controls[,c("SNP", "P")]
  names(temp_hwe_all)[names(temp_hwe_all)=="P"] <- "p_hwe_all"
  names(temp_hwe_cases)[names(temp_hwe_cases)=="P"] <- "p_hwe_cases"
  names(temp_hwe_controls)[names(temp_hwe_controls)=="P"] <- "p_hwe_controls"
  temp_missrate$callrate_total <- 1-as.numeric(temp_missrate$F_MISS)
  temp_missrate <- temp_missrate[,c("SNP", "callrate_total")]
  temp_missrate_cases$callrate_cases <- 1-as.numeric(temp_missrate_cases$F_MISS)
  temp_missrate_cases <- temp_missrate_cases[,c("SNP", "callrate_cases")]
  temp_missrate_controls$callrate_controls <- 1-as.numeric(temp_missrate_controls$F_MISS)
  temp_missrate_controls <- temp_missrate_controls[,c("SNP", "callrate_controls")]
  temp_freq_cases <- temp_freq_cases[,c("SNP", "MAF")]
  names(temp_freq_cases)[names(temp_freq_cases)=="MAF"] <- "maf_cases"
  temp_freq_controls <- temp_freq_controls[,c("SNP", "MAF")]
  names(temp_freq_controls)[names(temp_freq_controls)=="MAF"] <- "maf_controls"
  names(temp_diffmiss)[names(temp_diffmiss) != "SNP"] <- str_c(substudy, "_", names(temp_diffmiss)[names(temp_diffmiss) != "SNP"])
  names(temp_hwe_all)[names(temp_hwe_all) != "SNP"] <- str_c(substudy, "_", names(temp_hwe_all)[names(temp_hwe_all) != "SNP"])
  names(temp_hwe_cases)[names(temp_hwe_cases) != "SNP"] <- str_c(substudy, "_", names(temp_hwe_cases)[names(temp_hwe_cases) != "SNP"])
  names(temp_hwe_controls)[names(temp_hwe_controls) != "SNP"] <- str_c(substudy, "_", names(temp_hwe_controls)[names(temp_hwe_controls) != "SNP"])
  names(temp_missrate)[names(temp_missrate) != "SNP"] <- str_c(substudy, "_", names(temp_missrate)[names(temp_missrate) != "SNP"])
  names(temp_missrate_cases)[names(temp_missrate_cases) != "SNP"] <- str_c(substudy, "_", names(temp_missrate_cases)[names(temp_missrate_cases) != "SNP"])
  names(temp_missrate_controls)[names(temp_missrate_controls) != "SNP"] <- str_c(substudy, "_", names(temp_missrate_controls)[names(temp_missrate_controls) != "SNP"])
  names(temp_freq_cases)[names(temp_freq_cases) != "SNP"] <- str_c(substudy, "_", names(temp_freq_cases)[names(temp_freq_cases) != "SNP"])
  names(temp_freq_controls)[names(temp_freq_controls) != "SNP"] <- str_c(substudy, "_", names(temp_freq_controls)[names(temp_freq_controls) != "SNP"])
  
  qc_masterlist_substudies <- qc_masterlist_substudies %>% left_join(temp_diffmiss, by = "SNP") %>%
    left_join(temp_hwe_all, by = "SNP") %>%
    left_join(temp_hwe_cases, by = "SNP") %>%
    left_join(temp_hwe_controls, by = "SNP") %>%
    left_join(temp_missrate, by = "SNP") %>%
    left_join(temp_missrate_cases, by = "SNP") %>%
    left_join(temp_missrate_controls, by = "SNP") %>%
    left_join(temp_freq_cases, by = "SNP") %>%
    left_join(temp_freq_controls, by = "SNP")
}  

write.table(qc_masterlist_substudies, OUTPUT_FILENAME_QC_MASTERLIST_SUBSTUDIES, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

samplesizes$column_name <- samplesizes$study %>% str_replace("case", "mean_depth_cases") %>% str_replace("control", "mean_depth_controls")
samplesizes$column_name_callrate <- samplesizes$study %>% str_replace("case", "callrate_cases") %>% str_replace("control", "callrate_controls")
for(i in 1:dim(samplesizes)[1]){
  samplesizes$column_number[i] <- which(names(qc_masterlist_substudies) == samplesizes$column_name[i])
  samplesizes$column_number_callrate[i] <- which(names(qc_masterlist_substudies) == samplesizes$column_name_callrate[i])
}
samplesizes_cases <- subset(samplesizes, str_detect(study, "case"))
samplesizes_controls <- subset(samplesizes, str_detect(study, "control"))

calculate_overall_mean_depth_cases <- function(data){
  substudy_depths_cases <- as.numeric(data[samplesizes_cases$column_number])
  substudy_callrates_cases <- as.numeric(data[samplesizes_cases$column_number_callrate])
  samplesize_variant_cases <- as.integer((substudy_callrates_cases*2*samplesizes_cases$n)/2)
  multiplied_depths_cases <- substudy_depths_cases*samplesize_variant_cases
  summed_depths_cases <- sum(multiplied_depths_cases, na.rm=TRUE)
  mean_depth_cases <- summed_depths_cases/sum(samplesize_variant_cases, na.rm=TRUE)
  return(mean_depth_cases)
}

calculate_overall_mean_depth_controls <- function(data){
  substudy_depths_controls <- as.numeric(data[samplesizes_controls$column_number])
  substudy_callrates_controls <- as.numeric(data[samplesizes_controls$column_number_callrate])
  samplesize_variant_controls <- as.integer((substudy_callrates_controls*2*samplesizes_controls$n)/2)
  multiplied_depths_controls <- substudy_depths_controls*samplesize_variant_controls
  summed_depths_controls <- sum(multiplied_depths_controls, na.rm=TRUE)
  mean_depth_controls <- summed_depths_controls/sum(samplesize_variant_controls, na.rm=TRUE)
  return(mean_depth_controls)
}

mean_depths_cases <- apply(qc_masterlist_substudies, 1, calculate_overall_mean_depth_cases)
mean_depths_controls <- apply(qc_masterlist_substudies, 1, calculate_overall_mean_depth_controls)

qc_masterlist_substudies$mean_depth_cases <- mean_depths_cases
qc_masterlist_substudies$mean_depth_controls <- mean_depths_controls

qc_names <- names(qc_masterlist_substudies)
mean_depth_cases_cols <- qc_names[str_detect(qc_names, "_mean_depth_cases")]
mean_depth_controls_cols <- qc_names[str_detect(qc_names, "_mean_depth_controls")]
callrate_controls_cols <- qc_names[str_detect(qc_names, "_callrate_controls")]
callrate_cases_cols <- qc_names[str_detect(qc_names, "_callrate_cases")]
hwe_controls_cols <- qc_names[str_detect(qc_names, "_p_hwe_controls")]
hwe_cases_cols <- qc_names[str_detect(qc_names, "_p_hwe_cases")]
hwe_all_cols <- qc_names[str_detect(qc_names, "_p_hwe_all")]
diffdepth_ttest_cols <- qc_names[str_detect(qc_names, "_p_diffdepth$")]
diffdepth_wilcox_cols <- qc_names[str_detect(qc_names, "_p_diffdepth_wilcox")]
diffmiss_cols <- qc_names[str_detect(qc_names, "_p_diffmiss")]

get_mins <- function(X){
  min_hwe_controls <- min(as.numeric(as.character(X[hwe_controls_cols])))
  min_hwe_cases <- min(as.numeric(as.character(X[hwe_cases_cols])))
  min_hwe_all <- min(as.numeric(as.character(X[hwe_all_cols])))
  min_diffdepth_ttest <- min(as.numeric(as.character(X[diffdepth_ttest_cols])))
  min_diffdepth_wilcox <- min(as.numeric(as.character(X[diffdepth_wilcox_cols])))
  min_diffmiss <- min(as.numeric(as.character(X[diffmiss_cols])))
  min_depth_cases <- min(as.numeric(X[mean_depth_cases_cols]))
  min_depth_controls <- min(as.numeric(X[mean_depth_controls_cols]))
  min_callrate_cases <- min(as.numeric(X[callrate_cases_cols]))
  min_callrate_controls <- min(as.numeric(X[callrate_controls_cols]))
  variant_id <- X["SNP"]
  collect_output <- c(variant_id, min_hwe_controls, min_hwe_cases, min_hwe_all, min_diffdepth_ttest, min_diffdepth_wilcox, min_diffmiss, min_depth_cases, min_depth_controls, min_callrate_cases, min_callrate_controls)
  return(collect_output)
}


mins <- as.data.frame(t(apply(qc_masterlist_substudies, 1, get_mins)))
names(mins) <- c("SNP", "min_hwe_controls_persubstudy", "min_hwe_cases_persubstudy", "min_hwe_all_persubstudy", "min_diffdepth_ttest_persubstudy", "min_diffdepth_wilcox_persubstudy", "min_diffmiss_persubstudy", "min_depth_cases_persubstudy", "min_depth_controls_persubstudy", "min_callrate_cases_persubstudy", "min_callrate_controls_persubstudy")
qc_masterlist_substudies <- left_join(qc_masterlist_substudies, mins, by = "SNP")

keep_for_covariates <- subset(qc_masterlist_substudies, adsp_intersection == TRUE &
                                star_alleles == FALSE &
                                multiallelic == FALSE &
                                as.numeric(as.character(callrate_cases)) >= 0.85 &
                                as.numeric(as.character(callrate_controls)) >= 0.85 &
                                as.numeric(as.character(mean_depth_cases)) >= 15 &
                                as.numeric(as.character(mean_depth_controls)) >= 15 &
                                as.numeric(as.character(p_diffmiss)) >= 5e-25 &
                                as.numeric(as.character(C1)) != 0 &
                                as.numeric(as.character(C2)) != 0 & 
                                as.numeric(as.character(min_callrate_cases_persubstudy)) >= 0.6 &
                                as.numeric(as.character(min_callrate_controls_persubstudy)) >= 0.6
)

keep_for_pca <- subset(qc_masterlist_substudies, adsp_intersection == TRUE &
                                star_alleles == FALSE &
                                multiallelic == FALSE &
                                as.numeric(as.character(callrate_cases)) >= 0.85 &
                                as.numeric(as.character(callrate_controls)) >= 0.85 &
                                as.numeric(as.character(mean_depth_cases)) >= 15 &
                                as.numeric(as.character(mean_depth_controls)) >= 15 &
                                as.numeric(as.character(p_diffmiss)) >= 5e-25 &
                                as.numeric(as.character(C1)) != 0 &
                                as.numeric(as.character(C2)) != 0 & 
                                as.numeric(as.character(min_callrate_cases_persubstudy)) >= 0.6 &
                                as.numeric(as.character(min_callrate_controls_persubstudy)) >= 0.6
                                # as.numeric(as.character(p_diffdepth_wilcox)) >= 5e-25
                         )

write.table(qc_masterlist_substudies, OUTPUT_FILENAME_QC_MASTERLIST, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(keep_for_covariates$SNP, OUTPUT_FILENAME_COVARIATES, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(keep_for_pca$SNP, OUTPUT_FILENAME_PCA, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
