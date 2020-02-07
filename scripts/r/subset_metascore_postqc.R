library(tidyverse)
library(stringr)
QC_MASTERLIST <- Sys.getenv("QC_MASTERLIST")
qc_masterlist <- read_delim(QC_MASTERLIST, delim="\t", col_types=cols(.default="c"), col_names=TRUE)
samplesizes <- read_delim("clustering_results_keepmono_50depth_50anc_2to1/samplesizes.txt", delim="\t")
list_of_substudies <- c("majamaa", "oxford", "UCL1", "UCL2")
samplesizes$column_name <- samplesizes$study %>% str_replace("case", "mean_depth_cases") %>% str_replace("control", "mean_depth_controls")
samplesizes$column_name_callrate <- samplesizes$study %>% str_replace("case", "callrate_cases") %>% str_replace("control", "callrate_controls")
samplesizes$study_name <- samplesizes$study %>% str_replace("_case", "") %>% str_replace("_control", "")
for(i in 1:dim(samplesizes)[1]){
  samplesizes$column_number[i] <- which(names(qc_masterlist) == samplesizes$column_name[i])
  samplesizes$column_number_callrate[i] <- which(names(qc_masterlist) == samplesizes$column_name_callrate[i])
}
samplesizes_cases <- subset(samplesizes, str_detect(study, "case"))
samplesizes_controls <- subset(samplesizes, str_detect(study, "control"))
calculate_overall_mean_depth_cases_dropstudy <- function(data){
  substudy_depths_cases <- as.numeric(data[samplesizes_cases_studies_to_keep$column_number])
  substudy_callrates_cases <- as.numeric(data[samplesizes_cases_studies_to_keep$column_number_callrate])
  samplesize_variant_cases <- as.integer((substudy_callrates_cases*2*samplesizes_cases_studies_to_keep$n)/2)
  multiplied_depths_cases <- substudy_depths_cases*samplesize_variant_cases
  summed_depths_cases <- sum(multiplied_depths_cases, na.rm=TRUE)
  mean_depth_cases <- summed_depths_cases/sum(samplesize_variant_cases, na.rm=TRUE)
  return(mean_depth_cases)
}
calculate_overall_mean_depth_controls_dropstudy <- function(data){
  substudy_depths_controls <- as.numeric(data[samplesizes_controls_studies_to_keep$column_number])
  substudy_callrates_controls <- as.numeric(data[samplesizes_controls_studies_to_keep$column_number_callrate])
  samplesize_variant_controls <- as.integer((substudy_callrates_controls*2*samplesizes_controls_studies_to_keep$n)/2)
  multiplied_depths_controls <- substudy_depths_controls*samplesize_variant_controls
  summed_depths_controls <- sum(multiplied_depths_controls, na.rm=TRUE)
  mean_depth_controls <- summed_depths_controls/sum(samplesize_variant_controls, na.rm=TRUE)
  return(mean_depth_controls)
}
get_min_callrates_cases_dropstudy <- function(data){
  substudy_callrates_cases <- as.numeric(data[samplesizes_cases_studies_to_keep$column_number_callrate])
  min_callrate_cases <- min(substudy_callrates_cases)
  return(min_callrate_cases)
}
get_min_callrates_controls_dropstudy <- function(data){
  substudy_callrates_controls <- as.numeric(data[samplesizes_controls_studies_to_keep$column_number_callrate])
  min_callrate_controls <- min(substudy_callrates_controls)
  return(min_callrate_controls)
}
substudies_UCLonly <- c("UCL1", "UCL2")
substudies_majamaa <- "majamaa"
substudies_oxford <- "oxford"
substudies_leftovers <- c("cambridge", "Baylor", "oslo", "tuebingen")
substudies_UCL1 <- "UCL1"
substudies_UCL2 <- "UCL2"
list_of_groups <- c("substudies_UCLonly", "substudies_leftovers")
collected_mean_depths_dropstudies <- NULL
collected_min_callrates_dropstudies <- NULL
for(i in seq_along(list_of_groups)){
  group <- list_of_groups[i]
  substudies_to_keep <- eval(parse(text = group))
  samplesizes_cases_studies_to_keep <- subset(samplesizes_cases, study_name%in%substudies_to_keep)
  samplesizes_controls_studies_to_keep <- subset(samplesizes_controls, study_name%in%substudies_to_keep)
  mean_depths_cases_dropstudy <- apply(qc_masterlist, 1, calculate_overall_mean_depth_cases_dropstudy)
  print(str_c("Calculated mean depths for cases in study ", i, "!"))
  mean_depths_controls_dropstudy <- apply(qc_masterlist, 1, calculate_overall_mean_depth_controls_dropstudy)
  print(str_c("Calculated mean depths for controls in study ", i, "!"))
  mean_depths_cases_dropstudy_df <- as.data.frame(mean_depths_cases_dropstudy)
  mean_depths_controls_dropstudy_df <- as.data.frame(mean_depths_controls_dropstudy)
  names(mean_depths_cases_dropstudy_df) <- str_c("mean_depth_cases_", group)
  names(mean_depths_controls_dropstudy_df) <- str_c("mean_depth_controls_", group)
  collected_mean_depths_dropstudies <- bind_cols(collected_mean_depths_dropstudies, mean_depths_cases_dropstudy_df, mean_depths_controls_dropstudy_df)
  min_callrate_cases_dropstudy <- apply(qc_masterlist, 1, get_min_callrates_cases_dropstudy)
  print(str_c("Calculated min callrate for cases in study ", i, "!"))
  min_callrate_controls_dropstudy <- apply(qc_masterlist, 1, get_min_callrates_controls_dropstudy)
  print(str_c("Calculated min callrate for controls in study ", i, "!"))
  min_callrate_cases_dropstudy_df <- as.data.frame(min_callrate_cases_dropstudy)
  min_callrate_controls_dropstudy_df <- as.data.frame(min_callrate_controls_dropstudy)
  names(min_callrate_cases_dropstudy_df) <- str_c("min_callrate_cases_", group)
  names(min_callrate_controls_dropstudy_df) <- str_c("min_callrate_controls_", group)
  collected_min_callrates_dropstudies <- bind_cols(collected_min_callrates_dropstudies, min_callrate_cases_dropstudy_df, min_callrate_controls_dropstudy_df)
  print(str_c("Finished study ", i, "!"))
}
write.table(collected_mean_depths_dropstudies, "clustering_results_keepmono_50depth_50anc_2to1/collected_mean_depths_grouped.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(collected_min_callrates_dropstudies, "clustering_results_keepmono_50depth_50anc_2to1/collected_min_callrates_grouped.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

qc_masterlist_mean_depths_min_callrates_dropstudy <- qc_masterlist %>% bind_cols(collected_mean_depths_dropstudies) %>% bind_cols(collected_min_callrates_dropstudies)

dodgy_variants_substudies_UCLonly <- subset(qc_masterlist_mean_depths_min_callrates_dropstudy, (as.numeric(UCL1_maf_cases) > 0.01 & as.numeric(UCL1_maf_controls) == 0) |
							(as.numeric(UCL1_maf_controls) > 0.01 & as.numeric(UCL1_maf_cases) == 0) |
							(as.numeric(UCL1_maf_cases) < 0.99 & as.numeric(UCL1_maf_controls) == 1) |
							(as.numeric(UCL1_maf_controls) < 0.99 & as.numeric(UCL1_maf_cases) == 1) |
							(as.numeric(UCL2_maf_cases) > 0.01 & as.numeric(UCL2_maf_controls) == 0) |
							(as.numeric(UCL2_maf_controls) > 0.01 & as.numeric(UCL2_maf_cases) == 0) |
							(as.numeric(UCL2_maf_cases) < 0.99 & as.numeric(UCL2_maf_controls) == 1) |
							(as.numeric(UCL2_maf_controls) < 0.99 & as.numeric(UCL2_maf_cases) == 1))$SNP

dodgy_variants_substudies_leftovers <- subset(qc_masterlist_mean_depths_min_callrates_dropstudy, (as.numeric(cambridge_maf_cases) > 0.01 & as.numeric(cambridge_maf_controls) == 0) |
							(as.numeric(cambridge_maf_controls) > 0.01 & as.numeric(cambridge_maf_cases) == 0) |
							(as.numeric(cambridge_maf_cases) < 0.99 & as.numeric(cambridge_maf_controls) == 1) |
							(as.numeric(cambridge_maf_controls) < 0.99 & as.numeric(cambridge_maf_cases) == 1) |
							(as.numeric(Baylor_maf_cases) > 0.01 & as.numeric(Baylor_maf_controls) == 0) |
							(as.numeric(Baylor_maf_controls) > 0.01 & as.numeric(Baylor_maf_cases) == 0) |
							(as.numeric(Baylor_maf_cases) < 0.99 & as.numeric(Baylor_maf_controls) == 1) |
							(as.numeric(Baylor_maf_controls) < 0.99 & as.numeric(Baylor_maf_cases) == 1) |
							(as.numeric(oslo_maf_cases) > 0.01 & as.numeric(oslo_maf_controls) == 0) |
							(as.numeric(oslo_maf_controls) > 0.01 & as.numeric(oslo_maf_cases) == 0) |
							(as.numeric(oslo_maf_cases) < 0.99 & as.numeric(oslo_maf_controls) == 1) |
							(as.numeric(oslo_maf_controls) < 0.99 & as.numeric(oslo_maf_cases) == 1) |
							(as.numeric(tuebingen_maf_cases) > 0.01 & as.numeric(tuebingen_maf_controls) == 0) |
							(as.numeric(tuebingen_maf_controls) > 0.01 & as.numeric(tuebingen_maf_cases) == 0) |
							(as.numeric(tuebingen_maf_cases) < 0.99 & as.numeric(tuebingen_maf_controls) == 1) |
							(as.numeric(tuebingen_maf_controls) < 0.99 & as.numeric(tuebingen_maf_cases) == 1))$SNP


for(i in seq_along(list_of_groups)){
group <- list_of_groups[i]
keep_variants <- subset(qc_masterlist_mean_depths_min_callrates_dropstudy, adsp_intersection == "TRUE" &
                                star_alleles == "FALSE" &
                                multiallelic == "FALSE" &
								!SNP%in%eval(parse(text = str_c("dodgy_variants_", group))) &
                                as.numeric(as.character(eval(parse(text = str_c("mean_depth_cases_", group))))) >= 15 &
                                as.numeric(as.character(eval(parse(text = str_c("mean_depth_controls_", group))))) >= 15 &
                                as.numeric(as.character(eval(parse(text = str_c("min_callrate_cases_", group))))) >= 0.6 &
                                as.numeric(as.character(eval(parse(text = str_c("min_callrate_controls_", group))))) >= 0.6 &
                                as.numeric(as.character(p_diffmiss)) >= 5e-25
								)
keep_variants_positions <- keep_variants$SNP
write.table(keep_variants_positions, str_c("group2_postqc_keep_variants_2to1_depth_anc_hwe5_persubstudy50_", group, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

list_of_groups_rest <- c("majamaa", "oxford", "UCL1", "UCL2")
for(i in seq_along(list_of_groups_rest)){
group <- list_of_groups_rest[i]
dodgy_variants <- subset(qc_masterlist_mean_depths_min_callrates_dropstudy, (as.numeric(as.character(eval(parse(text = str_c(group, "_maf_cases"))))) > 0.01 & as.numeric(as.character(eval(parse(text = str_c(group, "_maf_controls"))))) == 0) |
							(as.numeric(as.character(eval(parse(text = str_c(group, "_maf_controls"))))) > 0.01 & as.numeric(as.character(eval(parse(text = str_c(group, "_maf_cases"))))) == 0) |
							(as.numeric(as.character(eval(parse(text = str_c(group, "_maf_cases"))))) < 0.99 & as.numeric(as.character(eval(parse(text = str_c(group, "_maf_controls"))))) == 1) |
							(as.numeric(as.character(eval(parse(text = str_c(group, "_maf_controls"))))) < 0.99 & as.numeric(as.character(eval(parse(text = str_c(group, "_maf_cases"))))) == 1))$SNP

keep_variants <- subset(qc_masterlist_mean_depths_min_callrates_dropstudy, adsp_intersection == "TRUE" &
                                star_alleles == "FALSE" &
                                multiallelic == "FALSE" &
 								!SNP%in%dodgy_variants &
                                as.numeric(as.character(eval(parse(text = str_c(group, "_mean_depth_cases"))))) >= 15 &
                                as.numeric(as.character(eval(parse(text = str_c(group, "_mean_depth_controls"))))) >= 15 &
                                as.numeric(as.character(p_diffmiss)) >= 5e-25
								)
keep_variants_positions <- keep_variants$SNP
print(length(keep_variants_positions))
write.table(keep_variants_positions, str_c("group2_postqc_keep_variants_2to1_depth_anc_hwe5_persubstudy50_substudies_", group, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

