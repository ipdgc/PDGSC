library(tidyverse)
library(stringr)
setwd("clustering_results_keepmono_50depth_50anc_2to1")
ped <- read_delim("caseControlAnalysisData_matched_group2_2to1_depth_ancestry_210818.ped", delim="\t", col_types=cols(.default="c"))
ped <- ped %>% separate(GROUP_matched_ucl, c("STUDY_matched", "case_matched"), "_")
studies <- unique(ped$STUDY_matched)
cases <- subset(ped, case2_cont1 == "2")
controls <- subset(ped, case2_cont1 == "1")
samplelist_matched <- ped$iid
caselist_matched <- cases$iid
controllist_matched <- controls$iid
phenofile_matched <- ped[,c("fid", "iid", "case2_cont1")]

write.table(samplelist_matched, "samplelist_matched_group2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(caselist_matched, "caselist_matched_group2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(controllist_matched, "controllist_matched_group2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(phenofile_matched, "phenofile_matched_group2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(studies, "list_of_substudies.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

for(i in seq_along(studies)){
  study <- studies[i]
  ped_study <- subset(ped, STUDY_matched == study)
  cases_study <- subset(cases, STUDY_matched == study)
  controls_study <- subset(controls, STUDY_matched == study)
  samplelist_study <- ped_study$iid
  caselist_study <- cases_study$iid
  controllist_study <- controls_study$iid
  phenofile_study <- ped_study[,c("fid", "iid", "case2_cont1")]

  write.table(samplelist_study, str_c("samplelist_matched_", study, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(caselist_study, str_c("caselist_matched_", study, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(controllist_study, str_c("controllist_matched_", study, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(phenofile_study, str_c("phenofile_matched_", study, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
