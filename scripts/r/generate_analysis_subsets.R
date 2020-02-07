library(tidyverse)
library(stringr)
CASECONTROLPEDFILE <- Sys.getenv("CASECONTROLPEDFILE")
casecontrol <- read_delim(CASECONTROLPEDFILE, delim="\t", col_types=cols(.default="c"))
substudies <- unique(casecontrol$STUDY_matched)
substudies_UCLonly <- c("UCL1", "UCL2")
substudies_majamaa <- "majamaa"
substudies_oxford <- "oxford"
substudies_leftovers <- c("cambridge", "Baylor", "oslo", "tuebingen")
substudies_UCL1 <- "UCL1"
substudies_UCL2 <- "UCL2"
substudies_cambridge <- "cambridge"
substudies_Baylor <- "Baylor"
substudies_oslo <- "oslo"
substudies_tuebingen <- "tuebingen"
pedfiles <- list()
list_of_groups <- c("substudies_UCLonly", "substudies_majamaa", "substudies_oxford", "substudies_leftovers", "substudies_UCL1", "substudies_UCL2", "substudies_cambridge", "substudies_Baylor", "substudies_oslo", "substudies_tuebingen")
for(i in seq_along(list_of_groups)){
    group <- list_of_groups[i]
    substudy <- eval(parse(text = group))
    substudy_colname <- str_c("STUDY_matched_", substudy)
    pedfiles[[i]] <- subset(casecontrol, STUDY_matched%in%substudy)
    substudy_colnumber <- which(names(pedfiles[[i]])%in%substudy_colname)
    keep_studycovar <- pedfiles[[i]][,substudy_colnumber]
    remove_studycovar <- which(str_detect(names(pedfiles[[i]]), "STUDY_matched_"))
    pedfiles[[i]] <- pedfiles[[i]][,-remove_studycovar]
    pedfiles[[i]] <- bind_cols(pedfiles[[i]], keep_studycovar)
    names(pedfiles)[i] <- str_c("caseControlAnalysisData040918_group2_matched_2to1_depth_anc_", group, ".ped")
    write.table(pedfiles[[i]], names(pedfiles)[i], sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}
write.table(list_of_groups, "list_of_casecontrol_subsets.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)