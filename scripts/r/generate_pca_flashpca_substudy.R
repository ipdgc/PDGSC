library(tidyverse)
library(flashpcaR)
library(dummies)
options(scipen=999)
covars_input <- Sys.getenv("COVARS_INPUT")
samplelist_input <- Sys.getenv("SAMPLELIST")
casecontrol_input <- Sys.getenv("CASECONTROL_INPUT")
fam_input <- Sys.getenv("FAM_INPUT")
output_filename <- Sys.getenv("OUTPUT_FILENAME")
output_filename_eigenvalues <- Sys.getenv("OUTPUT_FILENAME_EIGENVALUES")
covars <- read_delim(covars_input, "\t", escape_double = FALSE, col_types = cols(.default="c")) 
casecontrol <- read_delim(casecontrol_input, "\t", escape_double = FALSE, col_types = cols(.default="c"))
samplelist <- read_delim(samplelist_input, "\t", escape_double = FALSE, col_types = cols(.default="c"), col_names=FALSE)
names(samplelist) <- "iid"
casecontrol <- inner_join(casecontrol, samplelist, by ="iid")
casecontrol$mean_depth <- NULL
casecontrol$n_sites <- NULL
fam <- read_delim(str_c(fam_input, ".fam"), delim=" ", col_names=FALSE)
pca <- flashpca(fam_input, ndim = 20, stand = "binom2")
pca_vectors <- as.data.frame(pca$vectors)
pca_values <- as.data.frame(pca$values)
names(pca_vectors) <- str_c("PC", as.character(seq(1:20)))
fam$X2 <- NULL
fam$X3 <- NULL
fam$X4 <- NULL
fam$X5 <- NULL
fam$X6 <- NULL
names(fam)[1] <- "iid"
pca_results <- cbind(fam, pca_vectors)
casecontrol$STUDY_corrected_DZNE <- NULL 
casecontrol$STUDY_corrected_IPDGC <- NULL 
casecontrol$STUDY_corrected_PPMI <- NULL 
casecontrol$STUDY_corrected_umiami <- NULL 
casecontrol$STUDY_corrected_ADSP <- NULL 
casecontrol$STUDY_corrected_Baylor <- NULL 
casecontrol$STUDY_corrected_cambridge <- NULL 
casecontrol$STUDY_corrected_Luebeck <- NULL 
casecontrol$STUDY_corrected_majamaa <- NULL 
casecontrol$STUDY_corrected_oslo <- NULL 
casecontrol$STUDY_corrected_oxford <- NULL 
casecontrol$STUDY_corrected_tuebingen <- NULL 
casecontrol$STUDY_corrected_ucl <- NULL 
casecontrol$STUDY_corrected_UCL <- NULL
casecontrol_pc_cols <- which(names(casecontrol)%in%c(str_c("PC", as.character(seq(1:20)))))
casecontrol_nopc <- casecontrol[,-casecontrol_pc_cols]
casecontrol_newpc <- left_join(casecontrol_nopc, pca_results, by = "iid")
ped <- left_join(casecontrol_newpc, covars, by = "iid")
ped <- ped %>% separate(GROUP_matched_ucl, into = c("STUDY_matched", "status"), sep = "_", remove = FALSE)
ped_dummies <- dummy(ped$STUDY_matched, sep = "_", drop = TRUE, fun = as.integer, verbose = FALSE)
ped_names_dummies <- colnames(ped_dummies)
ped_dummies_df <- as.data.frame(ped_dummies)
names(ped_dummies_df) <- ped_names_dummies
ped_dummies_df_merged <- bind_cols(ped, ped_dummies_df)
write.table(ped_dummies_df_merged, output_filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(pca_values, output_filename_eigenvalues, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
