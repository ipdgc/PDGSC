### TI/TV ratio
setwd("c:/Users/dkia/Documents/PDGSC_QC/reqc/exclusions/")
library(readr)
pdgsc_postqc_210817 <- read_delim("~/PDGSC_QC/pdgsc_postqc_210817.istats", "\t", escape_double = FALSE, trim_ws = TRUE)
pdgsc_allsamples_preqc <- read_delim("~/PDGSC_QC/meta/meta_qc_metrics/pdgsc_allsamples_preqc.istats", "\t", escape_double = FALSE, trim_ws = TRUE)
summary(pdgsc_allsamples_preqc$TITV)
sd(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)
mean(pdgsc_allsamples_preqc$TITV)+3*sd(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)
mean(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)+3*sd(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)
mean(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)-3*sd(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)

summary(pdgsc_postqc_210817$TITV)
sd(pdgsc_postqc_210817$TITV, na.rm = TRUE)
mean(pdgsc_postqc_210817$TITV)+3*sd(pdgsc_postqc_210817$TITV, na.rm = TRUE)
mean(pdgsc_postqc_210817$TITV, na.rm = TRUE)+3*sd(pdgsc_postqc_210817$TITV, na.rm = TRUE)
mean(pdgsc_postqc_210817$TITV, na.rm = TRUE)-3*sd(pdgsc_postqc_210817$TITV, na.rm = TRUE)


### DP_ALT FILTERING ###
dp_alt_exclude <- subset(pdgsc_postqc_210817, DP/NALT < 25)
dp_alt_exclude_preqc <- subset(pdgsc_allsamples_preqc, DP/NALT < 25)
write.table(dp_alt_exclude, "C:/Users/dkia/Documents/PDGSC_QC/meta/alt_dp_exclusions_200917.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE)
merck_dp_alt_exclude <- subset(Merck_QC_iStats, DP/NALT < 25)
write.table(merck_dp_alt_exclude, "C:/Users/dkia/Documents/PDGSC_QC/meta/merck_alt_dp_exclusions_200917.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE)

### NOW HAVE A LOOK BACK AT THE QC
checkhet_remove_hard <- read_delim("~/PDGSC_QC/reqc/exclusions/checkhet_remove_hard.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
sexproblem <- read_delim("~/PDGSC_QC/reqc/exclusions/sexproblem.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
ADSP_freemix_0_05 <- read_delim("~/PDGSC_QC/meta/rareMETALS2_working_scripts/ADSP_freemix_0.05.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ibd_exclusion_postplinkrel <- read_delim("~/PDGSC_QC/reqc/exclusions/ibd_exclusion_postplinkrel.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
istats_exclusions <- read_delim("~/PDGSC_QC/reqc/exclusions/istats_exclusions.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
pca_exclusion_idlist <- read_delim("~/PDGSC_QC/reqc/exclusions/pca_exclusion_idlist.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
underage <- read_delim("~/PDGSC_QC/reqc/exclusions/underage.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
qc_mismatch_mikes_reasons <- read_delim("~/PDGSC_QC/meta/qc_mismatch_mikes_reasons.txt", "\t", col_types="cc", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
qc_mismatch_reasons <- read_delim("~/PDGSC_QC/meta/qc_mismatch_reasons.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
checkhet_remove_hard$X1 <- NULL
het_bringback <- subset(checkhet_remove_hard, X7 <= 0.1 & X7 >= -0.1)
sex_inspect <- subset(sexproblem, X7 <= 0.65 & X7 >= 0.35)
sex_inspect$X5 <- ifelse(test = sex_inspect$X7 > 0.5, 1, 2)
sex_bringback <- subset(sex_inspect, X4==X5)

het_bringback$X1 <- NULL
sex_bringback$X1 <- NULL

head(sex_bringback)
head(het_bringback)

sex_bringback[sex_bringback$X2%in%ADSP_freemix_0_05$IID,]
het_bringback[het_bringback$X2%in%ADSP_freemix_0_05$IID,]

sex_bringback_freemix <- sex_bringback[!sex_bringback$X2%in%ADSP_freemix_0_05$IID,]
het_bringback_freemix <- het_bringback[!het_bringback$X2%in%ADSP_freemix_0_05$IID,]

ibd_exclusion_postplinkrel[ibd_exclusion_postplinkrel$X1%in%ADSP_freemix_0_05$IID,]
ibd_bringback_freemix <- ibd_exclusion_postplinkrel[!ibd_exclusion_postplinkrel$X1%in%ADSP_freemix_0_05$IID,]

sex_bringback[sex_bringback$X2%in%istats_exclusions$X1,]
het_bringback[het_bringback$X2%in%istats_exclusions$X1,]
ibd_bringback_freemix[ibd_bringback_freemix$X1%in%istats_exclusions$X1,]

sex_bringback[sex_bringback$X2%in%pca_exclusion_idlist$X2,]
het_bringback[het_bringback$X2%in%pca_exclusion_idlist$X2,]
ibd_bringback_freemix[ibd_bringback_freemix$X1%in%pca_exclusion_idlist$X2,]

sex_bringback[sex_bringback$X2%in%underage$X1,]
het_bringback[het_bringback$X2%in%underage$X1,]
ibd_bringback_freemix[ibd_bringback_freemix$X1%in%underage$X1,]
het_bringback_freemix_age <- het_bringback_freemix[!het_bringback_freemix$X2%in%underage$X1,]
het_bringback_freemix_age_pca <- het_bringback_freemix_age[!het_bringback_freemix_age$X2%in%pca_exclusion_idlist$X2,]

sex_bringback[sex_bringback$X2%in%qc_mismatch_mikes_reasons$X1,]
het_bringback[het_bringback$X2%in%qc_mismatch_mikes_reasons$X1,]
ibd_bringback_freemix[ibd_bringback_freemix$X1%in%qc_mismatch_mikes_reasons$X1,]

sex_bringback[sex_bringback$X2%in%qc_mismatch_reasons$X1,]
het_bringback[het_bringback$X2%in%qc_mismatch_reasons$X1,]

sex_bringback_freemix_mike <- sex_bringback_freemix[sex_bringback_freemix$X2%in%qc_mismatch_reasons$X1,]
het_bringback_freemix_age_pca[het_bringback_freemix_age_pca$X2%in%qc_mismatch_reasons$X1,]
ibd_bringback_freemix[ibd_bringback_freemix$X1%in%qc_mismatch_reasons$X1,]

bringback_all <- c(sex_bringback_freemix_mike$X2,het_bringback_freemix_age_pca$X2,ibd_bringback_freemix$X1)
bringback_unique <- unique(bringback_all)

bringback_unique[bringback_unique%in%postqc$ID]
bringback_unique <- as.data.frame(bringback_unique)
bringback_istats <- merge(bringback_unique, pdgsc_allsamples_preqc, by.x="bringback_unique", by.y = "ID")

bringback_post_altdp <- bringback_unique[!bringback_unique%in%dp_alt_exclude_preqc$ID]

postqc_titv <- subset(postqc, TITV > 3.62 | TITV < 3.19)
postqc_titv_wholesample <- subset(postqc, TITV < mean(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)-3*sd(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)|TITV > mean(pdgsc_allsamples_preqc$TITV, na.rm = TRUE)+3*sd(pdgsc_allsamples_preqc$TITV, na.rm = TRUE))
postqc_titv_postqc <- subset(postqc, TITV < mean(postqc$TITV, na.rm = TRUE)-3*sd(postqc$TITV, na.rm = TRUE)|TITV > mean(postqc$TITV, na.rm = TRUE)+3*sd(postqc$TITV, na.rm = TRUE))

names(qc_keep)[1] <- "ID"

pdgsc_postqc_210917 <- merge(qc_keep, pdgsc_allsamples_preqc, by="ID")
write.table(pdgsc_postqc_210917, "C:/Users/dkia/Documents/PDGSC_QC/meta/meta_qc_metrics_210917/pdgsc_postqc_210917.istats", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

hist(pdgsc_diffmiss_original_no_midp_missing$P)
hist(differential_depth_pvalues$P_DIFF_DEPTH)

summary(differential_depth_pvalues$P_DIFF_DEPTH)

hist(pdgsc_postqc_210917$DP/pdgsc_postqc_210917$NALT)

head(updatedSampleInfoJanuary23rd2017)

samples <- read.table("C:/Users/dkia/Documents/PDGSC_QC/updatedSampleInfoJanuary23rd2017.tab", sep="\t", header=TRUE, stringsAsFactors = FALSE)
samples$CASE <- as.factor(samples$CASE)
samples$male1_Female2 <- as.factor(samples$male1_Female2)
samples$STUDY <- as.factor(samples$STUDY)


postqc_list <- pdgsc_postqc_210917$ID
preDPALT_list <- caseControlAnalysisDataSeptember5th2017$iid
samples_pdgsc <- samples[samples$ID%in%pdgsc_allsamples_preqc$ID,]
samples_pdgsc$casecontrol <- as.factor(ifelse(samples_pdgsc$CASE==1, "CASE", "CONTROL"))
samples_postDPALT <- samples_pdgsc[samples_pdgsc$ID%in%postqc_list,]
samples_preDPALT <- samples_pdgsc[samples_pdgsc$ID%in%preDPALT_list,]

write.table(samples_postDPALT$ID, "C:/Users/dkia/Documents/PDGSC_QC/meta/pdgsc_postDPALT_samplelist.txt", sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(samples_preDPALT$ID, "C:/Users/dkia/Documents/PDGSC_QC/meta/pdgsc_preDPALT_samplelist.txt", sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)


plot(samples_preDPALT$MEAN_DEPTH ~ samples_preDPALT$casecontrol, ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH before DP_ALT exclusions")
plot(samples_postDPALT$MEAN_DEPTH ~ samples_postDPALT$casecontrol, ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH after DP_ALT exclusions")

plot(samples_preDPALT$MEAN_DEPTH[samples_preDPALT$MEAN_DEPTH < 100] ~ samples_preDPALT$casecontrol[samples_preDPALT$MEAN_DEPTH < 100], ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH before DP_ALT exclusions\nAfter excluding samples with MEAN_DEPTH > 100 (n=191)")
plot(samples_postDPALT$MEAN_DEPTH[samples_postDPALT$MEAN_DEPTH < 100] ~ samples_postDPALT$casecontrol[samples_postDPALT$MEAN_DEPTH < 100], ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH after DP_ALT exclusions\nAfter excluding samples with MEAN_DEPTH > 100 (n=190)")

plot(samples_preDPALT$MEAN_DEPTH[samples_preDPALT$MEAN_DEPTH < 50] ~ samples_preDPALT$casecontrol[samples_preDPALT$MEAN_DEPTH < 50], ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH before DP_ALT exclusions\nAfter excluding samples with MEAN_DEPTH > 50 (n=514)")
plot(samples_postDPALT$MEAN_DEPTH[samples_postDPALT$MEAN_DEPTH < 50] ~ samples_postDPALT$casecontrol[samples_postDPALT$MEAN_DEPTH < 50], ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH after DP_ALT exclusions\nAfter excluding samples with MEAN_DEPTH > 50 (n=225)")


sum(samples_preDPALT$MEAN_DEPTH>50)

highdepth <- subset(samples_preDPALT, MEAN_DEPTH>50)
lowdepth_post <- subset(samples_postDPALT, MEAN_DEPTH<50)
lowdepth_pre <- subset(samples_preDPALT, MEAN_DEPTH<50)

test_cleandepth_pre <- subset(samples_preDPALT, MEAN_DEPTH<40 & MEAN_DEPTH>25)
test_cleandepth_post <- subset(samples_postDPALT, MEAN_DEPTH<40 & MEAN_DEPTH>25)

summary(as.factor(samples_preDPALT$casecontrol))
summary(as.factor(samples_postDPALT$casecontrol))
summary(as.factor(test_cleandepth_pre$casecontrol))
summary(as.factor(test_cleandepth_post$casecontrol))
summary(test_cleandepth_pre$MEAN_DEPTH[test_cleandepth_pre$casecontrol=="CASE"])
summary(test_cleandepth_pre$MEAN_DEPTH[test_cleandepth_pre$casecontrol=="CONTROL"])
summary(test_cleandepth_post$MEAN_DEPTH[test_cleandepth_post$casecontrol=="CASE"])
summary(test_cleandepth_post$MEAN_DEPTH[test_cleandepth_post$casecontrol=="CONTROL"])


summary(samples_preDPALT$MEAN_DEPTH[samples_preDPALT$casecontrol=="CASE"])
summary(samples_preDPALT$MEAN_DEPTH[samples_preDPALT$casecontrol=="CONTROL"])
summary(samples_postDPALT$MEAN_DEPTH[samples_postDPALT$casecontrol=="CASE"])
summary(samples_postDPALT$MEAN_DEPTH[samples_postDPALT$casecontrol=="CONTROL"])


summary(as.factor(test_cleandepth_post$MEAN_DEPTH))

plot(test_cleandepth_pre$MEAN_DEPTH ~ test_cleandepth_pre$casecontrol, ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH before DP_ALT exclusions\nAfter excluding samples with MEAN_DEPTH > 40 OR < 25 (n=225)")
plot(test_cleandepth_post$MEAN_DEPTH ~ test_cleandepth_post$casecontrol, ylab = NULL, xlab = NULL, main="PDGSC MEAN_DEPTH after DP_ALT exclusions\nAfter excluding samples with MEAN_DEPTH > 40 OR < 25 (n=225)")

summary(highdepth$MEAN_DEPTH)
summary(highdepth$N_SITES)

summary(lowdepth$MEAN_DEPTH)
summary(lowdepth$N_SITES)

write.table(lowdepth_pre$ID, "predpalt_lowdepth_keep.txt", sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(lowdepth_post$ID, "postdpalt_lowdepth_keep.txt", sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)

head(test_cleandepth_post)
dim(test_cleandepth_post)
dim(test_cleandepth_pre)

write.table(test_cleandepth_post$ID, "postdpalt_meandepthfiltered_keep.txt", sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(test_cleandepth_pre$ID, "predpalt_meandepthfiltered_keep.txt", sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)
