library(readr)
caseControlAnalysisDataApril13th2017 <- read_delim("~/PDGSC_QC/caseControlAnalysisDataApril13th2017.ped", "\t", escape_double = FALSE, col_types = c("ccccccccccccccccccccccccccccccccc"), trim_ws = TRUE)
ADSP_freemix_0_05 <- read_delim("~/PDGSC_QC/meta/rareMETALS2_working_scripts/ADSP_freemix_0.05.txt", "\t", escape_double = FALSE, col_types = cols(FREEMIX = col_character()), trim_ws = TRUE)
caseControlAnalysisDataAugust21st2017 <- subset(caseControlAnalysisDataApril13th2017, !iid%in%ADSP_freemix_0_05$IID)
caseOnlyAnalysisDataAugust21st2017 <- subset(caseControlAnalysisDataAugust21st2017, case2_cont1=="2")

write.table(caseControlAnalysisDataAugust21st2017, "C://Users/dkia/Documents/PDGSC_QC/meta/caseControlAnalysisDataAugust21st2017", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
write.table(caseOnlyAnalysisDataAugust21st2017, "C://Users/dkia/Documents/PDGSC_QC/meta/caseOnlyAnalysisDataAugust21st2017", sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

