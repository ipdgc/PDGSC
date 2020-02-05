library("stringr")
library("plyr")
sample_info <- Sys.getenv("SAMPLE_INFO")
fam_file <- Sys.getenv("FAM_FILE")
updated_fam_file <- Sys.getenv("UPDATED_FAM_FILE")
missing_clinical_data <- Sys.getenv("MISSING_CLINICAL")

clinical <- read.table(sample_info, sep="\t", header=TRUE, stringsAsFactors=FALSE)
fam <- read.table(fam_file, sep=" ", header=FALSE, stringsAsFactors=FALSE)
missing_pheno <- fam[!fam$V1%in%clinical$ID,]
names(fam)[2] <- "ID"
newfam <- join(fam, clinical, by="ID")
newfam$V6 <- newfam$CASE
newfam$V5 <- newfam$male1_Female2
newfam_rightcolumns <- newfam[,c(1:6,10)]
newfam_rightcolumns$V6 <- newfam_rightcolumns$V6+1
write.table(newfam_rightcolumns, updated_fam_file, row.names=FALSE, col.names=FALSE, quote=FALSE, na="0", sep="\t")
write.table(missing_pheno, missing_clinical_data, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
