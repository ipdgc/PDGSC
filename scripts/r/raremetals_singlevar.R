# Import libraries
# library(rareMETALS)
library(rareMETALS2)

# meta_QC_define.R (not needed anymore)

# Import variables from the environment into the R session
anno_london_casecontrol_assoc <- Sys.getenv("anno_london_casecontrol_assoc")
anno_london_caseonly_assoc <- Sys.getenv("anno_london_caseonly_assoc")
anno_merck_casecontrol_assoc <- Sys.getenv("anno_merck_casecontrol_assoc")
anno_merck_caseonly_assoc <- Sys.getenv("anno_merck_caseonly_assoc")
london_casecontrol_covar <- Sys.getenv("london_casecontrol_covar")
london_caseonly_covar <- Sys.getenv("london_caseonly_covar")
merck_casecontrol_covar <- Sys.getenv("merck_casecontrol_covar")
merck_caseonly_covar <- Sys.getenv("merck_caseonly_covar")
CHROM <- Sys.getenv("i")
genes <- Sys.getenv("genelist")

# Make the score and cov files
cov_files_casecontrol <- c(london_casecontrol_covar, merck_casecontrol_covar)
score_stats_files_casecontrol <- c(anno_london_casecontrol_assoc, anno_merck_casecontrol_assoc)
cov_files_caseonly <- c(london_caseonly_covar, merck_caseonly_covar)
score_stats_files_caseonly <- c(anno_london_caseonly_assoc, anno_merck_caseonly_assoc)

# Make annotations for gene-based tests
AllVars_annotations <- 'gene'
Lof_annotations <- 'Frameshift|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
LofPlusIndel_annotations <- 'Insertion|Deletion|Frameshift|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
ProCo_annotations <- 'Nonsynonymous|CodonGain|CodonLoss|Frameshift|Normal_Splice_Site|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
ProCoPlusIndel_annotations <- 'Insertion|Deletion|Nonsynonymous|CodonGain|CodonLoss|Frameshift|Normal_Splice_Site|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'

genelist_df <- read.table(genes, header=FALSE, stringsAsFactors=FALSE, sep="")
genelist <- genelist_df$V1
gene <- genelist

range <- paste(CHROM, ":1-1000000000", sep="")

# Single variant meta
res_single_casecontrol <- rareMETALS2.single(score_stats_files_casecontrol,range=range,alternative="two.sided",ix.gold=1,callrate.cutoff=0,hwe.cutoff=0,hwe.ctrl.cutoff=0)
res_single_casecontrol$integratedData <- NULL
res_single_casecontrol$raw.data <- NULL
res_single_casecontrol$clean.data <- NULL
res_single_casecontrol_df <- do.call(cbind, lapply(res_single_casecontrol, data.frame, stringsAsFactors=FALSE))
names(res_single_casecontrol_df) <- c("p.value", "ref", "alt", "statistic", "direction.by.study", "anno", "af", "afCase", "afCtrl", "QC.by.study", "no.sample", "no.case", "no.ctrl", "beta1.est", "beta1.sd", "hsq.est", "nearby", "pos")
write.table(res_single_casecontrol_df, paste("results/chr", CHROM, "/res_single_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

res_single_caseonly <- rareMETALS2.single(score_stats_files_caseonly,range=range,alternative="two.sided",ix.gold=1,callrate.cutoff=0,hwe.cutoff=0)
res_single_caseonly$integratedData <- NULL
res_single_caseonly$raw.data <- NULL
res_single_caseonly$clean.data <- NULL
res_single_caseonly_df <- do.call(cbind, lapply(res_single_caseonly, data.frame, stringsAsFactors=FALSE))
names(res_single_caseonly_df) <- c("p.value", "ref", "alt", "statistic", "direction.by.study", "anno", "af", "afCase", "afCtrl", "QC.by.study", "no.sample", "no.case", "no.ctrl", "beta1.est", "beta1.sd", "hsq.est", "nearby", "pos")
write.table(res_single_caseonly_df, paste("results/chr", CHROM, "/res_single_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)