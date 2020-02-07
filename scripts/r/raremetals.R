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

########## MAF 0.01 ##############

generesults <- list()

# AllVars
generesults$AllVars_maf01_WSS_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene=gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_WSS_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_GRANVIL_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_SKAT_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_SKAT_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_VT_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf01_VT_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# ProCo
generesults$ProCo_maf01_WSS_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_WSS_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_GRANVIL_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_SKAT_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_SKAT_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_VT_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf01_VT_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# Lof
generesults$Lof_maf01_WSS_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_WSS_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_GRANVIL_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_SKAT_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_SKAT_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_VT_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf01_VT_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# ProCoPlusIndel
generesults$ProCoPlusIndel_maf01_WSS_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_WSS_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_GRANVIL_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_SKAT_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_SKAT_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_VT_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf01_VT_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# LofPlusIndel
generesults$LofPlusIndel_maf01_WSS_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_WSS_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_GRANVIL_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_SKAT_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_SKAT_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_VT_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf01_VT_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.01,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)



########## MAF 0.05 ##############

# AllVars
generesults$AllVars_maf05_WSS_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_WSS_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_GRANVIL_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_SKAT_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_SKAT_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_VT_casecontrol <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$AllVars_maf05_VT_caseonly <- rareMETALS.gene(ANNO=AllVars_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# ProCo
generesults$ProCo_maf05_WSS_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_WSS_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_GRANVIL_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_SKAT_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_SKAT_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_VT_casecontrol <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCo_maf05_VT_caseonly <- rareMETALS.gene(ANNO=ProCo_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# Lof
generesults$Lof_maf05_WSS_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_WSS_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_GRANVIL_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_SKAT_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_SKAT_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_VT_casecontrol <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$Lof_maf05_VT_caseonly <- rareMETALS.gene(ANNO=Lof_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# ProCoPlusIndel
generesults$ProCoPlusIndel_maf05_WSS_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_WSS_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_GRANVIL_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_SKAT_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_SKAT_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_VT_casecontrol <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$ProCoPlusIndel_maf05_VT_caseonly <- rareMETALS.gene(ANNO=ProCoPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# LofPlusIndel
generesults$LofPlusIndel_maf05_WSS_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_WSS_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='WSS',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_GRANVIL_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_GRANVIL_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='GRANVIL',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_SKAT_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_SKAT_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='SKAT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_VT_casecontrol <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_casecontrol,cov_files_casecontrol,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)
generesults$LofPlusIndel_maf05_VT_caseonly <- rareMETALS.gene(ANNO=LofPlusIndel_annotations,score_stats_files_caseonly,cov_files_caseonly,gene,test='VT',maf.cutoff=0.05,no.boot=0,alternative='two.sided',alpha=0,ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,gene.file="refFlat_hg19.txt.gz",max.VT=NULL)

# Collect the heterogeneity estimates from the list element into the output table
for (i in 1:length(generesults)){
	hsq <- NULL
	for (j in 1:nrow(generesults[[i]]$res.out)){
		hsq[j] <- generesults[[i]]$res.list[[j]]$hsq
	}
	generesults[[i]]$res.final <- cbind(generesults[[i]]$res.out, hsq)
	hsq <- NULL
}

# Now write out all the outputs
write.table(generesults$AllVars_maf01_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf01_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_WSS_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf01_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf01_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf01_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf01_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf01_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_VT_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf01_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf01_VT_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf01_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf01_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_WSS_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf01_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf01_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf01_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf01_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf01_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_VT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf01_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf01_VT_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf01_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf01_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_WSS_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf01_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf01_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf01_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf01_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf01_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_VT_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf01_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf01_VT_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf01_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_WSS_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_VT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf01_VT_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf01_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_WSS_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_VT_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf01_VT_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf01_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf05_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_WSS_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf05_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf05_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf05_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf05_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf05_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_VT_casecontrol$res.final, paste("results/chr", CHROM, "/AllVars_maf05_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$AllVars_maf05_VT_caseonly$res.final, paste("results/chr", CHROM, "/AllVars_maf05_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf05_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_WSS_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf05_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf05_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf05_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf05_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf05_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_VT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCo_maf05_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCo_maf05_VT_caseonly$res.final, paste("results/chr", CHROM, "/ProCo_maf05_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf05_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_WSS_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf05_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf05_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf05_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf05_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf05_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_VT_casecontrol$res.final, paste("results/chr", CHROM, "/Lof_maf05_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$Lof_maf05_VT_caseonly$res.final, paste("results/chr", CHROM, "/Lof_maf05_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_WSS_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_VT_casecontrol$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$ProCoPlusIndel_maf05_VT_caseonly$res.final, paste("results/chr", CHROM, "/ProCoPlusIndel_maf05_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_WSS_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_WSS_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_WSS_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_WSS_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_GRANVIL_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_GRANVIL_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_GRANVIL_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_GRANVIL_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_SKAT_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_SKAT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_SKAT_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_SKAT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_VT_casecontrol$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_VT_casecontrol_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(generesults$LofPlusIndel_maf05_VT_caseonly$res.final, paste("results/chr", CHROM, "/LofPlusIndel_maf05_VT_caseonly_chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)