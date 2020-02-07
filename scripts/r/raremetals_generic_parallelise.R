# Import libraries
library(rareMETALS)
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
ANNOTYPE <- Sys.getenv("ANNOTYPE")
STUDYTYPE <- Sys.getenv("STUDYTYPE")
TEST <- Sys.getenv("TEST")
MAF_CUTOFF <- as.numeric(Sys.getenv("MAF_CUTOFF"))
resultsdir <- Sys.getenv("resultsdir")

# Make the score and cov files
cov_files_casecontrol <- c(london_casecontrol_covar, merck_casecontrol_covar)
score_stats_files_casecontrol <- c(anno_london_casecontrol_assoc, anno_merck_casecontrol_assoc)
cov_files_caseonly <- c(london_caseonly_covar, merck_caseonly_covar)
score_stats_files_caseonly <- c(anno_london_caseonly_assoc, anno_merck_caseonly_assoc)

# Make annotations for gene-based tests
AllVars <- 'gene'
Lof <- 'Frameshift|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
LofPlusIndel <- 'Insertion|Deletion|Frameshift|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
ProCo <- 'Nonsynonymous|CodonGain|CodonLoss|Frameshift|Normal_Splice_Site|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
ProCoPlusIndel <- 'Insertion|Deletion|Nonsynonymous|CodonGain|CodonLoss|Frameshift|Normal_Splice_Site|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'

genelist_df <- read.table(genes, header=FALSE, stringsAsFactors=FALSE, sep="")
gene <- genelist_df$V1
rm("genelist_df")
gc()

genepos <- find.gene.chrpos(gene)$tabix
drop_index <- NULL
for (i in 1:length(genepos)){
    if (!isTabixRange(genepos[i])){
        drop_index <- c(drop_index, i)
    }
}

if (!is.null(drop_index)){
    genepos <- genepos[-drop_index]
    gc()
    gene <- gene[-drop_index]
    gc()
    write.table(gene[drop_index], paste(resultsdir, "/dropped/", CHROM, "_", ANNOTYPE, "_", MAF_CUTOFF, "_", TEST, "_", STUDYTYPE, ".tab",  sep=""), sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

#counter <- NULL
#generesults <- list()
#chr17 (KCs just for VT)
badgenes <- which(gene=="MIR4315-1" | gene=="MIR4315-2" | gene=="KCNJ12" | gene=="KCNJ18" | gene=="TRIM37")
#chr19 (SNARs for both SKAT & VT, rest for just VT)
#badgenes <- which(gene=="SNAR-A10" | gene=="SNAR-A11" | gene=="SNAR-A14" | gene=="SNAR-A3" | gene=="SNAR-A4" | gene=="SNAR-A5" | gene=="SNAR-A6" | gene=="SNAR-A7" | gene=="SNAR-A8" | gene=="SNAR-A9" | gene=="MUC16" | gene=="CIC" | gene=="SHANK1")
#VT - chr1-16
#badgenes <- which(gene=="ANKRD26" | gene=="CD163L1" | gene=="COG5" | gene=="CORIN" | gene=="EPS8" | gene=="ERCC8" | gene=="LOC442028" | gene=="MUC6" | gene=="PIK3R4" | gene=="PLA2G4E" | gene=="PPRC1" | gene=="PRH1-PRR4" | gene=="RIPK2" | gene=="ZBTB11" | gene=="GXYLT1" | gene=="MEMO1" | gene=="SMYD1" | gene=="TAS2R10" | gene=="TAS2R13" | gene=="TAS2R14" | gene=="TAS2R19" | gene=="TAS2R20" | gene=="TAS2R30" | gene=="TAS2R31" | gene=="TAS2R42" | gene=="TAS2R43" | gene=="TAS2R46" | gene=="TAS2R50" | gene=="TAS2R7" | gene=="TAS2R8" | gene=="TAS2R9" | gene=="ZFC3H1" | gene=="SP100" | gene=="TEKT4")
#chr20 (VT)
#badgenes <- which(gene=="GPCPD1")
#chr22 (VT)
#badgenes <- which(gene=="NEFH")

if (!is.null(badgenes)){
    genepos <- genepos[-badgenes]
    gc()
    gene <- gene[-badgenes]
    gc()
}

stopifnot(length(gene)==length(genepos))

########## GENE-BASED ##############

# AllVars
if(STUDYTYPE=="casecontrol"){
    generesult <- rareMETALS2.range(ANNO=eval(parse(text = ANNOTYPE)),eval(parse(text = paste("score_stats_files_", STUDYTYPE ,sep=""))),eval(parse(text = paste("cov_files_", STUDYTYPE ,sep=""))),range=genepos,range.name=gene,test = TEST,maf.cutoff = MAF_CUTOFF,alternative ="two.sided",ix.gold = 1,callrate.cutoff = 0,hwe.cutoff = 0,hwe.ctrl.cutoff=0,max.VT = NULL)
    hsq <- NULL
    for (i in 1:nrow(generesult$res.out)){
	    hsq[i] <- generesult$res.list[[i]]$hsq
    }
    generesult$res.final <- cbind(generesult$res.out, hsq)
    write.table(generesult$res.final, paste(resultsdir, "/chr", CHROM, "/", ANNOTYPE, "_", MAF_CUTOFF, "_", TEST, "_", STUDYTYPE, "_", "chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

if(STUDYTYPE=="caseonly"){
    generesult <- rareMETALS.range(ANNO=eval(parse(text = ANNOTYPE)),eval(parse(text = paste("score_stats_files_", STUDYTYPE ,sep=""))),eval(parse(text = paste("cov_files_", STUDYTYPE ,sep=""))),range=genepos,range.name=gene,test = TEST,maf.cutoff = MAF_CUTOFF,alternative ="two.sided",ix.gold = 1,callrate.cutoff = 0,hwe.cutoff = 0,max.VT = NULL)
    hsq <- NULL
    for (i in 1:nrow(generesult$res.out)){
	    hsq[i] <- generesult$res.list[[i]]$hsq
    }
    generesult$res.final <- cbind(generesult$res.out, hsq)
    write.table(generesult$res.final, paste(resultsdir, "/chr", CHROM, "/", ANNOTYPE, "_", MAF_CUTOFF, "_", TEST, "_", STUDYTYPE, "_", "chr", CHROM, ".tab", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}
# Collect the heterogeneity estimates from the list element into the output table

print(paste("Mission complete: Chromosome", CHROM, ANNOTYPE, TEST, STUDYTYPE, "MAF", MAF_CUTOFF, sep=" "))
