bim_file <- Sys.getenv("BIM_FILE")
bim <- read.table(bim_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
palindromic <- bim[(bim$V5=="A" & bim$V6=="T") | (bim$V5=="T" & bim$V6=="A") | (bim$V5=="C" & bim$V6=="G") | (bim$V5=="G" & bim$V6=="C"),]
write.table(palindromic$V2, "palindromic_snps.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="")