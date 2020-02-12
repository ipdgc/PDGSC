library(stringr)
library(tidyverse)
library(rareMETALS)
GENELIST <- Sys.getenv("GENELIST")
GENELIST_DIR <- Sys.getenv("GENELIST_DIR")
genelist <- read.table(GENELIST, sep="", header=FALSE, stringsAsFactors=FALSE)
names(genelist)[1] <- "gene"
genelist$range <- find.gene.chrpos(genelist$gene)$tabix
genelist <- genelist %>% separate(range, into = c("chrom", "range"), sep = ":") %>%
    separate(range, into = c("chromStart", "chromEnd"), sep = "-")
chromosomes <- unique(genelist$chrom)
for (i in seq_along(chromosomes)){
    chromosome <- chromosomes[i]
    bed_chr <- subset(genelist, chrom==chromosome)
    bed_chr$gene <- NULL
    write.table(bed_chr, str_c(GENELIST_DIR, "/genelist_chr", chromosome, ".bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}
write.table(chromosomes, str_c(GENELIST_DIR, "/chromosome_list.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
