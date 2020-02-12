library(tidyverse)
POSTQC_VARIANTLIST <- Sys.getenv("POSTQC_VARIANTLIST")
postqc_variantlist <- read_delim(POSTQC_VARIANTLIST, delim="\t", col_types = cols(.default="c"), col_names=FALSE)
names(postqc_variantlist)[1] <- "variant"
postqc_variantlist_chrpos <- postqc_variantlist %>% separate(variant, c("CHROM", "POS"), sep=":", remove=FALSE)
chromosomes <- unique(postqc_variantlist_chrpos$CHROM)
lengthsum <- NULL
for (i in seq_along(chromosomes)){
    chromosome <- chromosomes[i]
    keepvariants_group1_chr <- subset(postqc_variantlist_chrpos, CHROM==chromosome)
    keepvariants_group1_chr$variant <- NULL
    lengthsum[i] <- dim(keepvariants_group1_chr)[1]
    write.table(keepvariants_group1_chr, str_c("keep_variants_chr.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
