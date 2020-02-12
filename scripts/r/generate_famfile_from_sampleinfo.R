library(stringr)
library(tidyverse)
sampleinfo <- Sys.getenv("SAMPLEINFO")
sampleinfo_famfile <- Sys.getenv("SAMPLEINFO_FAMFILE")
samples <- read_delim(sampleinfo, delim="\t", col_types=cols(.default="c"))
fam <- as.data.frame(samples$ID)
names(fam)[1] <- "fid"
fam$iid <- samples$ID
fam$mid <- 0
fam$pid <- 0
fam$sex <- samples$male1_Female2
fam$pheno <- as.numeric(samples$CASE)+1
write.table(fam, sampleinfo_famfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
