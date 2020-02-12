library(tidyverse)
freq <- read_delim("miss_freq_depth_statsplink_${GENE}/pdgsc_${GENE}_postqc_${GENE}.frq.counts.tab", delim="\t", col_types=cols(.default="c"))
write.table(freq$SNP, "${GENE}_variants.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
