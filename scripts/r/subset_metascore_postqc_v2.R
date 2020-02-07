library(tidyverse)
library(stringr)
METASCORE <- Sys.getenv("METASCORE")
POSTQC_VARIANTS <- Sys.getenv("POSTQC_VARIANTS")
OUTPUT_FILENAME <- METASCORE %>% str_replace(".MetaScore", "_hwe5_persubstudy50_noheader.MetaScore")
metascore <- read.table(METASCORE, sep="\t", header=TRUE, stringsAsFactors=FALSE)
postqc_variants <- read_delim(POSTQC_VARIANTS, delim="\t", col_types=cols(.default="c"), col_names=FALSE)$X1
metascore$position <- str_c(metascore$CHROM, ":", metascore$POS)
metascore_postqc <- subset(metascore, position%in%postqc_variants)
metascore_postqc <- metascore_postqc %>% separate(HWE_PVALUE, c("p_hwe_all", "p_hwe_cases", "p_hwe_controls"), sep=":", remove = FALSE)
metascore_postqc_hwe <- subset(metascore_postqc, as.numeric(p_hwe_all) > 1e-5 & as.numeric(p_hwe_controls) > 1e-5 & as.numeric(p_hwe_cases) > 1e-5)
metascore_postqc_hwe$p_hwe_all <- NULL
metascore_postqc_hwe$p_hwe_cases <- NULL
metascore_postqc_hwe$p_hwe_controls <- NULL
metascore_postqc_hwe <- metascore_postqc_hwe %>% separate(CALL_RATE, c("callrate_all", "callrate_cases", "callrate_controls"), sep=":", remove = FALSE)
metascore_postqc_hwe_callrate <- subset(metascore_postqc_hwe, as.numeric(callrate_cases) > 0.85 & as.numeric(callrate_controls) > 0.85)
metascore_postqc_hwe_callrate$callrate_all <- NULL
metascore_postqc_hwe_callrate$callrate_cases <- NULL
metascore_postqc_hwe_callrate$callrate_controls <- NULL
metascore_postqc_hwe_callrate <- metascore_postqc_hwe_callrate %>% separate(INFORMATIVE_ALT_AC, c("ac_all", "ac_cases", "ac_controls"), sep=":", remove = FALSE)
metascore_postqc_hwe_callrate_nomono <- subset(metascore_postqc_hwe_callrate, as.numeric(ac_all) > 0)
metascore_postqc_hwe_callrate_nomono$ac_all <- NULL
metascore_postqc_hwe_callrate_nomono$ac_cases <- NULL
metascore_postqc_hwe_callrate_nomono$ac_controls <- NULL
write.table(metascore_postqc_hwe_callrate_nomono, OUTPUT_FILENAME, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
print("Mission complete!")
