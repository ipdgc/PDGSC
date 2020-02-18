pca_vectors <- Sys.getenv("PCA_VECTORS")
fam_file <- Sys.getenv("FAM_FILE")

library("ggplot2")
library("plyr")
library("RColorBrewer")

pca <- read.table(pca_vectors, header=TRUE, sep="\t", stringsAsFactors=FALSE)
fam <- read.table(fam_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
fam$V6[fam$V6==1] <- "Control"
fam$V6[fam$V6==2] <- "Case"
fam$V6[fam$V6==0] <- "Missing"
fam$V6 <- as.factor(fam$V6)
hapfam <- read.table("hapmap3_r3_b36_fwd.consensus.qc.poly.fam", sep=" ", stringsAsFactors=FALSE, header=FALSE)
hapmappop <- read.table("relationships_w_pops_041510.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
names(hapfam)[1] <- "FID"
hapfampop <- join(hapfam, hapmappop, by="FID")
hapfampop_trimmed <- hapfampop[,c(1:5,12)]
names(hapfampop_trimmed)[c(1,6)] <- c("V1","V6")
hapfampop_trimmed$V6 <- as.factor(hapfampop_trimmed$V6)
totalfam <- rbind(fam, hapfampop_trimmed)
names(totalfam)[1] <- "IID"
total_pca_data <- join(pca, totalfam, by="IID")

p <- qplot(PC1, PC2, data=total_pca_data, colour=V6)
p <- p + labs(title = "Principal component analysis with HapMap phase 3 data")
p <- p + labs(x= "Principal component 1")
p <- p + labs(y= "Principal component 2")
p <- p + scale_colour_brewer(palette="Dark2")

ggsave("pca_brew.pdf")
