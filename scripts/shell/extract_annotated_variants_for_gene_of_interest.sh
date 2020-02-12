### This is a script for extracting an annotated list of all variants in a gene(s) of interest, followed by generating some common burden results

# Let's define some variables and make a directory for our genes
export GENE="type in name of your gene/gene group here"
export GENELIST=genelist.txt
export GENELIST_DIR=genelists
export SAMPLEINFO="type the name of the sample info file here"
mkdir $GENELIST_DIR

# Now type in the genes that are part of your analysis
echo "Type your genes of interest here line by line" > $GENELIST

# Let's now extract the co-ordinates for our genes of interest
Rscript generate_gene_coordinates.R

# Generate a fam file
export SAMPLEINFO_FAMFILE=pdgsc_preqc.fam
Rscript generate_famfile_from_sampleinfo.R


while read CHROM
do
    echo "vcftools --gzvcf pdgscOnlyFeb6th2017_chr${CHROM}_variantid.vcf.gz --bed $GENELIST_DIR/genelist_chr${CHROM}.bed --out $GENELIST_DIR/${GENE}_chr${CHROM}_indqc --keep samplelist_${GENE}.txt --recode --recode-INFO-all
./anno/executable/anno -i $GENELIST_DIR/${GENE}_chr${CHROM}_indqc.recode.vcf -o $GENELIST_DIR/${GENE}_chr${CHROM}_indqc_anno.vcf.gz -r anno/resources/hs37d5.fa -g anno/resources/refFlat_hg19.txt.gz -p anno/priority.txt -c anno/codon.txt --indexOutput" > $GENELIST_DIR/annotate_vcf_${GENE}_chr${CHROM}.sh
    nohup bash $GENELIST_DIR/annotate_vcf_${GENE}_chr${CHROM}.sh &> $GENELIST_DIR/nohup_annotate_vcf_${GENE}_chr${CHROM}.log &
done < $GENELIST_DIR/chromosome_list.txt
wait


while read CHROM
do
    echo "vcftools --keep samplelist_${GENE}.txt --gzvcf pdgscOnlyFeb6th2017_VEP_annotated_chr${CHROM}.vcf.gz --bed $GENELIST_DIR/genelist_chr${CHROM}.bed --out $GENELIST_DIR/genelist_chr${CHROM} --recode --recode-INFO-all
java -jar ./snpEff/SnpSift.jar caseControl -tfam pdgsc_preqc.fam $GENELIST_DIR/genelist_chr${CHROM}.recode.vcf > $GENELIST_DIR/genelist_chr${CHROM}_casecontrol.vcf
./gatk-4.0.5.1/gatk VariantsToTable -R human_g1k_v37.fasta -V $GENELIST_DIR/genelist_chr${CHROM}_casecontrol.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F CSQ -F Cases -F Controls -F CC_DOM  -F CC_REC  -F CC_ALL  -F CC_GENO -F CC_TREND -O $GENELIST_DIR/genelist_chr${CHROM}_casecontrol.tab" > $GENELIST_DIR/vcftools_filter_chr${CHROM}.sh
    nohup bash $GENELIST_DIR/vcftools_filter_chr${CHROM}.sh &> $GENELIST_DIR/nohup_vcftools_filter_chr${CHROM}.log &
done < $GENELIST_DIR/chromosome_list.txt
wait

CHROM=14
./plink --vcf $GENELIST_DIR/${GENE}_chr${CHROM}_indqc_anno.vcf.gz --double-id --keep-allele-order --make-bed --out $GENELIST_DIR/${GENE}_indqc
./plink --bfile $GENELIST_DIR/${GENE}_indqc --missing --keep-fam samplelist_${GENE}.txt --out pdgsc_${GENE}_postqc_${GENE}
./plink --bfile $GENELIST_DIR/${GENE}_indqc --missing --keep-fam controllist_${GENE}.txt --out pdgsc_${GENE}_postqc_controls_${GENE}
./plink --bfile $GENELIST_DIR/${GENE}_indqc --missing --keep-fam caselist_${GENE}.txt --out pdgsc_${GENE}_postqc_cases_${GENE}
./plink --bfile $GENELIST_DIR/${GENE}_indqc --freq counts --keep-fam samplelist_${GENE}.txt --out pdgsc_${GENE}_postqc_${GENE}
./plink --bfile $GENELIST_DIR/${GENE}_indqc --keep-fam samplelist_${GENE}.txt --pheno phenofile_${GENE}.txt --test-missing --out pdgsc_${GENE}_${GENE} --allow-no-sex
./plink --bfile $GENELIST_DIR/${GENE}_indqc --freq --keep-fam controllist_${GENE}.txt --out pdgsc_${GENE}_postqc_controls_${GENE}
./plink --bfile $GENELIST_DIR/${GENE}_indqc --freq --keep-fam caselist_${GENE}.txt --out pdgsc_${GENE}_postqc_cases_${GENE}

mkdir miss_freq_depth_statsplink_${GENE}

ls pdgsc_${GENE}* > list_format_files.txt

### Clean up files to make easily readable in R
while read FILE
do
    cat ${FILE} | tr ' ' '\t' > ${FILE}.tsv
    cat ${FILE}.tsv | tr -s ' ' '\t' > miss_freq_depth_statsplink_${GENE}/${FILE}.tab
done < list_format_files.txt

# Now generate a snp list
Rscript generate_snplist.R

### Calculate differential & absolute depth stats from the matrices (with called genotypes only)
export CHROM=14
export AD_FILE=PDGSC_AD_matrix_calledGenotypesOnly_chr14.tab

### scratch here

export POSTQC_VARIANTLIST="keep_variants_${GENE}.txt"

Rscript keepvariants_chr.R

vcftools --gzvcf $GENELIST_DIR/${GENE}_chr14_indqc_anno.vcf.gz --out $GENELIST_DIR/${GENE}_anno_indqc_snpqc --keep samplelist_${GENE}.txt --positions keep_variants_${GENE}_chr.txt --recode --recode-INFO-all
bgzip -c $GENELIST_DIR/${GENE}_anno_indqc_snpqc.recode.vcf > $GENELIST_DIR/${GENE}_anno_indqc_snpqc.vcf.gz
tabix -p vcf $GENELIST_DIR/${GENE}_anno_indqc_snpqc.vcf.gz

vcftools --keep samplelist_${GENE}.txt --positions keep_variants_${GENE}_chr.txt --gzvcf pdgscOnlyFeb6th2017_VEP_annotated_chr14.vcf.gz --out $GENELIST_DIR/${GENE}_genelist_indqc_snpqc --recode --recode-INFO-all
java -jar ./snpEff/SnpSift.jar caseControl -tfam ${GENE}.ped $GENELIST_DIR/${GENE}_genelist_indqc_snpqc.recode.vcf > $GENELIST_DIR/${GENE}_genelist_indqc_snpqc_casecontrol.vcf
./gatk-4.0.5.1/gatk VariantsToTable -R human_g1k_v37.fasta -V $GENELIST_DIR/${GENE}_genelist_indqc_snpqc_casecontrol.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F CSQ -F Cases -F Controls -F CC_DOM  -F CC_REC  -F CC_ALL  -F CC_GENO -F CC_TREND -O $GENELIST_DIR/${GENE}_genelist_indqc_snpqc_casecontrol.tab

wget http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat_hg19.txt.gz

export PEDFILE=caseControlAnalysisDataFebruary0218_rollback_dpalt_mean_depth_sing_rate_newcovars_diffdepth25_ttest_dodgyvariants_studydummy.ped
export INVCF=$GENELIST_DIR/${GENE}_anno_indqc_snpqc.vcf.gz
export ANNO_LOF='Frameshift|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'
export ANNO_PROCO='Nonsynonymous|CodonGain|CodonLoss|Frameshift|Normal_Splice_Site|Essential_Splice_Site|Start_Loss|Start_Gain|Stop_Loss|Stop_Gain'


nohup ./rvtests_june/rvtests/executable/rvtest --geneFile refFlat_hg19.txt.gz --numThread 8 --freqUpper 0.01 --annoType ${ANNO_LOF} --out ${GENELIST_DIR}/${GENE}_lof_maf01_june --burden cmc,zeggini,mb,fp --vt price,analytic --kernel skat,skato --inVcf ${INVCF} --pheno ${PEDFILE} --pheno-name case2_cont1 --covar ${PEDFILE} --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sexCov,AGE,mean_depth --impute mean &> ${GENELIST_DIR}/nohup_${GENE}_lof_maf01.log &
nohup ./rvtests_june/rvtests/executable/rvtest --geneFile refFlat_hg19.txt.gz --numThread 8 --freqUpper 0.05 --annoType ${ANNO_LOF} --out ${GENELIST_DIR}/${GENE}_lof_maf05_june --burden cmc,zeggini,mb,fp --vt price,analytic --kernel skat,skato --inVcf ${INVCF} --pheno ${PEDFILE} --pheno-name case2_cont1 --covar ${PEDFILE} --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sexCov,AGE,mean_depth --impute mean &> ${GENELIST_DIR}/nohup_${GENE}_lof_maf05.log &

nohup ./rvtests_june/rvtests/executable/rvtest --geneFile refFlat_hg19.txt.gz --numThread 8 --freqUpper 0.01 --annoType ${ANNO_PROCO} --out ${GENELIST_DIR}/${GENE}_proco_maf01_june --burden cmc,zeggini,mb,fp --vt price,analytic --kernel skat,skato --inVcf ${INVCF} --pheno ${PEDFILE} --pheno-name case2_cont1 --covar ${PEDFILE} --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sexCov,AGE,mean_depth --impute mean &> ${GENELIST_DIR}/nohup_${GENE}_proco_maf01.log &
nohup ./rvtests_june/rvtests/executable/rvtest --geneFile refFlat_hg19.txt.gz --numThread 8 --freqUpper 0.05 --annoType ${ANNO_PROCO} --out ${GENELIST_DIR}/${GENE}_proco_maf05_june --burden cmc,zeggini,mb,fp --vt price,analytic --kernel skat,skato --inVcf ${INVCF} --pheno ${PEDFILE} --pheno-name case2_cont1 --covar ${PEDFILE} --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sexCov,AGE,mean_depth --impute mean &> ${GENELIST_DIR}/nohup_${GENE}_proco_maf05.log &