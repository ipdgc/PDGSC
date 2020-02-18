### START QC PIPELINE ###
export ORIGINAL_VCF= # name of the original VCF
export VCF_STEM=${ORIGINAL_VCF//.recode.vcf/}
export SAMPLE_INFO= # sample sheet
export AD_MEANS="PDGSC_AD_means_nsites_calledGenoOnly.tab"
export DPALT_MEANS="PDGSC_DPALT_full_means.tab"
export FREEMIX="adsp.WES_target_intersect.txt"

# Run Pseq to obtain individual level summary metrics
nohup ./plinkseq-0.10/pseq ${ORIGINAL_VCF} i-stats &> pdgsc_allsamples_preqc.istats &

# Convert to Plink and filter for call rate & MAF
./plink --vcf ${ORIGINAL_VCF} --double-id --maf 0.05 --geno 0.05 --make-bed --out ${VCF_STEM}_maf05_geno05

# Update the phenotypes in the .fam file
export FAM_FILE=${VCF_STEM}_maf05_geno05.fam
export UPDATED_FAM_FILE=${VCF_STEM}_newfam.fam
export MISSING_CLINICAL=${VCF_STEM}_missingpheno.txt
Rscript ./pdgsc/scripts/r/update_famfile.R
cp ${UPDATED_FAM_FILE} ${FAM_FILE}

# Prune
./plink --bfile ${VCF_STEM}_maf05_geno05 --set-missing-var-ids @:# --make-bed --out ${VCF_STEM}_maf05_geno05_missid
./plink --bfile ${VCF_STEM}_maf05_geno05_missid --indep-pairwise 50 5 .5 --out prune4QC
./plink --bfile ${VCF_STEM}_maf05_geno05_missid --extract prune4QC.prune.in --make-bed --out ${VCF_STEM}_maf05_geno05_pruned

# Sex check
./plink --bfile ${VCF_STEM}_maf05_geno05_pruned --out ${VCF_STEM}_maf05_geno05_pruned_splitsex --make-bed --split-x b37
./plink --bfile ${VCF_STEM}_maf05_geno05_pruned_splitsex --check-sex 0.5 0.5 --out pdgsc_sexcheck

# Heterozygosity
./plink --bfile ${VCF_STEM}_maf05_geno05_pruned --het --out pdgsc_hetcheck

# Population structure: download and unzip Hapmap data
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/relationships_w_pops_041510.txt"
gunzip hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz
gunzip hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz
./plink --file hapmap3_r3_b36_fwd.consensus.qc.poly --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf05 --make-bed --maf 0.05

# Extract list of all SNP names
awk '{print $2}' hapmap3_r3_b36_fwd.consensus.qc.poly_maf05.bim > hapmap_snp_list
awk '{print $2}' ${VCF_STEM}_maf05_geno05.bim > ${VCF_STEM}_maf05_geno05_snplist

# Split the X-chromosome
./plink --bfile ${VCF_STEM}_maf05_geno05 --out ${VCF_STEM}_maf05_geno05_splitsex --make-bed --split-x b37
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap3_r3_b36_fwd.consensus.qc.poly_splitsex --make-bed --split-x b36

# Only keep consensus list of variants between PDGSC and Hapmap
./plink --bfile ${VCF_STEM}_maf05_geno05 --extract hapmap_snp_list --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly_maf05 --extract ${VCF_STEM}_maf05_geno05_snplist --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps

# Merge Hapmap and PDGSC
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps --bmerge hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged
# Flip mismatches
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps --flip ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged-merge.missnp --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped
# Merge again
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps --bmerge hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged
# Now exclude mismatches
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped --exclude ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged-merge.missnp --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped_nomismatch
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps --exclude ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged-merge.missnp --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_nomismatch
# Update chromosomes and positions
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped_nomismatch --update-chr ${VCF_STEM}_maf05_geno05_hapmapsnps_nomismatch.bim 1 2 --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped_nomismatch_updatedchr
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped_nomismatch_updatedchr --update-map ${VCF_STEM}_maf05_geno05_hapmapsnps_nomismatch.bim 4 2 --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped_nomismatch_updatedchrpos
# Now do final merge
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps_nomismatch --bmerge hapmap3_r3_b36_fwd.consensus.qc.poly_maf05_pdgscsnps_flipped_nomismatch_updatedchrpos --allow-no-sex --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged

# Make a list of palindromic variants
export BIM_FILE=${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged.bim
Rscript ./pdgsc/scripts/r/palindrome.R

# Exclude palindromic variants
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged --exclude palindromic_snps.txt --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome

# Only keep variants with total MAF > 0.05 across Hapmap and PDGSC
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome --maf 0.05 --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome_maftotal05

# Prune for PCA
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome_maftotal05 --indep-pairwise 50 5 0.5 --out prune_prepca

# Exclude samples that failed heterozygosity, sex check, or istats metrics, and only keep pruned variants (for PCA)
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome_maftotal05 --extract prune_prepca.in --remove-fam checkhet_sexcheck_istats_exclusion.txt --make-bed --out ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome_maftotal05_pruned_prepca

# PCA plot
./plink --bfile ${VCF_STEM}_maf05_geno05_hapmapsnps_hapmapmerged_nopalindrome_maftotal05_pruned_prepca --pca header tabs --out pdgsc_hapmap_pca
export PCA_VECTORS=pdgsc_hapmap_pca.eigenvec
Rscript ./pdgsc/scripts/r/plot_pca.R


# Now collate all results to exclude samples that fail freemix, depth statistics (obtained using get_ad_info.sh) istats metrics, sex check, heterozygosity, are population outliers or are age under 18. The cutoffs for these were determined by looking at histograms.
export ISTATS=pdgsc_allsamples_preqc.istats
export SEXCHECK=pdgsc_sexcheck.sexcheck
export HETEROZYGOSITY= pdgsc_hetcheck.het
export PCA_OUTLIERS=pca_outliers.tab
export IND_QC_OUTPUT=mean_depth_n_sites_dpalt_freemix_istats_het_sex_pca_underage_exclusions.list
# First convert the Plink spacing format into something readable in R
for FILE in $SEXCHECK $HETEROZYGOSITY
do
    cat $FILE | tr -s ' ' '\t' > $FILE.tsv
    cat $FILE.tsv | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > $FILE.tab
done
Rscript ./pdgsc/scripts/r/collect_individual_qc.R

# Now use Plink to exclude the samples that have failed QC and finish by filtering for relatedness
./plink --bfile ${VCF_STEM}_maf05_geno05_pruned --remove-fam mean_depth_n_sites_dpalt_freemix_istats_het_sex_pca_underage_exclusions.list --rel-cutoff 0.125 --out ${VCF_STEM}_indqc --make-bed

# And now we have our final list of samples that survived QC:
awk '{print $1}' ${VCF_STEM}_indqc > pdgsc_postqc_samples.txt