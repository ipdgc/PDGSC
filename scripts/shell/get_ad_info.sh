##### This script generates depth metrics for each sample at each variant by adding up the AD fields (addresses a bug in the VCF whereby all non-reference DP calls are set as missing)

export BUCKET= # bucket location
export ORIGINAL_VCF= # name of the original VCF
export DIR_VCF_BY_CHROM= # directory within bucket for the original VCF split by chromosome and gzipped
export DIR_VCF_CHRPOS_VARIANT_ID= # directory within bucket for the original VCF split by chromosome and gzipped, with variant ID column replaced to be in chr:pos format (instead of rsid)
export VCF_STEM=${ORIGINAL_VCF//.recode.vcf/}
export SAMPLE_INFO= # name of the sample sheet

### Use GATK VariantsToTable to extract AD fields
echo '#!/bin/bash/
# for CHROM in {1..22} X
# do
CHROM=21
    ./gatk-4.1.4.1/gatk -R human_g1k_v37.fasta -T VariantsToTable -V ${VCF_STEM}_chr${CHROM}_variantid.vcf.gz -F ID -F ALT -GF AD -o ${VCF_STEM}_chr${CHROM}_variantid_AD_fields.tab
# done' > GATK_get_AD.sh
nohup bash GATK_get_AD.sh &> nohup_GATK_get_AD.log &
wait

### Now create subsetted files of each chromosome, where each file has 10,000 rows. This is to make things run more efficiently on smaller files. Also extract header.
mkdir GATK_AD_files_split
for CHROM in {1..22} X
do
    split -l 10000 ${VCF_STEM}_chr${CHROM}_variantid_AD_fields.tab GATK_AD_files_split/split_AD_fields_chr${CHROM}
done
head -n 1 ${VCF_STEM}_chr1_variantid_AD_fields.tab > AD_FIELDS_HEADER.txt

### Now make a list of all the subsetted files (so we can loop through them). Also make a list of the first subset of each chromosome (as these already contain a header).
mkdir ad_split_lists
cd GATK_AD_files_split
ls split_AD_fields_chr* > ../ad_split_lists/ad_split_chrALL.txt
ls split_AD_fields_chr*aa_header.tab > ../ad_split_lists/ad_split_firsts.txt
cd

### Now loop though the subsets and add header to each
while read SUBSET
do
    cat AD_FIELDS_HEADER.txt GATK_AD_files_split/$SUBSET > GATK_AD_files_split/${SUBSET}_header.tab
    rm GATK_AD_files_split/$SUBSET
done < ad_split_lists/ad_split_chrALL.txt

### Remove the first line from the first subset of each chromosome as these now have two header lines.
while read FILE
do
    sed -i '1d' GATK_AD_files_split/$FILE
done < ad_split_lists/ad_split_firsts.txt

### Now we're subsetting the list of all the subsets into separate lists each containing 16 subsets. This is to make sure we only run 16 subsets at a time, as otherwise overload memory.
cd ../ad_split_lists
split -l 16 ad_split_chrALL.txt ad_split_chrALL_subsetted
ls ad_split_chrALL_subsetted* > AD_split_list_of_lists.txt
cd

### Now replace the commas in the AD fields separating the depths of each allele into +'s, so that we can sum them. We do  these 16 at a time to avoid memory overload.
mkdir split_processed_depth_files
while read LIST
do
    while read SUBSET
    do
        export FILE=${SUBSET}
        nohup Rscript ./pdgsc/scripts/r/process_ad_depths.R &> nohup_process_ad_depths_${FILE}.log &
    done < ad_split_lists/$LIST
wait
done < ad_split_lists/AD_split_list_of_lists.txt

### Now we make a list of file identifiers for the matrices above, so that we can loop through them. We also split these into groups of 16 to run 16 at a time.
mkdir split_processed_depth_files_summed_DP_ALT
cd split_processed_depth_files
ls AD_depths_* > ../ad_split_lists/processed_split_list.txt
sed -i 's/AD_depths_sum_//g' ../ad_split_lists/processed_split_list.txt
cd ../ad_split_lists/
split -l 16 processed_split_list.txt processed_split_list_subset
ls processed_split_list_subset* > processed_split_list_of_lists.txt

### Now we use generate_final_AD_matrices.R to convert these matrices into numeric matrices with the depth fields summed to get the full actual depths.
cd
while read LIST
do
    while read SUBSET
    do
        export FILE=${SUBSET}
        nohup Rscript ./pdgsc/scripts/r/generate_final_AD_matrices.R &> nohup_generate_final_AD_matrices_${FILE}.log &
    done < ad_split_lists/${LIST}
wait
done < ad_split_lists/processed_split_list_of_lists.txt

### Now attach headers to each of these matrices, and create a list of them to loop through later.
mkdir per_chromosome_AD_DPALT_files
mkdir per_chromosome_file_construction_lists
for CHROM in {1..22} X
do
    cat split_processed_depth_files_summed_DP_ALT/HEADER.txt split_processed_depth_files_summed_DP_ALT/AD_matrix_split_chr${CHROM}a*.csv > per_chromosome_AD_DPALT_files/AD_matrix_chr${CHROM}.csv
    ls split_processed_depth_files_summed_DP_ALT/AD_matrix_split_chr${CHROM}a*.csv > per_chromosome_file_construction_lists/AD_matrix_chr${CHROM}_list.txt
done

### Now use cut/sed to extract the genotypes from the VCFs. This requires unzipping the VCFs first, so can take a bit of space and it's not ideal. However, GATK's VariantsToTable for some reason only outputs the actual alleles (A/T/C/G) when extracting the genotype fields, and there doesn't seem to be an option to extract them as 0/0, 0/1, 1/1 etc (at least in the last version I checked). Having them as 1s & 0s makes extracting called genotypes and genotypes with alternative allele a lot easier, so I prefer this. I did also try the vcfR package, but found that too clumsy and vulnerable to memory crashing. I also wrote a script to use VariantsToTable for this, whereby the alternative allele is compared to the alleles at each entry to see if the genotype has an alternative allele, and am happy to send that around if needed. I find this method more reliable, however.
mkdir allsamples_genotypes
for CHROM in {1..22} X
do
    echo "rm ${VCF_STEM}_chr${CHROM}_variantid.vcf
    bgzip -c -d ${VCF_STEM}_chr${CHROM}_variantid.vcf.gz > ${VCF_STEM}_chr${CHROM}_variantid.vcf
    grep -v ^## ${VCF_STEM}_chr${CHROM}_variantid.vcf | cut -f1,2,10- | sed 's/:\S*//g' > allsamples_genotypes/genotypes_chr${CHROM}.tab
    rm ${VCF_STEM}_chr${CHROM}_variantid.vcf" > get_genotypes_chr${CHROM}.sh
	nohup bash get_genotypes_chr${CHROM}.sh &> nohup_get_genotypes_chr${CHROM}.log &
done
wait

### Now generate the final AD and DP_ALT files
mkdir final_AD_DPALT_files
for CHROM in {1..22} X
do
    export GT_FILE=genotypes_chr${CHROM}.tab
    export AD_FILE=AD_matrix_chr${CHROM}.csv
    nohup Rscript ./pdgsc/scripts/r/finalise_AD_DPALT.R &> nohup_finalise_AD_DPALT_chr${CHROM}.log &
done
wait

### Now generate subsetted matrices (each one containing 750 samples) to make working with these matrices to calculate individual level statistics easier.
mkdir final_AD_DPALT_files_per_chromosome_sample_subsets
for CHROM in {1..22} X
do
    export AD_FILE=PDGSC_AD_matrix_chr${CHROM}.tab
    export DPALT_FILE=PDGSC_DPALT_matrix_chr${CHROM}.tab
    nohup Rscript ./pdgsc/scripts/r/split_ad_files_samplesubsets.R &> nohup_split_ad_files_samplesubsets_chr${CHROM}.log &
done
wait

### Now generate a list of file identifiers for these subsetted matrices, so we can loop through them later.
ls final_AD_DPALT_files_per_chromosome_sample_subsets/PDGSC_AD_matrix_chr10_samplesubset_* > AD_matrix_samplesubset_list.txt
sed -i 's/final_AD_DPALT_files_per_chromosome_sample_subsets\/PDGSC_AD_matrix_chr10_//g' AD_matrix_samplesubset_list.txt

### Now bind all the chromosomes for each of these sample subset matrices into one (as we need the full genome for each sample to generate sample summaries).
mkdir final_AD_DPALT_files_sample_subsets
while read SUBSET
do
    export AD_SUBSET=$SUBSET
	nohup Rscript merge_per_chromosome_AD_DPALT_sample_subsets.R &> nohup_merge_per_chromosome_AD_DPALT_sample_subsets_chr${SUBSET}.log &
done < AD_matrix_samplesubset_list.txt
wait

### Now generate depth matrices with called genotypes only (uncalled genotypes set as NA)
mkdir final_depths_called_geno_only_nsites
for CHROM in {1..22} X
do
    export GT_FILE=genotypes_chr${CHROM}.tab
    export AD_FILE=PDGSC_AD_matrix_chr${CHROM}.tab
    nohup Rscript ./pdgsc/scripts/r/generate_ad_matrices_called_geno_only.R &> nohup_generate_ad_matrices_called_geno_only_chr${CHROM}.log &
done
wait

### Now make sample subsets (of 750 of each) of these depth matrices with called genotypes only
mkdir AD_calledGenoOnly_per_chromosome_sample_subsets
for CHROM in {1..22} X
do
    export AD_FILE=PDGSC_AD_matrix_calledGenotypesOnly_chr${CHROM}.tab
    nohup Rscript ./pdgsc/scripts/r/split_ad_calledGenoOnly_samplesubsets.R &> nohup_split_ad_calledGenoOnly_samplesubsets_chr${CHROM}.log &
done
wait

### Now merge all the chromosomes for each sample subset. First create a list of file identifiers so that can loop through them.
ls AD_calledGenoOnly_per_chromosome_sample_subsets/PDGSC_AD_matrix_calledGenotypesOnly_chr10_samplesubset_* > AD_matrix_samplesubset_list.txt
sed -i 's/AD_calledGenoOnly_per_chromosome_sample_subsets\/PDGSC_AD_matrix_calledGenotypesOnly_chr10_//g' AD_matrix_samplesubset_list.txt

mkdir final_AD_calledGenoOnly_files_sample_subsets
while read SUBSET
do
    export AD_SUBSET=$SUBSET
	nohup ./pdgsc/scripts/r/Rscript merge_per_chromosome_AD_calledGenoOnly_sample_subsets.R &> nohup_merge_per_chromosome_AD_calledGenoOnly_sample_subsets_chr${SUBSET}.log &
done < AD_matrix_samplesubset_list.txt
wait


############# Generate depth stats to get distributions and remove outliers ###################


### To calculate mean depth per sample (done with all variants here, but if want to only include post-QC variants, adapt and uncomment the lines in the script that subset the matrix to only include rows of QC'd variants
while read SUBSET
do
    export AD_SUBSET=PDGSC_AD_matrix_calledGenoOnly_chrALL_$SUBSET
    nohup Rscript ./pdgsc/scripts/r/calculate_mean_depths_calledGenoOnly.R &> nohup_calculate_mean_depths_calledGenoOnly_chr${CHROM}.log &
done < AD_matrix_samplesubset_list.txt
wait

### Then can plot and remove outliers

########### To calculate per sample DP_ALT stats (done with all variants here, but could also be done with post-QC variants). If want to do with QC'd variants, then edit the R script to subset the matrices first to only include variants that have passed QC before calculating means:
mkdir DPALT_per_sample_means
while read SUBSET
do
    export DPALT_SUBSET=PDGSC_DPALT_matrix_chrALL_${SUBSET}
    nohup Rscript ./pdgsc/scripts/r/calculate_mean_dpalt.R &> nohup_calculate_mean_depths_subset_${SUBSET}.log &
done < AD_matrix_samplesubset_list.txt
wait

cat PDGSC_AD_means_nsites_calledGenoOnly* > PDGSC_AD_means_nsites_calledGenoOnly.tab
cat PDGSC_DPALT_full_means* > PDGSC_DPALT_full_means.tab

### Then can plot and remove outliers
