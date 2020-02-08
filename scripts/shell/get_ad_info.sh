##### GET AD INFO

java -jar picard.jar CreateSequenceDictionary REFERENCE=human_g1k_v37.fasta OUTPUT=human_g1k_v37.dict

### Set the prefix for the VCF. Need to change this and replace all file names as variables, as currently all scripts reference the actual file name and need to be adapted.
export VCF_PREFIX=pdgscOnlyFeb6th2017

### Use GATK VariantsToTable to extract AD fields (note, here I'm using pdgscOnlyFeb6th2017_chr${CHROM}_variantid.vcf.gz, which has variant ID's set as chr:pos. If the variant ID's haven't been set like this, probably better to extract the CHROM and POS columns with VariantsToTable for each variant instead of ID, then merge these into chr:pos
echo '#!/bin/bash/
for CHROM in {1..22} X
do
    tabix -p vcf ${VCF_PREFIX}_chr${CHROM}_variantid.vcf.gz
    java -jar GenomeAnalysisTK.jar -R human_g1k_v37.fasta -T VariantsToTable -V ${VCF_PREFIX}_chr${CHROM}_variantid.vcf.gz -F ID -F ALT -GF AD -o ${VCF_PREFIX}_chr${CHROM}_variantid_AD_fields.tab
done' > GATK_get_AD.sh
nohup bash GATK_get_AD.sh &> nohup_GATK_get_AD.log &
wait

### Now create subsetted files of each chromosome, where each file has 10,000 rows. This is to make things run more efficiently on smaller files, and might not be necessary for Merck data. Also extract header.
mkdir GATK_AD_files_split
for CHROM in {1..22} X
do
    split -l 10000 ${VCF_PREFIX}_chr${CHROM}_variantid_AD_fields.tab GATK_AD_files_split/split_AD_fields_chr${CHROM}
done
head -n 1 ${VCF_PREFIX}_chr1_variantid_AD_fields.tab > AD_FIELDS_HEADER.txt

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
        nohup Rscript process_ad_depths.R &> nohup_process_ad_depths_${FILE}.log &
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
        nohup Rscript generate_final_AD_matrices.R &> nohup_generate_final_AD_matrices_${FILE}.log &
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
    echo "rm ${VCF_PREFIX}_chr${CHROM}_variantid.vcf
    bgzip -c -d ${VCF_PREFIX}_chr${CHROM}_variantid.vcf.gz > ${VCF_PREFIX}_chr${CHROM}_variantid.vcf
    grep -v ^## ${VCF_PREFIX}_chr${CHROM}_variantid.vcf | cut -f1,2,10- | sed 's/:\S*//g' > allsamples_genotypes/genotypes_chr${CHROM}.tab
    rm ${VCF_PREFIX}_chr${CHROM}_variantid.vcf" > get_genotypes_chr${CHROM}.sh
	nohup bash get_genotypes_chr${CHROM}.sh &> nohup_get_genotypes_chr${CHROM}.log &
done
wait

### Now generate the final AD and DP_ALT files
mkdir final_AD_DPALT_files
for CHROM in {1..22} X
do
    export GT_FILE=genotypes_chr${CHROM}.tab
    export AD_FILE=AD_matrix_chr${CHROM}.csv
    nohup Rscript finalise_AD_DPALT.R &> nohup_finalise_AD_DPALT_chr${CHROM}.log &
done
wait

### Now generate subsetted matrices (each one containing 750 samples) to make working with these matrices to calculate individual level statistics easier. This might not be necessary with the Merck data, as the number of samples is smaller.
mkdir final_AD_DPALT_files_per_chromosome_sample_subsets
for CHROM in {1..22} X
do
    export AD_FILE=PDGSC_AD_matrix_chr${CHROM}.tab
    export DPALT_FILE=PDGSC_DPALT_matrix_chr${CHROM}.tab
    nohup Rscript split_ad_files_samplesubsets.R &> nohup_split_ad_files_samplesubsets_chr${CHROM}.log &
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
    nohup Rscript generate_ad_matrices_called_geno_only.R &> nohup_generate_ad_matrices_called_geno_only_chr${CHROM}.log &
done
wait

### Now make sample subsets (of 750 of each) of these depth matrices with called genotypes only
mkdir AD_calledGenoOnly_per_chromosome_sample_subsets
for CHROM in {1..22} X
do
    export AD_FILE=PDGSC_AD_matrix_calledGenotypesOnly_chr${CHROM}.tab
    nohup Rscript split_ad_calledGenoOnly_samplesubsets.R &> nohup_split_ad_calledGenoOnly_samplesubsets_chr${CHROM}.log &
done
wait

### Now merge all the chromosomes for each sample subset. First create a list of file identifiers so that can loop through them.
ls AD_calledGenoOnly_per_chromosome_sample_subsets/PDGSC_AD_matrix_calledGenotypesOnly_chr10_samplesubset_* > AD_matrix_samplesubset_list.txt
sed -i 's/AD_calledGenoOnly_per_chromosome_sample_subsets\/PDGSC_AD_matrix_calledGenotypesOnly_chr10_//g' AD_matrix_samplesubset_list.txt

mkdir final_AD_calledGenoOnly_files_sample_subsets
while read SUBSET
do
    export AD_SUBSET=$SUBSET
	nohup Rscript merge_per_chromosome_AD_calledGenoOnly_sample_subsets.R &> nohup_merge_per_chromosome_AD_calledGenoOnly_sample_subsets_chr${SUBSET}.log &
done < AD_matrix_samplesubset_list.txt
wait


############# Generate depth stats to get distributions and remove outliers ###################


### To calculate mean depth per sample (done with all variants here, but if want to only include post-QC variants, adapt and uncomment the lines in the script that subset the matrix to only include rows of QC'd variants
for CHROM in {1..22} X
do
    export AD_FILE=PDGSC_AD_matrix_calledGenotypesOnly_chr${CHROM}.tab
    nohup Rscript calculate_mean_depths_calledGenoOnly.R &> nohup_calculate_mean_depths_calledGenoOnly_chr${CHROM}.log &
done
wait

### Then can plot and remove outliers

########### To calculate per sample DP_ALT stats (done with all variants here, but could also be done with post-QC variants). If want to do with QC'd variants, then edit the R script to subset the matrices first to only include variants that have passed QC before calculating means:
mkdir DPALT_per_sample_means
while read SUBSET
do
    export DPALT_SUBSET=PDGSC_DPALT_matrix_chrALL_${SUBSET}
    nohup Rscript calculate_mean_dpalt.R &> nohup_calculate_mean_depths_subset_${SUBSET}.log &
done < AD_matrix_samplesubset_list.txt
wait

### Then can plot and remove outliers
### Can also rerun pseq istats here (with all variants or post-QC variants) to look for outliers.