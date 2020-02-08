# Install required tools
sudo bash install_tools.sh

export SUBSTUDY=group2_matched_2to1_depth_anc

# group2_matched_all_depth_anc_chr1.singleVar_noMAC_03Sep_CaseControl_substudy_n_sites_substudies_UCLonly.MetaScore.assoc.gz
# group2_matched_2to1_depth_anc_chr1.singleVar_noMAC_04Sep_CaseControl_substudy_n_sites_substudies_UCLonly.MetaScore.assoc

while read CASECONTROLSUBSET
do
for CHROM in {1..22}
do
    export METASCORE_STEM=${SUBSTUDY}_chr${CHROM}.singleVar_noMAC_04Sep_CaseControl_n_sites_${CASECONTROLSUBSET}
    echo "bgzip -d -c ${METASCORE_STEM}.MetaScore.assoc.gz > ${METASCORE_STEM}.MetaScore.assoc" > unzip_metascore_no${CASECONTROLSUBSET}_chr${CHROM}.sh
    nohup bash unzip_metascore_no${CASECONTROLSUBSET}_chr${CHROM}.sh &> nohup_unzip_metascore_no${CASECONTROLSUBSET}_chr${CHROM}.log &
done
done < list_of_casecontrol_subsets.txt

while read CASECONTROLSUBSET
do
for CHROM in {1..22}
do
    export METASCORE_STEM=${SUBSTUDY}_chr${CHROM}.singleVar_noMAC_04Sep_CaseControl_n_sites_${CASECONTROLSUBSET}
	export POSTQC_VARIANTS=group2_postqc_keep_variants_2to1_depth_anc_hwe5_persubstudy50_${CASECONTROLSUBSET}.txt
    export METASCORE=${METASCORE_STEM}.MetaScore.assoc
	nohup Rscript subset_metascore_postqc.R &> nohup_subset_metascore_postqc_chr${CHROM}.log &
done
done < list_of_casecontrol_subsets.txt

export QC_MASTERLIST=qc_masterlist_${SUBSTUDY}.txt

while read CASECONTROLSUBSET
do
for CHROM in {1..22}
do
    export METASCORE_STEM=${SUBSTUDY}_chr${CHROM}.singleVar_noMAC_04Sep_CaseControl_n_sites_${CASECONTROLSUBSET}
    export METASCORE=${METASCORE_STEM}
	grep "^\#\#" ${METASCORE_STEM}.MetaScore.assoc > metascore_header_grouped${CASECONTROLSUBSET}_chr${CHROM}.txt
	cat metascore_header_grouped${CASECONTROLSUBSET}_chr${CHROM}.txt  ${METASCORE_STEM}_hwe5_persubstudy50_noheader.MetaScore.assoc >  ${METASCORE_STEM}_hwe5_persubstudy50.MetaScore.assoc
	bgzip -c ${METASCORE_STEM}_hwe5_persubstudy50.MetaScore.assoc > ${METASCORE_STEM}_hwe5_persubstudy50.MetaScore.assoc.gz
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${METASCORE_STEM}_hwe5_persubstudy50.MetaScore.assoc.gz
done
done < list_of_casecontrol_subsets.txt
