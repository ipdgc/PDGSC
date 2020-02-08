export SUBSTUDY=group2_matched_2to1_depth_anc
export resultsdir="results_${SUBSTUDY}"

export GROUP=group2
export SAMPLELIST_DIR=clustering_results_keepmono_50depth_50anc_2to1
export CASECONTROLPEDFILE=caseControlAnalysisData230818_${SUBSTUDY}.ped

export SAMPLELIST=clustering_results_keepmono_50depth_50anc_2to1/samplelist_matched_group2.txt
export CASELIST=clustering_results_keepmono_50depth_50anc_2to1/caselist_matched_group2.txt
export CONTROLLIST=clustering_results_keepmono_50depth_50anc_2to1/controllist_matched_group2.txt
export PHENOFILE=clustering_results_keepmono_50depth_50anc_2to1/phenofile_matched_group2.txt
export OUTPUT_STEM=pdgsc_group2_matched_2to1_depth_anc

cd

Rscript generate_analysis_subsets.R

while read CASECONTROLSUBSET
do
CASECONTROLPEDFILE=caseControlAnalysisData040918_group2_matched_2to1_depth_anc_${CASECONTROLSUBSET}.ped
    for CHROM in {1..22}
    do
        export INPUTVCF=pdgscOnlyFeb6th2017_chr${CHROM}.recode.vcf.gz
        export RVTESTS_OUTPUT_FILENAME_NOMAC=${SUBSTUDY}_chr${CHROM}.singleVar_noMAC_05Sep_CaseControl_${CASECONTROLSUBSET}
#        nohup ./rvtests_june/rvtests/executable/rvtest --rangeFile chr${CHROM}interval.txt --out ${RVTESTS_OUTPUT_FILENAME_NOMAC} --meta score,cov --inVcf ${INPUTVCF} --pheno ${CASECONTROLPEDFILE} --pheno-name case2_cont1 --covar ${CASECONTROLPEDFILE} --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sexCov,STUDY_matched_Baylor,STUDY_matched_cambridge,STUDY_matched_Luebeck,STUDY_matched_majamaa,STUDY_matched_oslo,STUDY_matched_oxford,STUDY_matched_tuebingen,STUDY_matched_UCL1,STUDY_matched_UCL2,mean_depth --impute mean &> nohup_${RVTESTS_OUTPUT_FILENAME_NOMAC}.log &
        nohup ./rvtests_june/rvtests/executable/rvtest --rangeFile chr${CHROM}interval.txt --out ${RVTESTS_OUTPUT_FILENAME_NOMAC} --meta score,cov --inVcf ${INPUTVCF} --pheno ${CASECONTROLPEDFILE} --pheno-name case2_cont1 --covar ${CASECONTROLPEDFILE} --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,sexCov,mean_depth --impute mean &> nohup_${RVTESTS_OUTPUT_FILENAME_NOMAC}.log &
    done
done < list_of_casecontrol_subsets_leftovers.txt
