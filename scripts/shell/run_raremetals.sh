export SUBSTUDY=group1_group2_matched_2to1_depth_anc_UCLonly_substudy_n_sites_leftovers_substudy_n_sites_majamaa_n_sites_oxford_n_sites_merck
export resultsdir="results_${SUBSTUDY}"

# rm -r ${resultsdir}
# mkdir ${resultsdir}
# mkdir ${resultsdir}/dropped


for i in 5
do
	# Export variable names so that R can recognise them
	export london_casecontrol_group1="group1_chr${i}.singleVar_noMAC_12June_CaseControl_forburden_no_n_sites_hwe5_persubstudy50"
	export london_casecontrol_group2="group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl__substudy_n_sites_substudies_UCLonly_hwe5_persubstudy50"
	export london_casecontrol_group3="group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl__substudy_n_sites_substudies_leftovers_hwe5_persubstudy50"
	export london_casecontrol_group4="group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl_n_sites_substudies_majamaa_hwe5_persubstudy50"
	export london_casecontrol_group5="group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl_n_sites_substudies_oxford_hwe5_persubstudy50"
	export merck_casecontrol="6Aug2018_Merck_070518_indqc_variantqc_nonsites_chr${i}"

	export london_casecontrol_covar_group1=${london_casecontrol_group1}.MetaCov.assoc.gz
	export london_casecontrol_covar_group2=group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl__substudy_n_sites_substudies_UCLonly.MetaCov.assoc.gz
	export london_casecontrol_covar_group3=group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl__substudy_n_sites_substudies_leftovers.MetaCov.assoc.gz
	export london_casecontrol_covar_group4=group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl_n_sites_substudies_majamaa.MetaCov.assoc.gz
	export london_casecontrol_covar_group5=group2_matched_2to1_depth_anc_chr${i}.singleVar_noMAC_04Sep_CaseControl_n_sites_substudies_oxford.MetaCov.assoc.gz
	export merck_casecontrol_covar=${merck_casecontrol}.MetaCov.assoc.gz

	export london_casecontrol_assoc_group1=${london_casecontrol_group1}.MetaScore.assoc.gz
	export london_casecontrol_assoc_group2=${london_casecontrol_group2}.MetaScore.assoc.gz
	export london_casecontrol_assoc_group3=${london_casecontrol_group3}.MetaScore.assoc.gz
	export london_casecontrol_assoc_group4=${london_casecontrol_group4}.MetaScore.assoc.gz
	export london_casecontrol_assoc_group5=${london_casecontrol_group5}.MetaScore.assoc.gz
	export merck_casecontrol_assoc=${merck_casecontrol}.MetaScore.assoc.gz

	Rscript seqmineranno_meta_london_group1_group2_group3_group4_group5_merck.R

	# Export variable names again and index the required files with Tabix
	export anno_london_casecontrol_assoc_group1=anno_${london_casecontrol_assoc_group1}
	export anno_london_casecontrol_assoc_group2=anno_${london_casecontrol_assoc_group2}
	export anno_london_casecontrol_assoc_group3=anno_${london_casecontrol_assoc_group3}
	export anno_london_casecontrol_assoc_group4=anno_${london_casecontrol_assoc_group4}
	export anno_london_casecontrol_assoc_group5=anno_${london_casecontrol_assoc_group5}
	export anno_merck_casecontrol_assoc=anno_${merck_casecontrol_assoc}
	
	export genelist="refFlat_hg19_genelist_chr${i}_unique.txt"
	
	sed -i '/^MIR/d' $genelist
	sed -i '/^SNAR/d' $genelist
	
	export i
	
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group1}
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group2}
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group3}
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group4}
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group5}
	/usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_merck_casecontrol_assoc}

	mkdir ${resultsdir}/chr${i}

### READ IN: ANNOTYPE,STUDYTYPE,TEST,MAF_CUTOFF

    for MAF_CUTOFF in 0.01 0.05
    do
        export MAF_CUTOFF=$MAF_CUTOFF
#        for ANNOTYPE in 
#        for ANNOTYPE in ProCo ProCoPlusIndel Lof LofPlusIndel
        for ANNOTYPE in ProCo ProCoPlusIndel
        do
            export ANNOTYPE=$ANNOTYPE
#            for STUDYTYPE in caseonly
            for STUDYTYPE in casecontrol
#            for STUDYTYPE in casecontrol caseonly
            do
                export STUDYTYPE=$STUDYTYPE
#                for TEST in SKAT_betakernel SKAT_linearkernel GRANVIL WSS VT
#                for TEST in GRANVIL WSS VT SKAT_betakernel
                for TEST in GRANVIL SKAT_betakernel WSS VT SKAT_linearkernel
#                 for TEST in  
                do
                    export TEST=$TEST
                    nohup Rscript --max-ppsize=500000 raremetals_generic_parallelise_london_meta_group1_group2_group3_group4_group5_merck.R &> nohup_raremetals_meta_${resultsdir}_${ANNOTYPE}_${MAF_CUTOFF}_chr${i}_${TEST}_${STUDYTYPE}.log &
#                    nohup Rscript --max-ppsize=500000 raremetals_generic_parallelise_london_meta_group1_group2_group3_group4_group5_merck.R &> nohup_raremetals_meta_${resultsdir}_${ANNOTYPE}_${MAF_CUTOFF}_chr${i}_${TEST}_${STUDYTYPE}.log &
#                wait
                done
#            wait
            done
#        wait
        done
    wait
    done
 wait
done

ls nohup_raremetals_meta_* | wc -l
grep "Mission" nohup_raremetals_meta_* | wc -l
grep -L "Mission" nohup_raremetals_meta_* | wc -l

export SUBSTUDY=group1_group2_matched_2to1_depth_anc_UCLonly_substudy_n_sites_leftovers_substudy_n_sites_majamaa_n_sites_oxford_n_sites_merck
export resultsdir="results_${SUBSTUDY}"
ls $resultsdir > results_directories.txt


wget http://launchpadlibrarian.net/74318164/tabix_0.2.5-1ubuntu1_amd64.deb
sudo dpkg -i tabix_0.2.5-1ubuntu1_amd64.deb
sudo apt-get -f install

export resultsdir="results_singlegroup_meta_group1_group2_nostudycovar_merck"
mkdir ${resultsdir}
mkdir ${resultsdir}/dropped

for i in {1..22} X
do
    # Export variable names so that R can recognise them
    export london_casecontrol_group1="london_chr${i}.singleVar_MAC3_25Feb_CaseControl_nomono_group1"
    export london_casecontrol_group2="london_chr${i}.singleVar_MAC3_25Feb_CaseControl_nomono_group2_nostudycovar"
	export merck_casecontrol="Merck_070518_indqc_variantqc_chr${i}_7May2018_MAC3"
    
    export london_casecontrol_assoc_group1=${london_casecontrol_group1}.MetaScore.assoc.gz
    export london_casecontrol_assoc_group2=${london_casecontrol_group2}.MetaScore.assoc.gz
    export merck_casecontrol_assoc=${merck_casecontrol}.MetaScore.assoc.gz
    
    Rscript seqmineranno_london_group1_group2_merck.R
    
    # Export variable names again and index the required files with Tabix
    export anno_london_casecontrol_assoc_group1=anno_${london_casecontrol_assoc_group1}
    export anno_london_casecontrol_assoc_group2=anno_${london_casecontrol_assoc_group2}
    export anno_merck_casecontrol_assoc=anno_${merck_casecontrol_assoc}
    
    export london_casecontrol_covar_group1=${london_casecontrol_group1}.MetaCov.assoc.gz
    export london_casecontrol_covar_group2=${london_casecontrol_group2}.MetaCov.assoc.gz
    export merck_casecontrol_covar=${merck_casecontrol}.MetaCov.assoc.gz

    export genelist="refFlat_hg19_genelist_chr${i}_unique.txt"
    export i
    
    /usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group1}
    /usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_london_casecontrol_assoc_group2}
    /usr/bin/tabix -f -s 1 -b 2 -e 2 -S 1 ${anno_merck_casecontrol_assoc}

    mkdir ${resultsdir}/chr${i}
    nohup Rscript --max-ppsize=500000 raremetals_singlevar_london_meta_group1_group2_merck.R &> nohup_raremetals_singlevar_grouped_chr${i}.log &
done


# gcloud compute instances delete leademismeta2
