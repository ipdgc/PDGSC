# Install dependencies
sudo apt-get update
sudo apt-get install gcc make zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev -y

# Install R separately to get version > 3.3
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/'
sudo apt-get update
sudo apt-get install r-base -y

# And htslib
wget "https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2"
tar -zxvf rvtests-20170418-b5169f-linux64-static.tar.gz
tar jvxf htslib-1.4.tar.bz2
cd htslib-1.4
sudo ./configure
sudo make
cd ..

sudo apt-get remove tabix

# Install older version of Tabix (0.2.5 or older) to be able to index the summary/covariance files
wget http://launchpadlibrarian.net/74318164/tabix_0.2.5-1ubuntu1_amd64.deb
sudo dpkg -i tabix_0.2.5-1ubuntu1_amd64.deb
sudo apt-get -f install

# Now get rareMETALS2 and the refFlat gene file (don't need to download rareMETALS2 as have the working modified version)
#wget "http://genome.sph.umich.edu/w/images/b/b2/RareMETALS_6.8.tar.gz"
# wget "http://genome.sph.umich.edu/w/images/b/b7/RareMETALS2_0.1.tar.gz"
wget "http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat_hg19.txt.gz"

# Install R dependencies for rareMETALS2
sudo su - -c "R -e \"install.packages('seqminer', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('mvtnorm', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('CompQuadForm', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('getopt', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('data.table', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('curl', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('devtools', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('stringr', repos='http://cran.rstudio.com/')\""

# Make changes in rareMETALS2 and repackage:
# tar -zxvf rareMETALS2_0.1.tar.gz
# rm rareMETALS2_0.1.tar.gz
# R CMD build rareMETALS2

# Now install the new rareMETALS2:
sudo R CMD INSTALL rareMETALS2_0.1_modified0141017.tar.gz
sudo R CMD INSTALL rareMETALS_6.8_modified041017.tar.gz

# unzip the refFlat file
gunzip -c refFlat_hg19.txt.gz > refFlat_hg19.txt
# bgzip -d -c refFlat_hg19.txt.gz > refFlat_hg19.txt

# awk '$3 == "chr1" {print $1"\t"$2"\t"$3}' refFlat_hg19.txt > refFlat_hg19_genelist_chr1.txt

# Extract genelist per chromosome into a separate file:
awk '$3 == "chr1" || $3 ~ /chr1_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr1.txt
awk '$3 == "chr2" || $3 ~ /chr2_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr2.txt
awk '$3 == "chr3" || $3 ~ /chr3_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr3.txt
awk '$3 == "chr4" || $3 ~ /chr4_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr4.txt
awk '$3 == "chr5" || $3 ~ /chr5_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr5.txt
awk '$3 == "chr6" || $3 ~ /chr6_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr6.txt
awk '$3 == "chr7" || $3 ~ /chr7_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr7.txt
awk '$3 == "chr8" || $3 ~ /chr8_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr8.txt
awk '$3 == "chr9" || $3 ~ /chr9_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr9.txt
awk '$3 == "chr10" || $3 ~ /chr10_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr10.txt
awk '$3 == "chr11" || $3 ~ /chr11_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr11.txt
awk '$3 == "chr12" || $3 ~ /chr12_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr12.txt
awk '$3 == "chr13" || $3 ~ /chr13_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr13.txt
awk '$3 == "chr14" || $3 ~ /chr14_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr14.txt
awk '$3 == "chr15" || $3 ~ /chr15_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr15.txt
awk '$3 == "chr16" || $3 ~ /chr16_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr16.txt
awk '$3 == "chr17" || $3 ~ /chr17_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr17.txt
awk '$3 == "chr18" || $3 ~ /chr18_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr18.txt
awk '$3 == "chr19" || $3 ~ /chr19_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr19.txt
awk '$3 == "chr20" || $3 ~ /chr20_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr20.txt
awk '$3 == "chr21" || $3 ~ /chr21_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr21.txt
awk '$3 == "chr22" || $3 ~ /chr22_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr22.txt
awk '$3 == "chrX" || $3 ~ /chrX_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chrX.txt

# Remove duplicated gene names from the gene list files so left with just unique names. Also create a text file with gene numbers per gene as a sanity check...
rm genenumbers.txt
for i in {1..22} X
do
	sort -u refFlat_hg19_genelist_chr${i}.txt > refFlat_hg19_genelist_chr${i}_unique.txt
	wc -l refFlat_hg19_genelist_chr${i}_unique.txt >> genenumbers.txt
done

# Download annotation resouurces for seqminer annotator
Rscript seqminer_download.R

### READ IN: ANNOTYPE,STUDYTYPE,TEST,MAF_CUTOFF
wget http://launchpadlibrarian.net/74318164/tabix_0.2.5-1ubuntu1_amd64.deb
sudo dpkg -i tabix_0.2.5-1ubuntu1_amd64.deb
sudo apt-get -f install



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
