### Define variables:
export BUCKET="bucket location"
export VCF_NAME="VCF filename stem"
export SAMPLE_INFO="updatedSampleInfoJanuary23rd2017.tab"
### copy the VCF from the bucket into the VM
gsutil -m cp -r ${BUCKET}/${VCF_NAME}.vcf .
gsutil -m cp -r ${BUCKET}/${SAMPLE_INFO} .

### install required tools
sudo bash install_tools.sh

### START QC PIPELINE ###
### Filter for call rate, MAF and HWE
./plink --vcf ${VCF_NAME}.vcf --double-id --maf 0.01 --geno 0.10 --make-bed --out ${VCF_NAME}_maf01_geno10

### Update the phenotypes in the .fam file
export FAM_FILE=${VCF_NAME}.fam
export UPDATED_FAM_FILE=${VCF_NAME}_newfam.fam
export MISSING_CLINICAL=${VCF_NAME}_missingpheno.txt

Rscript updatepheno.R

cp ${UPDATED_FAM_FILE} ${VCF_NAME}.fam

# Download and unzip Hapmap data
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/relationships_w_pops_041510.txt"
gunzip hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz
gunzip hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz

./plink --file hapmap3_r3_b36_fwd.consensus.qc.poly --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap3_r3_b36_fwd.consensus.qc.poly_maf01 --make-bed --maf 0.01
awk '{print $2}' hapmap3_r3_b36_fwd.consensus.qc.poly_maf01.bim > hapmap_snp_list
awk '{print $2}' pdgsc_joint.recode_maf01_geno10.bim > pdgsc_snp_list_maf01_geno10

./plink --bfile pdgsc_joint.recode_maf01_geno10 --out pdgsc_joint.recode_maf01_geno10_splitsex --make-bed --split-x b37
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap3_r3_b36_fwd.consensus.qc.poly_splitsex --make-bed --split-x b37

./plink --bfile pdgsc_joint.recode_maf01_geno10 --extract hapmap_snp_list --make-bed --out pdgsc_joint.recode_maf01_geno10_hapmapsnps
./plink --bfile hapmap3_r3_b36_fwd.consensus.qc.poly_maf01 --extract pdgsc_snp_list_maf01_geno10 --make-bed --out hapmap3_r3_b36_fwd.consensus.qc.poly_pdgscsnpsmaf01


