### Define variables:
export BUCKET= # bucket location
export ORIGINAL_VCF= # name of the original VCF
export DIR_VCF_BY_CHROM= # directory within bucket for the original VCF split by chromosome and gzipped
export DIR_VCF_CHRPOS_VARIANT_ID= # directory within bucket for the original VCF split by chromosome and gzipped, with variant ID column replaced to be in chr:pos format (instead of rsid)
export VCF_STEM=${ORIGINAL_VCF//.recode.vcf/}
export SAMPLE_INFO= # name of the sample sheet

### copy the VCFs from the bucket into the VM
gsutil -m cp ${BUCKET}/${ORIGINAL_VCF} .
gsutil -m cp ${BUCKET}/${DIR_VCF_BY_CHROM}/${VCF_STEM}_chr*.recode.vcf* .
gsutil -m cp ${BUCKET}/${DIR_VCF_CHRPOS_VARIANT_ID}/${VCF_STEM}_chr*_variantid.vcf* .
gsutil -m cp ${BUCKET}/${SAMPLE_INFO} .

### zip and index the VCF
bgzip -c ${ORIGINAL_VCF} > ${ORIGINAL_VCF}.gz
tabix -p vcf ${ORIGINAL_VCF}.gz