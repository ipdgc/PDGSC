##### Downloads & installs most commonly used tools and resources in the analysis pipeline

# Install some common utilities
sudo apt-get update
sudo apt-get install gcc make zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev autoconf g++ pkg-config unzip libxml2-dev libncurses5-dev libncursesw5-dev default-jre python3 python -y

# Install R separately to ensure version > 3.3
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/'
sudo apt-get update
sudo apt-get install r-base -y

# Install useful R packages
sudo bash ./pdgsc/scripts/shell/install_r_packages.sh

# Now install the new rareMETALS2/rareMETALS:
sudo R CMD INSTALL ./pdgsc/tools/rareMETALS2_0.1_modified.tar.gz
sudo R CMD INSTALL ./pdgsc/tools/rareMETALS_6.8_modified.tar.gz

# Download the hg19 refFlat file
wget "http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat_hg19.txt.gz"
gunzip -c refFlat_hg19.txt.gz > refFlat_hg19.txt

# Extract genelist per chromosome into a separate file:
# Also remove duplicated gene names from the gene list files so left with just unique names. Also create a text file with gene numbers per gene as a sanity check...
rm genenumbers.txt
for CHROM in {1..22} X
do
    CHROMOSOME=chr${CHROM}
    awk -v chrom="$CHROMOSOME" '$3 == chrom || $3 ~ /chrom_./ {print $1}' refFlat_hg19.txt > refFlat_hg19_genelist_chr${CHROM}.txt
    sort -u refFlat_hg19_genelist_chr${CHROM}.txt > refFlat_hg19_genelist_chr${CHROM}_unique.txt
    wc -l refFlat_hg19_genelist_chr${CHROM}_unique.txt >> genenumbers.txt
done

# Download annotation resources for seqminer annotator
Rscript ./pdgsc/scripts/r/seqminer_download.R

# Download exclusion regions for PCA
wget "https://github.com/gabraham/flashpca/blob/master/exclusion_regions_hg19.txt"

# Download & install Rvtest
wget https://github.com/zhanxw/rvtests/releases/download/v2.0.5/rvtests-20170613-444e2f-linux64-static.tar.gz
mkdir rvtests_june
tar zxvf rvtests-20170613-444e2f-linux64-static.tar.gz -C rvtests_june

### Install Vcftools, Bcftools, Htslib, Samtools etc:
wget "http://www.zlib.net/zlib-1.2.11.tar.gz"
tar -zxvf zlib-1.2.11.tar.gz
cd zlib-1.2.11 && sudo ./configure && sudo make && sudo make install && cd ..
git clone "https://github.com/vcftools/vcftools.git"
cd vcftools && sudo ./autogen.sh && sudo ./configure && sudo make && sudo make install && cd ..
git clone --branch=develop git://github.com/samtools/bcftools.git
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/samtools.git
cd htslib && sudo make && sudo make install && cd ..
cd bcftools && sudo make && sudo make install && cd ..
cd samtools && sudo make && sudo make install && cd ..

# Install ANNO
git clone "https://github.com/zhanxw/anno.git"
cd anno && sudo make && cd resources && ./download.sh && cd ../..

# Install Plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200121.zip
unzip plink_linux_x86_64*.zip

# Install Plinkseq
wget http://psychgen.u.hpc.mssm.edu/plinkseq_downloads/plinkseq-x86_64-latest.zip
unzip plinkseq-x86_64-latest.zip

### Download GATK & PICARD
wget "https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip"
wget "https://github.com/broadinstitute/picard/releases/download/2.16.0/picard.jar"
wget "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
wget "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai"
gunzip human_g1k_v37.fasta.gz
unzip gatk-4.1.4.1.zip
./gatk-4.1.4.1/gatk CreateSequenceDictionary -R human_g1k_v37.fasta -O human_g1k_v37.dict

### Download and install snpEff
wget "https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip"
unzip snpEff_latest_core.zip

# Install older version of Tabix (0.2.5 or older) to be able to index the summary/covariance files
sudo apt-get remove tabix -y
sudo rm -r /usr/local/bin/tabix
wget http://launchpadlibrarian.net/74318164/tabix_0.2.5-1ubuntu1_amd64.deb
sudo apt install ./tabix_0.2.5-1ubuntu1_amd64.deb
