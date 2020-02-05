##### Installs most commonly used tools in the analysis pipeline

sudo apt-get install gcc make zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev autoconf g++ pkg-config -y

### download plink
wget "https://www.cog-genomics.org/static/bin/plink170113/plink_linux_x86_64.zip"
sudo apt install unzip
unzip plink_linux_x86_64.zip

### download zlib (for vcftools)
wget "http://www.zlib.net/zlib-1.2.11.tar.gz"
tar -zxvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
sudo ./configure
sudo make
sudo make install
cd ..

### download vcftools
git clone "https://github.com/vcftools/vcftools.git"
cd vcftools
sudo ./autogen.sh
sudo ./configure
sudo make
sudo make install
cd ..

