#!/usr/bin/env bash

#
# FusionCatcher v1.33
#


cd "$(dirname "$(realpath -s "$0")")"

# required
wget https://github.com/BenLangmead/bowtie/releases/download/v1.2.3/bowtie-1.2.3-linux-x86_64.zip -O bowtie-1.2.3-linux-x86_64.zip --no-check-certificate
unzip bowtie-1.2.3-linux-x86_64.zip
ln -s bowtie-1.2.3-linux-x86_64 bowtie

# required
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip -O bowtie2-2.3.5.1-linux-x86_64.zip --no-check-certificate
unzip bowtie2-2.3.5.1-linux-x86_64.zip
ln -s bowtie2-2.3.5.1-linux-x86_64 bowtie2

# optional (only in case that SRA is given as input instead of FASTQ)
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-centos_linux64.tar.gz -O sratoolkit.2.9.6-centos_linux64.tar.gz --no-check-certificate
tar --overwrite -xvzf sratoolkit.2.9.6-centos_linux64.tar.gz
ln -s sratoolkit.2.9.6-centos_linux64 sratoolkit

# optional (it is faster than gzip though?)
wget http://zlib.net/pigz/pigz-2.4.tar.gz -O pigz-2.4.tar.gz --no-check-certificate
tar --overwrite -xvzf pigz-2.4.tar.gz
make -C pigz-2.4
ln -s pigz-2.4 pigz

# required
mkdir liftover
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -O liftover/liftOver --no-check-certificate
chmod +x liftover/liftOver

# required (if there are Blat license issues FusionCatcher can still work without it)
mkdir blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -O blat/blat --no-check-certificate
chmod +x blat/blat

# required
mkdir fatotwobit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit -O fatotwobit/faToTwoBit --no-check-certificate
chmod +x fatotwobit/faToTwoBit

# required
wget http://github.com/ndaniel/seqtk/archive/1.2-r101c.tar.gz -O 1.2-r101c.tar.gz --no-check-certificate
tar --overwrite -xvzf 1.2-r101c.tar.gz -C .
make -C seqtk-1.2-r101c
chmod +x seqtk-1.2-r101c/seqtk
ln -s seqtk-1.2-r101c seqtk

# required
wget https://github.com/alexdobin/STAR/archive/2.7.2b.tar.gz -O 2.7.2b.tar.gz --no-check-certificate
tar --overwrite -xvzf 2.7.2b.tar.gz -C .
cp -f STAR-2.7.2b/bin/Linux_x86_64_static/STAR STAR-2.7.2b/source/STAR
ln -s STAR-2.7.2b star

# required
wget https://sourceforge.net/projects/bbmap/files/BBMap_38.44.tar.gz -O BBMap_38.44.tar.gz --no-check-certificate
tar --overwrite -xvzf BBMap_38.44.tar.gz -C .
mv bbmap BBMap_38.44
ln -s BBMap_38.44 bbmap
chmod +x bbmap/*.sh

# optional
wget  https://ftp.gnu.org/gnu/parallel/parallel-20201222.tar.bz2 -O parallel-20201222.tar.bz2 --no-check-certificate
tar --overwrite -xvjf parallel-20201222.tar.bz2 -C .
ln -s parallel-20201222 parallel
cd parallel
./configure
make
cd ..

# required
wget https://github.com/ndaniel/fastqtk/archive/v0.27.tar.gz -O v0.27.tar.gz --no-check-certificate
tar --overwrite -xvzf v0.27.tar.gz -C .
make -C fastqtk-0.27
chmod +x fastqtk-0.27/fastqtk
ln -s fastqtk-0.27 fastqtk

# optional (only in case that BAM is given as input instead of FASTQ)
mkdir picard
wget https://github.com/broadinstitute/picard/releases/download/2.21.2/picard.jar -O picard/picard.jar --no-check-certificate
chmod +x picard/picard.jar

# required (it installs biopython; if BioPython is allready installed then this can be skipped safely)
pip install biopython=1.73






