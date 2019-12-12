
FusionCatcher
=============

Finder of somatic fusion-genes in RNA-seq data.


Download / Install / Update / Upgrade [FusionCatcher](http://github.com/ndaniel/fusioncatcher)
----------------------------------------------------------------------------------------------

Use this one-line command:

```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download
```

If one wants to have all the questions asked by boostrap.py answered automatically with yes then add `-y` to the 
command above. For more installing options, see:

```
bootstrap.py --help
```

On Ubuntu Linux running this command before installing FusionCatcher using `bootstrap.py` would help making the installation process smoother:

```
sudo apt-get install wget gawk gcc g++ make cmake automake curl unzip zip bzip2 tar gzip pigz parallel build-essential libncurses5-dev libc6-dev zlib1g zlib1g-dev libtbb-dev libtbb2 python python-dev python-numpy python-biopython python-xlrd python-openpyxl default-jdk
```

FusionCatcher can be installed also using `conda`, as follows:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n fusioncatcher fusioncatcher
source activate fusioncatcher
download-human-db.sh
```

FusionCatcher can be installed also from GitHub, as follows:
```
git clone https://github.com/ndaniel/fusioncatcher
cd fusioncatcher/tools/
./install_tools.sh
cd ../data
./download-human-db.sh
```
NOTE: Here it is assumed that Python 2.7.x, BioPython (>v1.5), and Java Runtime 
Environment 1.8 are already installed.



Description
-----------
FusionCatcher searches for novel/known somatic fusion genes, translocations, and
chimeras in RNA-seq data (paired-end or single-end reads from Illumina NGS platforms 
like Solexa/HiSeq/NextSeq/MiSeq/MiniSeq) from diseased samples.

The aims of FusionCatcher are:
 * very good detection rate for finding candidate somatic fusion
   genes (see somatic mutations; using a matched normal sample is
   optional; several databases of known fusion genes found in healthy
   samples are used as a list of known false positives; biological
   knowledge is used, like for example gene fusion between a gene and
   its pseudogene is filtered out),
 * very good RT-PCR validation rate of found candidate somatic fusion
   genes (this is very important for us),
 * very good detection of challenging fusion genes, like for example 
   IGH fusions, CIC fusions, DUX4 fusions, CRLF2 fusions, TCF3 fusions, etc.
 * very easy to use (i.e. no a priori knowledge of bioinformatic
   databases and bioinformatics is needed in order to run FusionCatcher BUT
   Linux/Unix knowledge is needed; it allows a very high level of control
   for expert users),
 * to be as automatic as possible (i.e. the FusionCatcher will choose
   automatically the best parameters in order to find candidate somatic
   fusion genes, e.g. finding automatically the adapters, quality trimming
   of reads, building the exon-exon junctions automatically based on the
   length of the reads given as input, etc. while giving also full control
   to expert users) while providing the best possible detection rate for
   finding somatic fusion genes (with a very low rate of false positives
   but a very good precision).


Manual
------
A detailed manual is available [here](doc/manual.md).


Forum
-----
A forum for FusionCatcher is available at 
[Google Groups](http://groups.google.com/d/forum/fusioncatcher).


Release history
---------------
A complete release history can be found [here](NEWS).


Official releases
-----------------
Old releases and the latest official release of FusionCatcher are on https://sourceforge.net/projects/fusioncatcher/files/


Citing
------
D. Nicorici, M. Satalan, H. Edgren, S. Kangaspeska, A. Murumagi, O. Kallioniemi,
S. Virtanen, O. Kilkku, FusionCatcher – a tool for finding somatic fusion genes
in paired-end RNA-sequencing data, bioRxiv, Nov. 2014, 
[DOI:10.1101/011650](http://dx.doi.org/10.1101/011650)


