# Table of Contents
1. [Introduction](#1---introduction)
2. [Hardware requirements and dependencies](#2---hardware-requirements-and-dependencies)
3. [FusionCatcher in scientific articles](#3---fusioncatcher-in-scientific-articles)
4. [Installation and usage examples](#4---installation-and-usage-examples)
5. [Quick installation](#5---quick-installation)
6. [Usage](#6---usage)
7. [Aligners](#7---aligners)
8. [Command line options](#8---command-line-options)
9. [Methods](#9---methods)
10. [Comparisons to other tools](#10---comparisons-to-other-tools)
11. [License](#11---license)
12. [Citing](#12---citing)
13. [Reporting bugs](#13---reporting-bugs)

---


# 1 - INTRODUCTION

*FusionCatcher* searchers for **somatic** novel/known **fusion genes**, **translocations** and/or **chimeras** in RNA-seq data (stranded/unstranded  **paired-end/single-end** reads FASTQ files produced by Illumina next-generation sequencing platforms like Illumina Solexa/`HiSeq`/`NextSeq`/`MiSeq`/`MiniSeq`) from _diseased_ samples.

The aims of *FusionCatcher* are:
  * very good detection rate for finding candidate **somatic fusion genes** (see [somatic mutations](http://en.wikipedia.org/wiki/Mutation#Somatic_mutations); using a matched **normal** sample is optional; several databases of known fusion genes found in healthy samples are used as a list of known false positives; biological knowledge is used, like for example gene fusion between a gene and its pseudogene is filtered out),
  * very good RT-PCR validation rate of found candidate somatic fusion genes (this is very important for us),
  * very easy to use (i.e. no _a priori_ knowledge of bioinformatic databases and bioinformatics is needed in order to run *FusionCatcher* BUT Linux/Unix knowledge is needed; it allows a very high level of control for expert users),
  * to be as automatic as possible (i.e. the *FusionCatcher* will choose automatically the best parameters in order to find candidate somatic fusion genes, e.g. finding automatically the adapters, quality trimming of reads, building the exon-exon junctions automatically based on the length of the reads given as input, etc. while giving also full control to expert users) while providing the best possible detection rate for finding somatic fusion genes (with a very low rate of false positives but a very good sensitivity).

*FusionCatcher* supports:
  * as input FASTQ and/or SRA file types (paired-end reads from stranded or strand-specific experiments, single-end reads when they are longer than 130bp),
  * five different methods (using Bowtie aligner and optionally BLAT, STAR, BOWTIE2, BWA aligners) for finding new fusion genes **BUT** by default only Bowtie, Blat, and STAR aligners will be used,
  * several eukaryotic organisms (which are in [Ensembl database](http://www.ensembl.org/index.html)), like for example, human, rat, mouse, dog, etc.

---


# 2 - HARDWARE REQUIREMENTS AND DEPENDENCIES

For running *FusionCatcher* it is needed a computer with:
  * 64-bit `*NIX` environment
  * minimum 24 GB of RAM (in many cases it might work even with 16GB of RAM for very small input FASTQ files in order of megabytes)
  * 1 CPU (minimum)
  * ~700 GB temporary disk space (needed just for temporary files)


## 2.1 - Required dependencies
  * **Linux/Unix** 64-bit (e.g. Ubuntu version 12.04/14.04 or newer)
  * **Python** version 2.7.6 (>=2.6.0 and < 3.0 is fine)
  * **BioPython** version 1.66 (>=1.50 is fine)
  * **Bowtie** 64-bit version 1.1.2 http://bowtie-bio.sourceforge.net/index.shtml (will be installed by `boostrap.py`)
  * **SeqTK** version 1.0-r82b-dirty  http://github.com/ndaniel/seqtk (will be installed by `boostrap.py`)
  * organism specific  data from [Ensembl](http://www.ensembl.org) database release 84 (all downloading and the necessary building process is handled automatically by the included/provided tool `fusioncatcher-build` and therefore no knowledge of Ensembl database or other databases is needed) (will be installed by `boostrap.py`)
  * **STAR** version 2.5.1b https://github.com/alexdobin/STAR . Executables are available at http://github.com/alexdobin/STAR/releases (will be installed by `boostrap.py`)
  * **BOWTIE2** version 2.2.9 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml (will be installed by `boostrap.py`)
  * **BWA** version 0.7.12 http://sourceforge.net/projects/bio-bwa/ (will be installed by `boostrap.py`)

## 2.2 - Optional dependencies

### 2.2.1 - Strongly recommended
These are expected by default to be installed but their use can be disabled by using the command line option '--skip-blat'.
  * **BLAT** version 0.35 http://users.soe.ucsc.edu/~kent/src/ . Executables are available at http://hgdownload.cse.ucsc.edu/admin/exe/ . Please, check the license to see if it allows you to run/use it! This is needed by *FusionCatcher* (hint: if you are a non-profit organization you should be fine) (will be installed by `boostrap.py`)
  * **faToTwoBit** http://users.soe.ucsc.edu/~kent/src/ . Executables are available at http://hgdownload.cse.ucsc.edu/admin/exe/ . Please, check the license to see if it allows you to run/use it! This is needed by *FusionCatcher* and `fusioncatcher-build` if one plans to use BLAT as a second (optional) alternative method for finding fusion genes! (required also by option `--blat-visualization`) (will be installed by `boostrap.py`)

Note: If one does not want to install BLAT (whilst installing *FusionCatcher* automatically thru `bootstrap.py`) and also not to use BLAT with *FusionCatcher* then using command line `-k` option of `bootstrap.py` will do that.

### 2.2.2 - Nice to have
  * **Velvet** (de novo assembler) version 1.2.10 http://www.ebi.ac.uk/~zerbino/velvet/ . This is needed if one plans to do _de novo_ assembly of the reads which support the candidate fusion genes. (required by option `--assembly` of *FusionCatcher*) (will be installed by `boostrap.py`)
  * **fastq-dump** version 2.6.2 from NCBI SRA Toolkit http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software . This is needed by *FusionCatcher* if one plans to use as input SRA files (will be installed by `boostrap.py`)
  * Python library **openpyxl** version 1.5.6 http://pypi.python.org/pypi/openpyxl (other versions might work but have not been tested). It is needed by `fusioncatcher-build` for parsing the [ConjoinG](http://metasystems.riken.jp/conjoing/) database.
  * Python library **xlrd** version 0.6.1 http://pypi.python.org/pypi/xlrd (other versions might work but have not been tested). It is needed by `fusioncatcher-build` for parsing the [ChimerDB](http://ercsb.ewha.ac.kr/FusionGene/) database.
  * **coreutils** version 8.25  http://ftp.gnu.org/gnu/coreutils for a newer version of SORT command which allows the use of several CPUs in parallel, that is '--parallel'  command line options (other older versions might also support this!) (will be installed by `boostrap.py`)
  * **pigz** version 2.3.1 http://zlib.net/pigz/ for using GZIP on several CPUs in parallel (other older versions might support this) (will be installed by `boostrap.py`)
  * **SAMTools** version 1.19 http://www.htslib.org/ (will be installed by `boostrap.py`)
  * **Picard tools** version 2.2.2 http://broadinstitute.github.io/picard/ (will be installed by `boostrap.py`)


## 2.3 - Genomic Databases
These are used (downloaded and parsed) automatically by `boostrap.py` of *FusionCatcher*:
  * **ENSEMBL** database http://www.ensembl.org/ (required)
  * **UCSC** database http://hgdownload.cse.ucsc.edu/downloads.html#human (required)
  * **RefSeq** database (thru **UCSC** database) (required)
  * **Viruses/bacteria/phages** genomes database (from the NCBI database) [ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/](ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) (required)
  * **COSMIC** database http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/ (optional)
  * **TICdb** database http://www.unav.es/genetica/TICdb/ (optional)
  * **ChimerDB 2.0** database (literature-based annotation) http://ercsb.ewha.ac.kr/FusionGene/ (optional)
  * **Cancer Genome Project** **(CGP)** translocations database http://www.sanger.ac.uk/genetics/CGP/Census/ (optional)
  * **ConjoinG** database http://metasystems.riken.jp/conjoing/ (optional)
  * **CACG** conjoined genes database http://cgc.kribb.re.kr/map/ (optional)
  * **DGD** database http://dgd.genouest.org/ (optional)
  * **Illumina BodyMap2 RNA-seq** database http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/

**NOTES**:
  * **ENSEMBL** database is used for finding novel/known fusions genes
  * **COSMIC**, **TICdb**, **ChimerDB**, **Cancer Genome Project**, **ConjoinG**, and manual curated fusion gene database are indexed and used further for annotating/labeling the found fusion genes for an easier visualization of **novel** genes (i.e. not published yet) found by *FusionCatcher*. For more information how this is used see [Tables 1,2,3](#output-data).
  * *FusionCatcher* can work just fine and is able to find fusion genes without any of the optional dependencies/tools/programs!
  * if **BLAT** is not installed (or one does not want to use it) please use option '--skip-blat' in order to let know *FusionCatcher* that it should not use it! Also, specifying that BLAT should be be used can be done by editing manually the line `aligners = blat,star` from file `fusioncatcher/etc/configuration.cfg` (that is remove the `blat` string).


---

# 3 - FusionCatcher in scientific articles

*FusionCatcher* has been used for finding novel and known fusion genes in the following articles:
  * S. Kangaspeska, S. Hultsch, H. Edgren, D. Nicorici, A. Murumägi, O.P. Kallioniemi, **Reanalysis of RNA-sequencing data reveals several additional fusion genes with multiple isoforms**, PLOS One, Oct. 2012. http://dx.doi.org/10.1371/journal.pone.0048745
  * H. Edgren, A. Murumagi, S. Kangaspeska, D. Nicorici, V. Hongisto, K. Kleivi, I.H. Rye, S. Nyberg, M. Wolf, A.L. Borresen-Dale, O.P. Kallioniemi, **Identification of fusion genes in breast cancer by paired-end RNA-sequencing**, Genome Biology, Vol. 12, Jan. 2011. http://dx.doi.org/10.1186/gb-2011-12-1-r6
  * JN. Honeyman, EP. Simon, N. Robine, R. Chiaroni-Clarke, DG. Darcy, I. Isabel, P. Lim, CE. Gleason, JM. Murphy, BR. Rosenberg, L. Teegan, CN. Takacs, S. Botero, R. Belote, S. Germer, A-K. Emde, V. Vacic, U. Bhanot, MP. LaQuaglia, and S.M. Simon, **Detection of a Recurrent DNAJB1-PRKACA Chimeric Transcript in Fibrolamellar Hepatocellular Carcinoma**, Science 343 (6174), Feb. 2014, pp. 1010-1014, http://dx.doi.org/10.1126/science.1249484
  * T. Pietsch, I. Wohlers, T. Goschzik, V. Dreschmann, D. Denkhaus, E. Dorner, S. Rahmann, L. Klein-Hitpass, **Supratentorial ependymomas of childhood carry C11orf95–RELA fusions leading to pathological activation of the NF-kB signaling pathway**, Acta Neuropathologica 127(4), Apr. 2014, pp. 609-611. http://dx.doi.org/10.1007/s00401-014-1264-4
  * M. Jimbo, K.E. Knudsen, J.R. Brody, **Fusing Transcriptomics to Progressive Prostate Cancer**, The American Journal of Pathology, 2014, http://dx.doi.org/10.1016/j.ajpath.2014.08.001
  * Y.P. Yu, Y. Ding, Z. Chen, S. Liu, A. Michalopoulos, R. Chen, Z. Gulzar, B. Yang, K.M. Cieply, A. Luvison, B.G. Ren, J.D. Brooks, D. Jarrard, J.B. Nelson. G.K. Michalopoulos, G.C. Tseng, J.H. Luo, **Novel fusion transcripts associate with progressive prostate cancer**,  The American Journal of Pathology, 2014, http://dx.doi.org/10.1016/j.ajpath.2014.06.025
  * C.M Lindqvist, J. Nordlund, D. Ekman, A. Johansson, B.T. Moghadam, A. Raine, E. Overnas, J. Dahlberg, P. Wahlberg, N. Henriksson, J. Abrahamsson, B.M. Frost, D. Grander, M. Heyman, Rolf Larsson, J. Palle, S. Soderhall, E. Forestier, G. Lonnerholm, A.C. Syvanen, E.C. Berglund, **The Mutational Landscape in Pediatric Acute Lymphoblastic Leukemia Deciphered by Whole Genome Sequencing**, Human Mutation, 2014, http://dx.doi.org/10.1002/humu.22719
  * I. Panagopoulos, L. Gorunova, B. Davidson, Sverre Heim, **Novel TNS3-MAP3K3 and ZFPM2-ELF5 fusion genes identified by RNA sequencing in multicystic mesothelioma with t(7;17)(p12;q23) and t(8;11)(q23;p13)**, Cancer Letters, Dec. 2014, http://dx.doi.org/10.1016/j.canlet.2014.12.002
  * J.C Lee, Y.M. Jeng, S.Y. Su, C.T Wu, K.S. Tsai, C.H. Lee, C.Y. Lin, J.M. Carter, J. W. Huang, S.H. Chen, S.R. Shih, A. Marino-Enriquez, C.C. Chen, A.L. Folpe, Y.L. Chang, C.W. Liang, **Identification of a novel FN1–FGFR1 genetic fusion as a frequent event in phosphaturic mesenchymal tumour**, The journal of Pathology, Jan. 2015, http://dx.doi.org/10.1002/path.4465
  * J. Nordlund, C.L. Backlin, V. Zachariadis, L. Cavelier, J. Dahlberg, I. Ofverholm, G. Barbany, A. Nordgren, E. Overnas, J. Abrahamsson, T. Flaegstad, M.M. Heyman, O.G. Jonsson, J. Kanerva, R. Larsson, J. Palle, K. Schmiegelow, M.G. Gustafsson, G. Lonnerholm, E. Forestier, A.C. Syvanen, **DNA methylation-based subtype prediction for pediatric acute lymphoblastic leukemia**, Clinical Epigenetics, 7:11, Feb. 2015,  http://dx.doi.org/10.1186/s13148-014-0039-z
  * J.H. Luo, S. Liu, Z.H. Zuo, R. Chen, G.C. Tseng, Y.P. Yu, **Discovery and Classification of Fusion Transcripts in Prostate Cancer and Normal Prostate Tissue**, The American Journal of Pathology, May 2015, http://dx.doi.org/10.1016/j.ajpath.2015.03.008
  * T. Meissner, K.M. Fisch, L. Gioia, **OncoRep: an n-of-1 reporting tool to support genome-guided treatment for breast cancer patients using RNA-sequencing**, BMC Medical Genomics, May 2015, http://dx.doi.org/10.1186/s12920-015-0095-z
  * S. Torkildsen, L. Gorunova, K. Beiske, G.E. Tjonnfjord, S. Heim, I. Panagopoulos, **Novel ZEB2-BCL11B Fusion Gene Identified by RNA-Sequencing in Acute Myeloid Leukemia with t(2;14)(q22;q32)**, PLOS One, July 2015, http://dx.doi.org/10.1371/journal.pone.0132736
  * M. Cieslik, R. Chugh, Y.M. Wu, M. Wu, C. Brennan, R. Lonigro, F. Su, R. Wang, J. Siddiqui, R. Mehra, X. Cao, D. Lucas, A.M. Chinnaiyan, D. Robinson, **The use of exome capture RNA-seq for highly degraded RNA with application to clinical cancer sequencing**, Genome Research, August 2015, http://dx.doi.org/10.1101/gr.189621.115
  * E.P. Simon, C.A. Freije, B.A. Farber, G. Lalazar, D.G. Darcy, J.N. Honeyman, R. Chiaroni-Clark, B.D. Dill, H. Molina, U.K. Bhanot, M.P. La Quaglia, B.R. Rosenberg, S.M. Simon, **Transcriptomic characterization of fibrolamellar hepatocellular carcinoma**, PNAS, October 2015, http://dx.doi.org/10.1073/pnas.1424894112 
  * Y. Marincevic-Zuniga, V. Zachariadis, L. Cavelier, A. Castor, G. Barbany, E. Forestier, L. Fogelstrand, M. Heyman, J. Abrahamsson, G. Lonnerholm, A. Nordgren, A.C. Syvanen, J. Nordlund, **PAX5-ESRRB is a recurrent fusion gene in B-cell precursor pediatric acute lymphoblastic leukemia**, Haematologica, October 2015, http://dx.doi.org/10.3324/haematol.2015.132332
  * M. Brenca, S. Rossi, M. Polano, D. Gasparotto, L. Zanatta, D. Racanelli, L. Valori, S. Lamon, A.P. Dei Tos, R. Maestro, **Transcriptome sequencing identifies ETV6-NTRK3 as a gene fusion involved in GIST**, The Journal of Pathology, November 2015, http://dx.doi.org/10.1002/path.4677
  * Roberts K.G., et al. **High Frequency and Poor Outcome of Ph-like Acute Lymphoblastic Leukemia in Adults**, Blood Journal, Vol. 126, December 2015, http://www.bloodjournal.org/content/126/23/2618
  * Kekeeva T., et al. **Novel fusion transcripts in bladder cancer identified by RNA-seq**, Cancer Letters, Feb. 2016, http://dx.doi.org/10.1016/j.canlet.2016.02.010
  * Panagopoulos I. et al. **Rare MLL-ELL fusion transcripts in childhood acute myeloid leukemia-association with young age and myeloid sarcomas?**, Experimental Hematology & Oncology, Vol. 5, March 2016, http://dx.doi.org/10.1186/s40164-016-0037-2 
  * Chang W. et al, **Multi-Dimensional ClinOmics for Precision Therapy of Children and Adolescent Young Adults with Relapsed and Refractory Cancer: A report from the Center for Cancer Research**, Clinical Cancer Research, March 2016, http://dx.doi.org/10.1158/1078-0432.CCR-15-2717
  * Spans L. et al., **Recurrent MALAT1-GLI1 oncogenic fusion and GLI1 upregulation define a subset of plexiform fibromyxoma**, The Journal of Pathology, April 2016, http://dx.doi.org/10.1002/path.4730
  * Micci F. et al., **Cytogenetic and Molecular Profile of Endometrial Stromal Sarcoma**, Genes Chromosomes & Cancer, May 2016, http://dx.doi.org/10.1002/gcc.22380
  * Olsen K.T. et al., **Novel Fusion Genes and Chimeric Transcripts in Ependymal Tumors**, Genes Chromosomes & Cancer, July 2016, http://dx.doi.org/10.1002/gcc.22392
  * Panagopoulos I. et al., **Recurrent fusion of the genes FN1 and ALK in gastrointestinal leiomyomas**, Modern Pathology, July 2016, http://dx.doi.org/10.1038/modpathol.2016.129
  * Barnes D.J. et al., **A germline mutation of CDKN2A and a novel RPLP1-C19MC fusion detected in a rare melanotic neuroectodermal tumor of infancy: a case report**, BMC Cancer, August 2016, http://dx.doi.org/10.1186/s12885-016-2669-3
  * Panagopoulos I. et al., **Gene fusions AHRR-NCOA2, NCOA2-ETV4, ETV4-AHRR, P4HA2-TBCK, and TBCK-P4HA2 resulting from the translocations t(5;8;17)(p15;q13;q21) and t(4;5)(q24;q31) in a soft tissue angiofibroma**, Oncology Reports, Sept. 2016, http://dx.doi.org/10.3892/or.2016.5096
  * Lang P.Y. et al., **ATR maintains chromosomal integrity during postnatal cerebellar neurogenesis and is required for medulloblastoma formation**, Development, Nov. 2016, http://dx.doi.org/10.1242/dev.139022
  * Gu Z. et al., **Genomic analyses identify recurrent MEF2D fusions in acute lymphoblastic leukaemia**, Nature Communications, Nov. 2016, http://dx.doi.org/10.1038/ncomms13331
  * Alaei-Mahabadi B. et al., **Global analysis of somatic structural genomic alterations and their impact on gene expression in diverse human cancers**, PNAS, Nov. 2016, http://dx.doi.org/10.1073/pnas.1606220113
  * Yap K.L. et al., **Diagnostic evaluation of RNA sequencing for the detection of genetic abnormalities associated with Ph-like acute lymphoblastic leukemia (ALL)**, Leukemia & Lymphoma, Nov. 2016, http://dx.doi.org/10.1080/10428194.2016.1219902
  
  
---

# 4 - INSTALLATION AND USAGE EXAMPLES

## 4.1 - Automatic installation {#automatic-installation}

This is an example of automatic installation of *FusionCatcher* (and it is installed here "~/fusioncatcher" if these are run in your home directory) and the required databases and indexes (which are downloaded instead of being built locally):
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download
```
where:
  * `wget http://sf.net/projects/fusioncatcher/files/bootstrap.py` downloads from internet the `bootstrap.py` which is the installation script (it is recommended to use the `boostrap.py` from `ttp://sf.net/projects/fusioncatcher/files/bootstrap.py` because it is more up to date)
  * `python bootstrap.py` runs using `python` the installation script `bootstrap.py` (here one may replace `python` with its own custom installation of python, like for example `/some/other/custom/python`)
  * `-t` installs the software tools (and their exact version) needed by *FusionCatcher*
  * `--download` forces the installation script `bootstrap.py` to download and install automatically also the databases needed by *FusionCatcher* (if this is not used the databases needed by *FusionCatcher* will not be installed and the user will have to build/install them manually later)

In case that there are several Python versions installed already then it is possible to point which one to use for installation and running *FusionCatcher*, as following (no required databases and indexes are installed automatically in this example):
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py

/some/other/python bootstrap.py
```

In case that one wants to install *FusionCatcher* here `/some/directory/fusioncatcher/`, then this shall be run (no required databases and indexes are installed automatically in this example):
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py

/your/favourite/python bootstrap.py --prefix=/some/directory/
```

In case that one wants to install *FusionCatcher* and download the databases directly and build locally the indexes, then this shall be run:
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py --build
```

This is an example of automatic installation of *FusionCatcher* and the required databases and indexes (which are downloaded instead of being built locally) while all the questions asked by the installation script are answered automatically with YES (WARNING: this might overwrite files/directories):
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py --download -y
```

In case that one has the admin/root rights then it is possible to install *FusionCatcher* as following (no required databases and indexes are installed automatically in this example):
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py
sudo python bootstrap.py
```

In case that one plans not to use at all BLAT with *FusionCatcher* then add the command line `-k` to the `boostrap.py`, as following:

```
python bootstrap.py -k
```

In case that you do not know which one to use from these examples, please use the first one!
Also, for more info about what options offered by `bootstrap.py`, please run
```
bootstrap.py --help
```

Please, do not forget to build/download the organism data after this is done running (please notice the last lines displayed by `bootstrap.py` after it finished running and execute the commands suggested there, e.g. use `download.sh`)!


## 4.2 - Manual installation {#manual-installation}

This is an example (or one of the many ways) for installing *FusionCatcher* on a **Ubuntu Linux 12.04/14.04 64-bit system** and the *FusionCatcher* and its dependencies are installed in `/apps`.

  * check that Python 2.6.X or 2.7.X is installed and working properly! If not then install it together with its development environment and other (probably) needed dependencies (required):
  ```
  sudo apt-get install build-essential
  sudo apt-get install python-dev
  ```
  
  and for [RedHat](http://www.redhat.com)/[CentOS](http://www.centos.org) this would be required
  ```
  sudo yum groupinstall "Development Tools"
  sudo yum install gcc
  sudo yum install ncurses-devel
  sudo yum install python-devel
  sudo yum install zlib-devel
  ```
  
  and for [OpenSUSE](http://www.opensuse.org) this would be required
  ```
  sudo zypper in --type pattern Basis-Devel
  sudo zypper in gcc
  sudo zypper in ncurses-devel
  sudo zypper in python-devel
  sudo zypper in zlib-devel
  ```
  
  * installing **BioPython** (required):
  ```
  sudo apt-get install python-numpy
  sudo apt-get install python-biopython
  ```
  
  * installing Python module **xlrd** (optional):
  ```
  sudo apt-get install python-xlrd
  ```
  
  * installing Python module **openpyxl** (optional):
  ```
  sudo apt-get install python-openpyxl
  ```
  
  * create the needed directories:
  ```
  mkdir -p /apps/fusioncatcher/tools
  mkdir -p /apps/fusioncatcher/data
  ```
  
  * installing **Bowtie** 64-bit version 1.1.2 (required)
  ```
  cd /apps/fusioncatcher/tools
  wget http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.1/bowtie-1.1.2-linux-x86_64.zip
  unzip bowtie-1.1.2-linux-x86_64.zip
  ln -s bowtie-1.1.2 bowtie
  ```
  
  * installing **Bowtie2** 64-bit version 2.2.9 (required)
  ```
  cd /apps/fusioncatcher/tools
  wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip
  unzip bowtie2-2.2.9-linux-x86_64.zip
  ln -s bowtie2-2.2.9-linux-x86_64 bowtie2
  ```
  
  * installing **BLAT** version 0.35 (optional; if **BLAT** is not installed please use option '--skip-blat' or remove the `blat` string from the `aligners =...` line of file `fusioncatcher/etc/configuration.cfg` in order to let know *FusionCatcher* that it should not use it)
  ```
  cd /apps/fusioncatcher/tools
  mkdir blat_v0.35
  cd blat_v0.35
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat
  chmod +x blat
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/faToTwoBit
  chmod +x faToTwoBit
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/liftOver
  chmod +x liftOver
  cd ..
  ln -s blat_v0.35 blat
  ```
  
  * installing **fastq-dump** version 2.6.2 from **NCBI SRA Toolkit** version 2.6.2 (optional but it is required in case that one wants to test the installation using the example given below)
  ```
  cd /apps/fusioncatcher/tools
  wget http://ftp-private.ncbi.nlm.nih.gov/sra/sdk/2.6.2/sratoolkit.2.6.2-centos_linux64.tar.gz
  tar zxvf sratoolkit.2.6.2-centos_linux64.tar.gz
  ln -s sratoolkit.2.6.2-centos_linux64 sratoolkit
  ```
  
  * installing **SeqTK** version 1.0-r82b (please note that this is a different development branch than the official development) (required)
  ```
  cd /apps/fusioncatcher/tools
  wget http://github.com/ndaniel/seqtk/archive/1.0-r82b.tar.gz -O 1.0-r82b.tar.gz
  tar zxvf 1.0-r82b.tar.gz
  cd seqtk-1.0-r82b
  make
  cd ..
  ln -s seqtk-1.0-r82b seqtk
  ```
  
  * installing **STAR** version 2.5.1b (required)
  ```
  cd /apps/fusioncatcher/tools
  wget http://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz -O 2.5.1b.tar.gz
  tar zxvf 2.5.1b.tar.gz
  cd 2.5.1b
  cd source
  rm -f STAR
  cp ../bin/Linux_x86_64_static/STAR .
  ```
  
  Try to run this command (if it fails please ignore the error messages and continue further; continue further either way)
  ```
  make
  ```
  
  and continue with
  ```
  cd ..
  ln -s 2.5.1b star
  ```
  
  * installing **Velvet** version 1.2.10 (optional)
  ```
  cd /apps/fusioncatcher/tools
  wget http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
  tar zxvf velvet_1.2.10.tgz
  cd velvet_1.2.10
  make
  cd ..
  ln -s velvet_1.2.10 velvet
  ```
  > *Note*: Velvet depends on zlib-dev which may be installed like this
  ```
  sudo apt-get install zlib-dev
  ```
  
  * installing **coreutils** version 8.25 for a newer version of **sort** command which allows use of several CPUs in parallel, that is the use of `--parallel` command line option (optional)
  ```
  cd /apps/fusioncatcher/tools
  wget http://ftp.gnu.org/gnu/coreutils/coreutils-8.25.tar.xz
  tar --xz -xvf coreutils-8.25.tar.xz
  cd coreutils-8.25
  ./configure
  make
  cd ..
  ln -s coreutils-8.25 coreutils
  ```
  
  * installing **pigz** version 1.2.10 (optional)
  ```
  cd /apps/fusioncatcher/tools
  wget http://zlib.net/pigz/pigz-2.3.1.tar.gz
  tar zxvf pigz-2.3.1.tar.gz
  cd pigz-2.3.1
  make
  cd ..
  ln -s pigz-2.3.1 pigz
  ```
  
  * installing **SAMTools** version 1.19 (optional)
  ```
  cd /apps/fusioncatcher/tools
  wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
  tar jxvf samtools-0.1.19.tar.bz2
  cd samtools-0.1.19
  make
  cd ..
  ln -s samtools-0.1.19 samtools
  ```
  
  * installing **Picard tools** version 2.2.2 (optional)
  ```
  cd /apps/fusioncatcher/tools
  wget http://github.com/broadinstitute/picard/releases/download/2.2.2/picard-tools-2.2.2.zip
  unzip picard-tools-2.2.2.zip
  ln -s picard-tools-2.2.2 picard
  ```
  
  * installing *FusionCatcher* version 0.99.6a (required)
  ```
  cd /apps/fusioncatcher
  wget http://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.6a.zip
  unzip fusioncatcher_v0.99.6a.zip
  cd fusioncatcher_v0.99.6a
  
  rm -rf ../bin
  rm -rf ../etc
  rm -rf ../doc
  rm -rf ../VERSION
  rm -rf ../NEWS
  rm -rf ../LICENSE
  rm -rf ../README.md
  rm -rf ../DEPENDENCIES
  
  ln -s $(pwd)/bin ../bin
  ln -s $(pwd)/etc ../etc
  ln -s $(pwd)/doc ../doc
  ln -s $(pwd)/test ../test
  ln -s $(pwd)/VERSION ../VERSION
  ln -s $(pwd)/NEWS ../NEWS
  ln -s $(pwd)/LICENSE ../LICENSE
  ln -s $(pwd)/README.md ../README.md
  ln -s $(pwd)/DEPENDENCIES ../DEPENDENCIES
  ```
  
  * specify the paths to the above tools such that *FusionCatcher* can find them. There are two choices.
   * *Choice A*: Edit the *FusionCatcher* configuration file `configuration.cfg` (type: `nano /apps/fusioncatcher/etc/configuration.cfg` at command line), which has priority over `PATH`,  and make sure that the *FusionCatcher*'s configuration file **'configuration.cfg'** looks like this:
   ```
   [paths]
   python = /usr/bin/
   data = /apps/fusioncatcher/data/current/
   scripts = /apps/fusioncatcher/bin/
   bowtie = /apps/fusioncatcher/tools/bowtie/
   blat = /apps/fusioncatcher/tools/blat/
   bowtie2 = /apps/fusioncatcher/tools/bowtie2/
   bwa = /apps/fusioncatcher/tools/bwa/
   star = /apps/fusioncatcher/tools/star/source/
   seqtk = /apps/fusioncatcher/tools/seqtk/
   velvet = /apps/fusioncatcher/tools/velvet/
   fatotwobit = /apps/fusioncatcher/tools/blat/
   liftover = /apps/fusioncatcher/tools/blat/
   sra = /apps/fusioncatcher/tools/sratoolkit/bin/
   numpy = /apps/fusioncatcher/tools/numpy/
   biopython = /apps/fusioncatcher/tools/biopython/
   xlrd = /apps/fusioncatcher/tools/xlrd/
   openpyxl = /apps/fusioncatcher/tools/openpyxl
   lzop = /apps/fusioncatcher/tools/lzop/src/
   coreutils = /apps/fusioncatcher/tools/coreutils/src/
   pigz = /apps/fusioncatcher/tools/pigz/
   samtools = /apps/fusioncatcher/tools/samtools/
   picard = /apps/fusioncatcher/tools/picard/
   parallel = /appsfusioncatcher/tools/paralell/src/
   pxz = /apps/fusioncatcher/tools/pxz/
   java = /usr/bin/
   [parameters]
   threads = 0
   aligners = blat,star
   [versions]
   fusioncatcher = 0.99.6a beta
   ```
   
   * *Choice B*: Add the paths for the needed tools to the `PATH` variable by editing, for example, the `.bashrc` file (type: `nano ~/.bashrc` at command line) and add the following lines at the end:
   ```
   export PATH=/apps/fusioncatcher/bin:$PATH
   export PATH=/apps/fusioncatcher/tools/bowtie:$PATH
   export PATH=/apps/fusioncatcher/tools/bowtie2:$PATH
   export PATH=/apps/fusioncatcher/tools/bwa:$PATH
   export PATH=/apps/fusioncatcher/tools/blat:$PATH
   export PATH=/apps/fusioncatcher/tools/star/source/:$PATH
   export PATH=/apps/fusioncatcher/tools/liftover:$PATH
   export PATH=/apps/fusioncatcher/tools/seqtk:$PATH
   export PATH=/apps/fusioncatcher/tools/sratoolkit/bin:$PATH
   export PATH=/apps/fusioncatcher/tools/velvet/:$PATH
   export PATH=/apps/fusioncatcher/tools/fatotwobit/:$PATH
   export PATH=/apps/fusioncatcher/tools/lzop/src/:$PATH
   export PATH=/apps/fusioncatcher/tools/coreutils/src/:$PATH
   export PATH=/apps/fusioncatcher/tools/pigz/:$PATH
   export PATH=/apps/fusioncatcher/tools/samtools/:$PATH
   export PATH=/apps/fusioncatcher/tools/picard/:$PATH
   ```
   
   > *Note 1*: If a different version of Python is used/needed by *FusionCatcher* than the standard `/usr/bin/env python` then also please make sure that that specific version of Python is added to the `PATH` variable by editing, for example, the `.bashrc` file (type: `nano ~/.bashrc` at command line) or add the following lines at the end:
   ```
   export PATH=/some/other/version/of/python:$PATH
   ```
   
   > *Note 2*: `fusioncatcher/etc/configuration.cfg` **has priority** over `$PATH`. 
   
   > *Note 3*: In some cases it might not be enough to change the Python's path in `.bashrc` file, like for example the case when *FusionCatcher* is run on a server which defaults to another Python than one used to install *FusionCatcher*. In this case it is required that one changes all the [shebangs](http://en.wikipedia.org/wiki/Shebang_(Unix)) of the all Python scripts which belong to *FusionCatcher*. In case that one uses the Python which has the following Python executable path `/some/other/python` than this can be done like this (it changes in place `/usr/bin/env python` into `/some/other/python` in all `/apps/fusioncatcher/bin/*.py`):
   ```
   sed -i 's/\/usr\/bin\/env\ python/\/some\/other\/python/g' /apps/fusioncatcher/bin/*.py
   ```
  
  * download/build the human build data from Ensembl database and other databases and build the necessary indexes and files (the latest release of Ensembl data is release 84 as May 2016 when this section was updated last time). There are two alternative ways to get the human **build data**. The recommended way is to use `fusioncatcher-build`.
   * Using direct download
   ```
   mkdir -p /apps/fusioncatcher/data
   cd /apps/fusioncatcher/data
   wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.aa
   wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.ab
   wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.ac
   wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.ad
   cat ensembl_v84.tar.gz.* | tar xz
   ln -s ensembl_v84 current
   ```
   
   * Using `fusioncatcher-build` -- It will takes several hours (e.g. 5-10 hours) and it depends highly on the bandwidth of your internet connection. One may find out what Ensembl database version is available at [www.ensembl.org] and what version has been downloaded by looking to the last three lines printed on the screen by `fusioncatcher-build`.
   ```
   mkdir -p /apps/fusioncatcher/data/ensembl_v84
   cd /apps/fusioncatcher/data/ensembl_v84
   /apps/fusioncatcher/bin/fusioncatcher-build -g homo_sapiens -o .
   cd ..
   ln -s ensembl_v84 current
   ```



## 4.3 - Semi-automatic installation

This is an example of semi-automatic installation of *FusionCatcher* (and it is installed here: `/some/server/apps/fusioncatcher`). This may be used when *FusionCatcher* should be installed on a computer without internet connection. Shortly, in this case all the software dependencies and indexes of databases need to be downloaded separately on another computer which has internet connection and from there they should be copied/moved to the computer without internet connection. Here are the steps for achieving these:

  * on computer A (which has internet connection):
   * create locally a folder named, for example, `fuscat`:
    ```
    mkdir fuscat
    cd fuscat
    ```
   * download `bootstrap.py`
    ```
    wget http://sf.net/projects/fusioncatcher/files/bootstrap.py
    ```
   * find out the dependencies needed to be downloaded and download them manually into folder `fuscat` (their URLs will be shown and the user needs to download them manually using wget or its favourite browser)
    ```
    python bootstrap.py --list-dependencies
    ...
    ```
   * copy/move (manually) the folder `fuscat` and all its content to computer B (which does not have internet connection) where one intends to install *FusionCatcher*
  * on the computer B (which does not have internet connection), where one intends to install *FusionCatcher*:
   * go to the folder `fuscat` and make sure that the downloaded files do **not** have their permissions set as executables (this might confuse `bootstrap.py`)
    ```
    cd fuscat
    chmod -x *
    ```
   * start the installing process of *FusionCatcher* using `bootstrap.py` (if one wishes to use another version of Python, like for example having the path `/some/other/python` then below please replace `python` with `/some/other/python`)
    ```
    python bootstrap.py --local .
    ```
   * for installing the pre-built index files for human, please run (or take a look for instructions to) `download.sh` (it should be in `bin` directory where *FusionCatcher* has been installed)
   * for installing the pre-built index files for other organisms than human please, use `fusioncatcher-build` according tot the manual

For more information regarding the installation settings and possibilities, run:
```
python bootstrap.py --help
```

## 4.4 - Testing installation

This test works only when human organism.

### 4.4.1 - Automatic

Here are the steps for testing the installation of *FusionCatcher* using human genome.
```
cd ~
/apps/fusioncatcher/test/test.sh
```
Afterwards a message will be shown at console if the installation test went fine or not.

### 4.4.1 - Manual

Here are the steps for testing the installation of *FusionCatcher* using human genome.

```
mkdir ~/test
cd ~/test

wget http://sourceforge.net/projects/fusioncatcher/files/test/reads_1.fq.gz
wget http://sourceforge.net/projects/fusioncatcher/files/test/reads_2.fq.gz

cd ..

/apps/fusioncatcher/bin/fusioncatcher \
-d /apps/fusioncatcher/data/current/ \
--input ~/test/ \
--output ~/test-results/
```

This should take around 5 minutes to run and the result file `~/test-results/final-list_candidates-fusion-genes.txt` should look like [this](http://sourceforge.net/projects/fusioncatcher/files/test/final-list_candidate-fusion-genes.txt).

This dataset contains a very small set of short reads covering 12 already known fusion genes from human tumor cell lines which have been RNA sequenced (for more see [here](http://sourceforge.net/projects/fusioncatcher/files/test/readme.txt)).

## 4.5 - Breast cancer cell line

This is an example of finding fusion genes in the BT474 cell line using the public available RNA-seq data (from SRA archive):
  * download the publicly available RNA-seq data for BT-474 tumor breast cell line published in article **H. Edgren, A. Murumagi, S. Kangaspeska, D. Nicorici, V. Hongisto, K. Kleivi, I.H. Rye, S. Nyberg, M. Wolf, A.L. Borresen-Dale, O.P. Kallioniemi, Identification of fusion genes in breast cancer by paired-end RNA-sequencing, Genome Biology, Vol. 12, Jan. 2011** http://genomebiology.com/2011/12/1/R6/ :
   ```
   mkdir -p ~/bt474
   cd ~/bt474
   wget http://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR064/SRR064438/SRR064438.sra
   wget http://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR064/SRR064439/SRR064439.sra
   ```
   
  * run *FusionCatcher* (it takes around ~2.5 hours):
   ```
   /apps/fusioncatcher/bin/fusioncatcher \
   -d /apps/fusioncatcher/data/current/ \
   -i ~/bt474/ \
   -o ~/bt474_fusions/
   ```
   
  * if the run was successful then there should be the (non-empty) files (for more information see [here](#output-data)):
   ```
   ~/bt474_fusions/final-list_candidate_fusion_genes.txt
   ~/bt474_fusions/summary_candidate_fusions.txt
   ~/bt474_fusions/viruses_bacteria_phages.txt
   ~/bt474_fusions/supporting-reads_gene-fusions_BOWTIE.zip
   ~/bt474_fusions/supporting-reads_gene-fusions_BLAT.zip
   ~/bt474_fusions/supporting-reads_gene-fusions_STAR.zip
   ~/bt474_fusions/info.txt
   ~/bt474_fusions/fusioncatcher.log
  ```
  
  and the file
  ```
  ~/bt474_fusions/final-list_candidate_fusion_genes.txt
  ```
  
  should contain almost all fusion genes which have been published here:
   * S. Kangaspeska, S. Hultsch, H. Edgren, D. Nicorici, A. Murumägi, O.P. Kallioniemi, Reanalysis of RNA-sequencing data reveals several additional fusion genes with multiple isoforms, PLOS One, Oct. 2012. http://dx.plos.org/10.1371/journal.pone.0048745
   * H. Edgren, A. Murumagi, S. Kangaspeska, D. Nicorici, V. Hongisto, K. Kleivi, I.H. Rye, S. Nyberg, M. Wolf, A.L. Borresen-Dale, O.P. Kallioniemi, Identification of fusion genes in breast cancer by paired-end RNA-sequencing, Genome Biology, Vol. 12, Jan. 2011. http://genomebiology.com/2011/12/1/R6

## 4.6 - Batch mode

This is an example of finding fusion genes in the [Illumina Body Map 2.0](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/) RNA-seq data which consists of 16 RNA samples from 16 different organs from healthy persons. For doing this the batch mode is used (where the input is [this file](http://sourceforge.net/projects/fusioncatcher/files/examples/illumina-bodymap2.txt)), as shown here:
```
mkdir -p ~/bodymap
cd ~/bodymap
wget http://sourceforge.net/projects/fusioncatcher/files/examples/illumina-bodymap2.txt
fusioncatcher-batch.py -i illumina-bodymap2.txt -o results
```
The input file for `fusioncatcher-batch.py` is a text tab-separated file with two columns and 16 lines (one line for each organ from Illumina Body Map 2.0). The first column contains the URLs for the input FASTQ files and the second column (which is optional) contains the name of the organ (which will be used to create a output directory later where the results will be). Therefore, *FusionCatcher* will be run automatically 16 times by the `fusioncatcher-batch.py`.

The fusion genes found in *Illumina Body Map 2.0* could be used later, for example, as a list of known false positives when looking for fusion genes in diseased/tumor samples.

## 4.7 - Matched normal sample

In case that there is available RNA-seq data from a tumor sample and its match normal sample then the somatic mode of *FusionCatcher* may be used.  By default *FusionCatcher* is using a background list of fusion genes which have been found previously in normal healthy samples (e.g. Illumina BodyMap2 , etc.).

For example, lets assume that (i) the BT-474 is the rumor sample from here, and (ii) the matched normal samples if the healthy breast sample from here.  In this case, in order to find the somatic fusion genes in the BT-474 (that are the fusion genes which are found in BT-474 and are not found in the healthy sample) *FusionCatcher* should be run as follows:

```
mkdir -p ~/bt474
cd ~/bt474
wget http://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR064/SRR064438/SRR064438.sra
wget http://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR064/SRR064439/SRR064439.sra


mkdir -p ~/healthy
cd ~/healthy
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR064/SRR064437/SRR064437.sra

cd ~

fusioncatcher.py
--input ~/bt474/
--normal ~/healthy
--output ~/results/
```

The somatic fusion genes for BT-474 will be found in `~/results/bt474/final-list_candidate_fusion_genes.txt` file. The fusion genes which marked as **matched-normal** (see column 3) in BT-474 (that is `~/results/bt474/final-list_candidate_fusion_genes.txt` file) have been found also in the healthy sample also and most likely they are not somatic.

In case that there are several tumor samples and their matched healthy samples then batch mode of *FusionCatcher* may be used, as follows:
```
fusioncatcher-batch.py
--input /some/path/tumor-file.txt
--normal /some/path/healthy-file.txt
--output /some/path/results/
```
where:
  * `/some/path/tumor-file.txt` is a text file containing on each line a path to FASTQ files belonging to the tumor cells (an example is [here](http://sourceforge.net/projects/fusioncatcher/files/examples/edgren.txt)),
  * `/some/path/healthy-file.txt` is a text file containing on each line a path to FASTQ files belonging to the healthy cells (an example is [here](http://sourceforge.net/projects/fusioncatcher/files/examples/illumina-bodymap2.txt)),
  * `/some/path/results` is the output directory where the results are placed.

## 4.8 - Edgren RNA-seq dataset

This is an example of finding fusion genes in the Edgren RNA-seq data (from SRA archive):
  * H. Edgren, A. Murumagi, S. Kangaspeska, D. Nicorici, V. Hongisto, K. Kleivi, I.H. Rye, S. Nyberg, M. Wolf, A.L. Borresen-Dale, O.P. Kallioniemi, **Identification of fusion genes in breast cancer by paired-end RNA-sequencing**, Genome Biology, Vol. 12, Jan. 2011. http://genomebiology.com/2011/12/1/R6
  * S. Kangaspeska, S. Hultsch, H. Edgren, D. Nicorici, A. Murumägi, O.P. Kallioniemi, **Reanalysis of RNA-sequencing data reveals several additional fusion genes with multiple isoforms**, PLOS One, Oct. 2012. http://dx.plos.org/10.1371/journal.pone.0048745
```
fusioncatcher-batch.py -i http://sourceforge.net/projects/fusioncatcher/files/examples/edgren.txt -o results/
```
NOTE: **DO NOT POOL** the samples from all these cell lines. **DO NOT** give at once all these SRA/FASTQ files as input to FusionCatcher! Run *FusionCatcher* separately for each cell line! It is ok the pool together the samples from the same cell line together (but still do not concatenate yourself the FASTQ files and let FusionCatcher do it for you)!


---


# 5 - QUICK INSTALLATION

## 5.1 - Getting executables

For a fully automatic installation (including the required indexes of databases) run:
```
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py && python bootstrap.py --download
```

In case of a manual installation, first please check that (i) the required dependencies are installed, and (ii) download the source files of *FusionCatcher*, like for example:
```
wget http://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.6a.zip 
unzip fusioncatcher_v0.99.6a.zip
```

For an example of:
  * fully automatic installation see [here](#automatic-installation),
  * manual installation see [here](#manual-installation), and
  * semi-automatic installation see  [here](#manual-installation).

## 5.2 - Organism's build data

First, it is needed to download data or build the necessary files/indexes for running the *FusionCatcher*. This process should be done once for every single organism or every time when the [Ensembl](http://www.ensembl.org) database is updated.


### 5.2.1 - Direct download of human build data

Here, in this example, the necessary data is downloaded and necessary files/indexes for the **human genome** are downloaded in the directory `/some/human/data/ensembl_v84/` which will be used later.

```
mkdir -p /some/human/data/
cd /some/human/data/
wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.aa
wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.ab
wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.ac
wget http://sourceforge.net/projects/fusioncatcher/files/data/ensembl_v84.tar.gz.ad
cat ensembl_v84.tar.gz.* | tar xz
ln -s ensembl_v84 current
```

If this works then it is not necessary to start building yourself the build data as shown below (which is **only** needed in case that the direct download for some reason does not work or one wishes to use the build data of another organism which is not available for download).

### 5.2.2 - Building yourself the organism's build data

Here, in this example, the necessary data is downloaded and necessary files/indexes are built for the **human genome** in the directory `/some/human/data/directory/` which will be used later.
```
fusioncatcher-build -g homo_sapiens -o /some/human/data/directory/
```
This takes around 5-10 hours (downloading, building indexes/databases, etc.).

In case that one wants to use a Ensembl server which is situated geographically closer, then one has:
  * Ensembl server in Europe (used by default):
```
fusioncatcher-build -g homo_sapiens -o /some/human/data/directory/
```
  * Ensembl server in East US:
```
fusioncatcher-build -g homo_sapiens -o /some/human/data/directory/ --web=useast.ensembl.org
```
  * Ensembl server in West US:
```
fusioncatcher-build -g homo_sapiens -o /some/human/data/directory/ --web=uswest.ensembl.org
```
  * Ensembl server in Asia:
```
fusioncatcher-build -g homo_sapiens -o /some/human/data/directory/ --web=asia.ensembl.org
```


In case, that it is not possible to use `fusioncatcher-build` for vary reasons (e.g. access to Ensembl website is very slow) then one may directly download the latest **human build data** (generated by `fusioncatcher-build` using Ensembl database release 84) necessary for running *FusionCatcher* from [here](http://sourceforge.net/projects/fusioncatcher/files/data/) (all files are needed and the total size is ~25 GB).

For rat genome, one has
```
fusioncatcher-build -g rattus_norvegicus -o /some/rat/data/directory/
```

For mouse genome, one has
```
fusioncatcher-build -g mus_musculus -o /some/mouse/data/directory/
```

**NOTE**: *FusionCatcher* version 0.99.6a needs a newer **build data** than the previous version (that is 0.99.5a) of 'fusioncatcher-build'.

---


# 6 - USAGE

Searching for fusion genes in a human organism, one has:
```
fusioncatcher \
-d /some/human/data/directory/ \
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/
```
where:
  * `/some/human/data/directory/` - contains the data and files generated by `fusioncatcher-build` (see [Get data](#output-data) section)
  * `/some/input/directory/containing/fastq/files/` - contains the input FASTQ (or SRA if NCBI SRA toolkit is installed) files (and not any other type of files which are not do not contain sequecing data, e.g. readme.txt)
  * `/some/output/directory/` - contains output files (for more information see [here](#output-data)):
    * `final-list_candidate_fusion_genes.txt`
    * `summary_candidate_fusions.txt`
    * `supporting-reads_gene-fusions_BOWTIE.zip`
    * `supporting-reads_gene-fusions_BLAT.zip`
    * `supporting-reads_gene-fusions_STAR.zip`
    * `viruses_bacteria_phages.txt`
    * `info.txt`
    * `fusioncatcher.log`

Searching for fusion genes in a rat organism, one has:
```
fusioncatcher \
-d /some/rat/data/directory/ \
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/
```

Searching for fusion genes in a mouse organism, one has:
```
fusioncatcher \
-d /some/mouse/data/directory/ \
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/
```



## 6.1 - Input data

The input data shall be:
  * a directory containing the input FASTQ/SRA files (this is highly recommended), or
  * several files separated by comma, e.g. `file01.fq,file02.fq` (there should be no blank before and after the comma!).

All types of raw FASTQ files produced by the Illumina Solexa and Illumina HiSeq platforms containing **paired-end** reads may be given as input. The FASTQ files shall:
  * come from an **RNA** sequencing experiment (i.e. the transcriptome/RNA is sequenced), and
  * contain **paired-end** reads, and
  * paired FASTQ files should be synchronized (i.e. reads which for a pair should be on the same line number if both FASTQ files).
  * paired-end reads which follow the suggested Illumina's sample preparation protocol, that is the two read-mates are: (i) from opposite strands, and (ii) opposite directions to one another (in other words, in order to 'bring' a read and its mate-read on the same strand then one needs to perform reverse-complement operation on only one of them)
  * paired-end reads shall come from a **stranded** (i.e. strand-specific) or **unstranded** sample preparation protocol (both are supported by *FusionCatcher*)!


It is **highly recommended** that:
  * the input FASTQ files contain the **raw** reads generated by the Illumina sequencers **without any additional trimming** (i.e. all reads from all files shall have the same length), and
  * every single input FASTQ file contains reads from only and only **one** sample/replicate (i.e. **do not concatenate** in one big FASTQ file several other FASTQ files; just give the input all FASTQ files and *FusionCatcher* will do the concatenation).

*FusionCatcher* will automatically pre-process the input reads, as follows:
  * trimming 3' end of the reads based on quality scores (default Q5),
  * removing automatically the adapter from the reads (it predicts the adapter sequence based on the reads which form a pair and also overlap and the non-overlapping parts are the predicted adapters),
  * trimming the poly A/C/G/T tails,
  * removing the reads which contain short tandem repeats (see: M. Gymrek, et al. lobSTR: A short tandem repeat profiler for personal genomes, Genome Res. 2012 Jun;22(6):1154-62, here http://genome.cshlp.org/content/22/6/1154.abstract ,
  * removing the reads which are marked as bad by Illumina sequencer,
  * removing the reads which are too short after the trimming,
  * removing the reads which map on ribosomal RNA,
  * removing the reads which map on genomes of bacteria/phages/viruses (from: [ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/](ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/) ).


The SRA files are accepted as input as long as the NCBI SRA toolkit (see dependencies section) is installed and available.

The following files are accepted/used as input:
  * `*.fq.zip`
  * `*.fq.gz`
  * `*.fq.bz2`
  * `*.fastq.zip`
  * `*.fastq.gz`
  * `*.fastq.bz2`
  * `*.txt.zip`
  * `*.txt.gz`
  * `*.txt.bz2`
  * `*.fq`
  * `*.fastq`
  * `*.sra`
and the `zip` and `gz` archives should contain only one file.

*FusionCatcher* also accepts as input also URLs (it shall start with ftp:// or http://), like for example
```
fusioncatcher -i ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030872/ERR030872_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030872/ERR030872_2.fastq.gz -o thyroid
```

The FASTQ input files should be the ones generated by the Illumina platform without any kind of additional processing/filtering/trimming. *FusionCatcher* is doing its own filtering (automatically identifies and trims the adapters, trimming if needed, quality filtering of reads, removal of rRNA reads, etc.).

*FusionCatcher* is using the file names of the input files in order to figure out which files are paired (which contains ´[R1](https://code.google.com/p/FusionCatcher/source/detail?r=1)´ reads and which contains the corresponding ´[R2](https://code.google.com/p/FusionCatcher/source/detail?r=2)´ reads). For this *FusionCatcher* is ordering alphabetically the file names (found in a input directory) and it considers that the first two form a pair, the next two are forming a pair and so on. Shortly, here the first two files, the next two files and so on should be synchronized. For example, considering the following input files:
```
L002_R1.fastq.gz
L002_R2.fastq.gz
L003_R1.fastq.gz
L003_R2.fastq.gz
```
*FusionCatcher* would automatically figure out correctly that the first two files form a pair and the next two would form another pair.

There are cases when ordering alphabetically the input file names would make *FusionCatcher* to pair wrongly the input files (i.e. the input FASTQ files are not synchronized). In this kind of cases, renaming the input files such that they fit this rule helps. For example, considering the following input files (where the file containing ´[R1](https://code.google.com/p/FusionCatcher/source/detail?r=1)´ reads is split into two files and the file containing ´[R2](https://code.google.com/p/FusionCatcher/source/detail?r=2)´ reads is split into another two files):
```
L002_R1_01.fastq.gz
L002_R1_02.fastq.gz
L002_R2_01.fastq.gz
L002_R2_02.fastq.gz
```
*FusionCatcher* would automatically figure out **wrongly** that the first two files form a pair and the next two would form another pair. In this case, renaming the files like this
```
mv L002_R1_01.fastq.gz 01_L002_R1_01.fastq.gz 
mv L002_R1_02.fastq.gz 03_L002_R1_02.fastq.gz 
mv L002_R2_01.fastq.gz 02_L002_R2_01.fastq.gz 
mv L002_R2_02.fastq.gz 04_L002_R2_02.fastq.gz 
```
would give this alphabetically order list:
```
01_L002_R1_01.fastq.gz
02_L002_R2_01.fastq.gz
03_L002_R1_02.fastq.gz
04_L002_R2_02.fastq.gz
```
and with this input files *FusionCatcher* would work correctly. Another way, around this would be to give the input files separated by comma (in the correct order and no blanks before and after the comma), like this
```
fusioncatcher -i L002_R1_01.fastq.gz,L002_R2_01.fastq.gz,L002_R1_02.fastq.gz,L002_R2_02.fastq.gz 
```

For example, this is a valid input:
```
01_L002_R1_01.fastq.gz
02_L002_R2_01.fastq.gz
03_L002_R1_02.fastq.gz
04_L002_R2_02.fastq.gz
05_L002_R1_03.fastq.gz
06_L002_R2_03.fastq.gz
07_L002_R1_04.fastq.gz
08_L002_R2_04.fastq.gz
09_L002_R1_05.fastq.gz
10_L002_R2_05.fastq.gz
11_L002_R1_06.fastq.gz
12_L002_R2_06.fastq.gz
```

For example, this is **NOT** a valid input:
```
10_L002_R2_05.fastq.gz
11_L002_R1_06.fastq.gz
12_L002_R2_06.fastq.gz
1_L002_R1_01.fastq.gz
2_L002_R2_01.fastq.gz
3_L002_R1_02.fastq.gz
4_L002_R2_02.fastq.gz
5_L002_R1_03.fastq.gz
6_L002_R2_03.fastq.gz
7_L002_R1_04.fastq.gz
8_L002_R2_04.fastq.gz
9_L002_R1_05.fastq.gz
```


**NOTE:**
  * In case that a directory is given as input, one shall make sure that the input directory does not contain files which do not contain reads sequences (e.g. readme.txt, info.txt, etc.)!
  * Please, let *FusionCatcher* _do_ the the concatenation of several FASTQ files (i.e. just put all the FASTQ files into one folder and give that folder as input to *FusionCatcher*) and do NOT do concatenate the FASTQ files yourself (e.g. using `cat`). This is because most likely different FASTQ files might have:
   * different adapter sequences (*FusionCatcher* is expecting that there are only one type of adapter, exactly like it comes directly from the Illumina sequencers),
   * different fragment sizes, and
   * different read lengths.
  * **DO NOT POOL** samples from different cell lines or from different patients! Run FusionCatcher separately with one sample at the time! It is ok the pool together the samples, which come from the (i) same cell line, or (ii) from the same patient (but still do not concatenate yourself the FASTQ files and let FusionCatcher do it for you)!

## 6.2 - Output data {#output-data}

*FusionCatcher* produces a list of candidate fusion genes using the given input data. It is recommended that this list of candidate of fusion genes is further validated in the wet-lab using for example PCR/FISH experiments.

The output files are:
  * `final-list_candidate_fusion_genes.txt` - final list with the newly found candidates fusion genes (it contains the fusion genes with their junction sequence and points); Starting with version 0.99.3c the coordinates of fusion genes are given here for human genome using **only** assembly **hg38/GRCh38**; See [Table 1](#output-data) for columns' descriptions;
  * `final-list_candidate_fusion_genes.GRCh37.txt` - final list with the newly found candidates fusion genes (it contains the fusion genes with their junction sequence and points); Starting with version 0.99.3d the coordinates of fusion genes are given here for human genome using assembly **hg19/GRCh37**; See [Table 1](#output-data) for columns' descriptions;
  * `summary_candidate_fusions.txt` - contains an executive summary (meant to be read directly by the medical doctors or biologist) of candidate fusion genes found;
  * `final-list_candidate_fusion_genes.caption.md.txt` - explains in detail the labels found in column `Fusion_description` of files `final-list_candidate_fusion_genes.txt` and `final-list_candidate_fusion_genes.GRCh37.txt`;
  * `supporting-reads_gene-fusions_BOWTIE.zip` - sequences of short reads supporting the newly found candidate fusion genes found using only and exclusively the Bowtie aligner;
  * `supporting-reads_gene-fusions_BLAT.zip` - sequences of short reads supporting the newly found candidate fusion genes found using Bowtie and Blat aligners;
  * `supporting-reads_gene-fusions_STAR.zip` - sequences of short reads supporting the newly found candidate fusion genes found using Bowtie and STAR aligners;
  * `supporting-reads_gene-fusions_BOWTIE2.zip` - sequences of short reads supporting the newly found candidate fusion genes found using Bowtie and Bowtie2 aligners;
  * `supporting-reads_gene-fusions_BWA.zip` - sequences of short reads supporting the newly found candidate fusion genes found using Bowtie and BWA aligners;`
  * `viruses_bacteria_phages.txt` - (non-zero) reads counts for each virus/bacteria/phage from NCBI database  ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/ 
  * `info.txt` - information regarding genome version, Ensembl database version, versions of tools used, read counts, etc.;
  * `fusioncatcher.log` -  log of the entire run (e.g. all commands/programs which have been run, command line arguments used, running time for each command, etc.).


FusionCatcher reports:
 * **multiple times** (up to four times) exactly the same candidate fusion gene, which has exactly the same fusion points/junction (i.e. FusionCatcher will output separately the fusions found for each of its four aligners/methods such that it is easy to see what method was used to find a fusion gene)
 * [reciprocal fusion genes](href='http://www.ncbi.nlm.nih.gov/pubmed/23340173) if they are found (e.g. geneA-geneB and also geneB-geneA)
 * every alternative splicing event found for each fusion gene (i.e. alternative fusion isoforms of the same fusion gene)


Table 1 - Columns description for file `final-list_candidate-fusion-genes.txt`

| **Column** | **Description** |
|:-----------|:----------------|
| **Gene\_1\_symbol(5end\_fusion\_partner)** | Gene symbol of the 5' end fusion partner |
| **Gene\_2\_symbol\_2(3end\_fusion\_partner)** | Gene symbol of the 3' end fusion partner |
| **Gene\_1\_id(5end\_fusion\_partner)** | Ensembl gene id of the 5' end fusion partner |
| **Gene\_2\_id(3end\_fusion\_partner)** | Ensembl gene id of the 3' end fusion partner |
| **Exon\_1\_id(5end\_fusion\_partner)** | Ensembl exon id of the 5' end fusion exon-exon junction |
| **Exon\_2\_id(3end\_fusion\_partner)** | Ensembl exon id of the 3' end fusion exon-exon junction |
| **Fusion\_point\_for\_gene\_1(5end\_fusion\_partner)** | Chromosomal position of the 5' end of fusion junction (chromosome:position:strand); 1-based coordinate |
| **Fusion\_point\_for\_gene\_2(3end\_fusion\_partner)** | Chromosomal position of the 3' end of fusion junction (chromosome:position:strand); 1-based coordinate |
| **Spanning\_pairs** | Count of pair-end reads supporting the fusion |
| **Spanning\_unique\_reads** | Count of unique reads (i.e. unique mapping positions) mapping on the fusion junction. Shortly, here are counted all the reads which map on fusion junction minus the PCR duplicated reads. |
| **Longest\_anchor\_found** | Longest anchor (hangover) found among the unique reads mapping on the fusion junction |
| **Fusion\_finding\_method** | Aligning method used for mapping the reads and finding the fusion genes. Here are two methods used which are: (i) **BOWTIE** = only Bowtie aligner is used for mapping the reads on the genome and exon-exon fusion junctions, (ii) **BOWTIE+BLAT** = Bowtie aligner is used for mapping reads on the genome and BLAT is used for mapping reads for finding the fusion junction,  (iii) **BOWTIE+STAR** = Bowtie aligner is used for mapping reads on the genome and STAR is used for mapping reads for finding the fusion junction, (iv) **BOWTIE+BOWTIE2** = Bowtie aligner is used for mapping reads on the genome and Bowtie2 is used for mapping reads for finding the fusion junction, and (v) **BOWTIE+BWA** = Bowtie aligner is used for mapping reads on the genome and Bowtie2 is used for mapping reads for finding the fusion junction. |
| **Fusion\_sequence** | The inferred fusion junction (the asterisk sign marks the junction point) |
| **Fusion\_description** | Type of the fusion gene (see the Table 2) |
| **Counts\_of\_common\_mapping\_reads** | Count of reads mapping simultaneously on both genes which form the fusion gene. This is an indication how similar are the DNA/RNA sequences of the genes forming the fusion gene (i.e. what is their homology because highly homologous genes tend to appear show as candidate fusion genes). In case of completely different sequences of the genes involved in forming a fusion gene then here it is expected to have the value zero. |
| **Predicted\_effect** | Predicted effect of the candidate fusion gene using the annotation from Ensembl database. This is shown in format **effect\_gene\_1**/**effect\_gene\_2**, where the possible values for effect\_gene\_1 or effect\_gene\_2 are: **intergenic**, **intronic**, **exonic(no-known-CDS)**, **UTR**, **CDS(not-reliable-start-or-end)**, **CDS(truncated)**, or **CDS(complete)**. In case that the fusion junction for both genes is within their CDS (coding sequence) then only the values **in-frame** or **out-of-frame** will be shown. |
| **Predicted\_fused\_transcripts** | All possible known fused transcripts in format ENSEMBL-TRANSCRIPT-1:POSITION-1/ENSEMBLE-TRANSCRIPT-B:POSITION-2, where are fused the sequence 1:POSITION-1 of transcript ENSEMBL-TRANSCRIPT-1 with sequence POSITION-2:END of transcript ENSEMBL-TRANSCRIPT-2 |
| **Predicted\_fused\_proteins** | Predicted amino acid sequences of all possible fused proteins (separated by ";").  |

Table 2 - Labels used to describe the found fusion genes (column *Fusion\_ description* from file `final-list_candidate-fusion-genes.txt`)

| **Fusion\_description** | **Description** |
|:------------------------|:----------------|
| **1000genomes**             | fusion gene has been seen in a healthy sample. It has been found in [RNA-seq data from some samples from 1000 genomes project](http://dx.doi.org/10.1371/journal.pone.0104567). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **18cancers**              | fusion gene found in a RNA-seq dataset of 18 types of cancers from 600 tumor samples published [here](http://dx.doi.org/10.1073/pnas.1606220113). |
| **adjacent**           | both genes forming the fusion are adjacent on the genome (i.e. same strand and there is no other genes situated between them on the same strand)|
| **antisense**           | one or both genes is a gene coding for [antisense RNA](http://en.wikipedia.org/wiki/Antisense_RNA)|
| **banned**              | fusion gene is on a list of known false positive fusion genes. These were found with very strong supporting data in healthy samples (i.e. it showed up in file final-list\_candidate\_fusion\_genes.txt). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **bodymap2**            | fusion gene is on a list of known false positive fusion genes. It has been found in healthy human samples collected from 16 organs from  [Illumina BodyMap2 RNA-seq database](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **cacg**                | known conjoined genes (that is fusion genes found in samples from healthy patients) from the [CACG](http://cgc.kribb.re.kr/map/) database (please see CACG database for more information). *A candidate fusion gene having this label has a very high probability of being a false positive in case that one looks for fusion genes specific to a disease.*|
| **cell\_lines**         | known fusion gene from paper: C. Klijn et al., A comprehensive transcriptional portrait of human cancer cell lines, Nature Biotechnology, Dec. 2014, [DOI:10.1038/nbt.3080](http://dx.doi.org/10.1038/nbt.3080) |
| **cgp**                 | known fusion gene from the [CGP](http://www.sanger.ac.uk/genetics/CGP/Census/) database |
| **chimerdb2**           | known fusion gene from the [ChimerDB 2](http://ercsb.ewha.ac.kr/FusionGene/) database|
| **chimerdb3kb**           | known fusion gene from the [ChimerDB 3 KB (literature curration)](http://ercsb.ewha.ac.kr/FusionGene/) database |
| **chimerdb3pub**           | known fusion gene from the [ChimerDB 3 PUB (PubMed articles)](http://ercsb.ewha.ac.kr/FusionGene/) database |
| **chimerdb3seq**           | known fusion gene from the [ChimerDB 3 SEQ (TCGA)](http://ercsb.ewha.ac.kr/FusionGene/) database |
| **conjoing**            | known conjoined genes (that is fusion genes found in samples from healthy patients) from the [ConjoinG](http://metasystems.riken.jp/conjoing/) database (please use ConjoinG database for more information regarding the fusion gene). *A candidate fusion gene having this label has a very high probability of being a false positive in case that one looks for fusion genes specific to a disease.* |
| **cosmic**              | known fusion gene from the [COSMIC](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/) database (please use COSMIC database for more information regarding the fusion gene) |
| **cta\_gene**           | one gene or both genes is CTA gene (that is that the gene name starts with **CTA-**). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **ctb\_gene**           | one gene or both genes is CTB gene (that is that the gene name starts with **CTB-**). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **ctc\_gene**           | one gene or both genes is CTC gene (that is that the gene name starts with **CTC-**). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **ctd\_gene**           | one gene or both genes is CTD gene (that is that the gene name starts with **CTD-**). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **distance1000bp**      | both genes are on the same strand and they are less than 1 000 bp apart. *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **distance100kbp**      | both genes are on the same strand and they are less than 100 000 bp apart. *A candidate fusion gene having this label has a higher probability than expected of being a false positive.* |
| **distance10kbp**       | both genes are on the same strand and they are less than 10 000 bp apart. *A candidate fusion gene having this label has a higher probability than expected of being a false positive.* |
| **duplicates**          | both genes involved in the fusion gene are [paralog](http://en.wikipedia.org/wiki/Paralog#Paralogy) for each other. For more see [Duplicated Genes Database (DGD) database](http://dgd.genouest.org/) . *A candidate fusion gene having this label has a higher probability than expected of being a false positive.* |
| **ensembl\_fully\_overlapping** | the genes forming the fusion gene are fully overlapping according to Ensembl database. *A candidate fusion gene having this label has a very high probability of being a false positive.*|
| **ensembl\_partially\_overlapping** | the genes forming the fusion gene are partially overlapping (on same strand or on different strands) according the Ensembl database. *A candidate fusion gene having this label has a good probability of being a false positive.</i> </font> |
| **ensembl\_same\_strand\_overlapping** | the genes forming the fusion gene are fully/partially overlapping and are both on the same strand according to Ensembl database. *A candidate fusion gene having this label has a very high probability of being a false positive (this is most likely and alternative splicing event).</i> </font> |
| **fragments** | the genes forming the fusion are supported by only and only one fragment of RNA. *A candidate fusion gene having this label has a medium probability of being a false positive.*|
| **gliomas**              | fusion gene found in a RNA-seq dataset of 272 glioblastoms published [here](http://dx.doi.org/10.1101/gr.165126.113). |
| **gtex**             | fusion gene has been seen in a healthy sample. It has been found in [GTEx database](http://www.gtexportal.org/home/) of healthy tissues (thru [FusionAnnotator](https://github.com/FusionAnnotator/FusionAnnotator)). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **healthy**             | fusion gene has been seen in a healthy sample. These have been found in healthy samples but the support for them is less strong (i.e. paired reads were found to map on both genes but no fusion junction was found) than in the case of **banned** label (i.e. it showed up in file preliminary list of candidate fusion genes). Also genes which have some degree of sequence similarity may show up marked like this.*A candidate fusion gene having this label has a small probability of being a false positive in case that one looks for fusion genes specific to a disease.* |
| **hpa**             | fusion gene has been seen in a healthy sample. It has been found in [RNA-seq database of 27 healthy tissues](http://dx.doi.org/10.1074/mcp.M113.035600). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **known**       | fusion gene which has been previously reported or published in scientific articles/reports/books/abstracts/databases indexed by [Google](http://www.google.com/), [Google Scholar](http://scholar.google.com/), [PubMed](http://www.ncbi.nlm.nih.gov/pubmed), etc. This label has only the role to answer with YES or NO the question "has ever before a given (candidate) fusion gene been published or reported?". This label does not have in anyway the role to provide the original references to the original scientific articles/reports/books/abstracts/databases for a given fusion gene. |
| **lincrna**             | one or both genes is a [lincRNA](http://en.wikipedia.org/wiki/Long_non-coding_RNA) |
| **matched-normal**      | candidate fusion gene (which is supported by paired reads mapping on both genes and also by reads mapping on the junction point) was found also in the matched normal sample given as input to the command line option '--normal' |
| **metazoa**             | one or both genes is a metazoa\_srp gene [Metazia\_srp](http://www.genecards.org/index.php?path=/Search/keyword/metazoa_srp) |
| **mirna**               | one or both genes is a [miRNA](http://en.wikipedia.org/wiki/MicroRNA) |
| **mt**                  | one or both genes are situated on [mitochondrion](http://en.wikipedia.org/wiki/Mitochondrion). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **non\_cancer\_tissues**   | fusion gene which has been previously reported/found in non-cancer tissues and cell lines in [Babiceanu et al, Recurrent chimeric fusion RNAs in non-cancer tissues and cells, Nucl. Acids Res. 2016](http://nar.oxfordjournals.org/content/early/2016/02/01/nar.gkw032.abstract). These are considered as non-somatic mutation and therefore they may be skipped and not reported. |
| **non\_tumor\_cells**   | fusion gene which has been previously reported/found in non-tumor cell lines, like for example HEK293. These are considered as non-somatic mutation and therefore may be skipped and not reported. |
| **no\_protein\_product** | one or both genes have no known protein product |
| **oncogene**            | one gene or both genes are a known [oncogene](http://en.wikipedia.org/wiki/Oncogene) |
| **pair\_pseudo\_genes** | one gene is the other's [pseudogene](http://en.wikipedia.org/wiki/Pseudogene). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **pancreases**           | known fusion gene found in pancreatic tumors from article: P. Bailey et al., Genomic analyses identify molecular subtypes of pancreatic cancer, Nature, Feb. 2016, http://dx.doi.org/110.1038/nature16965 |
| **paralogs**            | both genes involved in the fusion gene are  [paralog](http://en.wikipedia.org/wiki/Paralog#Paralogy) for each other (most likely this is a false positive fusion gene). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **partial-matched-normal** | candidate fusion gene (which is supported by paired reads mapping on both genes **but** _no_ reads were found which map on the junction point) was found also in the matched normal sample given as input to the command line option '--normal'. This is much weaker than **matched-normal**. |
| **prostates**           | known fusion gene found in 150 prostate tumors RNAs from paper: D. Robison et al, Integrative Clinical Genomics of Advanced Prostate Cancer, Cell, Vol. 161, May 2015, http://dx.doi.org/10.1016/j.cell.2015.05.001 |
| **pseudogene**          | one or both of the genes is a [pseudogene](http://en.wikipedia.org/wiki/Pseudogene) |
| **readthrough**         | the fusion gene is a readthrough event (that is both genes forming the fusion are on the same strand and there is no known gene situated in between); Please notice, that many of readthrough fusion genes might be false positive fusion genes due to errors in Ensembl database annotation (for example, one gene is annotated in Ensembl database as two separate genes). *A candidate fusion gene having this label has a high probability of being a false positive.* |
| **refseq\_fully\_overlapping** | the genes forming the fusion gene are fully overlapping according to RefSeq NCBI database. *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **refseq\_partially\_overlapping** | the genes forming the fusion gene are partially overlapping (on same strand or on different strands) according the RefSeq NCBI. *A candidate fusion gene having this label has a good probability of being a false positive.</i> </font> |
| **refseq\_same\_strand\_overlapping** | the genes forming the fusion gene are fully/partially overlapping and are both on the same strand according to RefSeq NCBI database. *A candidate fusion gene having this label has a very high probability of being a false positive (this is most likely and alternative splicing event).</i> </font> |
| **ribosomal\_protein**  | one or both gene is a gene encoding for [ribosomal protein](http://en.wikipedia.org/wiki/Ribosomal_protein) |
| **rp11\_gene**          | one gene or both genes is RP11 gene (that is that the gene name starts with **RP11-**). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **rp\_gene**            | one gene or both genes is RP?? gene (that is that the gene name starts with **RP??-**) where ? is a digit. *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **rrna**                | one or both genes is a [rRNA](http://en.wikipedia.org/wiki/Ribosomal_RNA).  *A candidate fusion gene having this label has a very high probability of being a false positive.*|
| **short\_distance**     | both genes are on the same strand and they are less than X bp apart, where X is set using the option '--dist-fusion' and by default it is 200 000 bp. *A candidate fusion gene having this label has a higher probability than expected of being a false positive.* |
| **similar\_reads**      | both genes have the same reads which map simultaneously on both of them (this is an indicator of how similar are the sequences of both genes; ideally this should be zero or as close to zero as possible for a real fusion). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **similar\_symbols**    | both genes have the same or very similar gene names (for example: RP11ADF.1 and RP11ADF.2). *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **snorna**              | one or both genes is a [snoRNA](http://en.wikipedia.org/wiki/Small_nucleolar_RNA) |
| **snrna**               | one or both genes is a [snRNA](http://en.wikipedia.org/wiki/Small_nuclear_RNA) |
| **tcga**                | known fusion gene from the [TCGA](https://tcga-data.nci.nih.gov/tcga/) database (please use Google for more information regarding the fusion gene) |
| **ticdb**               | known fusion gene from the [TICdb](http://www.unav.es/genetica/TICdb/) database (please use TICdb database for more information regarding the fusion gene) |
| **trna**                | one or both genes is a [tRNA](http://en.wikipedia.org/wiki/Transfer_RNA) |
| **ucsc\_fully\_overlapping** | the genes forming the fusion gene are fully overlapping according to UCSC database. *A candidate fusion gene having this label has a very high probability of being a false positive.* |
| **ucsc\_partially\_overlapping** | the genes forming the fusion gene are partially overlapping (on same strand or on different strands) according the UCSC database.  *A candidate fusion gene having this label has a good probability of being a false positive.</i> </font> |
| **ucsc\_same\_strand\_overlapping** | the genes forming the fusion gene are fully/partially overlapping and are both on the same strand according to UCSC database. *A candidate fusion gene having this label has a very high probability of being a false positive (this is most likely and alternative splicing event).</i> </font> |
| **yrna**                | one or both genes is a [Y RNA](http://en.wikipedia.org/wiki/Y_RNA) |


## 6.3 - Visualization
*FusionCatcher* outputs also the zipped FASTA files containing the reads which support the found candidate fusions genes. The files are:
  * `supporting-reads_gene-fusions_BOWTIE.zip`,
  * `supporting-reads_gene-fusions_BLAT.zip`,
  * `supporting-reads_gene-fusions_STAR.zip`,
  * `supporting-reads_gene-fusions_BOWTIE2.zip`, and
  * `supporting-reads_gene-fusions_BWA.zip`.

The reads which support the:
  * junction of the candidate fusion have their name ending with `_supports_fusion_junction`, and
  * candidate fusion (i.e. one reads map on one gene and the paired-read maps on the other fusion gene) have their name ending with `_supports_fusion_pair`.

These supporting reads (given as FASTA and FASTQ files) may be used for further visualization purposes. For example, one may use these supporting reads and align them himself/herself using his/her favourite:
  * aligner (e.g. `Bowtie/Bowtie2/TopHat/STAR/GSNAP/etc.`),
  * version/assembly of genome,
  * mapping format output (e.g. SAM/BAM), and
  * NGS visualizer (e.g. [IGV](http://www.broadinstitute.org/igv/)/[UCSC Genome Browser](http://genome.ucsc.edu/)/etc.)

### 6.3.1 - UCSC Genome Browser
For example, the sequences of supporting reads for a given candidate fusion gene may be visualized using [UCSC Genome Browser](http://genome.ucsc.edu/) by aligning them using the [UCSC Genome Browser](http://genome.ucsc.edu/)'s  BLAT aligner (i.e. copy and paste the reads here: [BLAT tool of UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgBlat?command=start) --> click the button **Submit** --> navigate into the [UCSC Genome Browser](http://genome.ucsc.edu/) to the genes that form the fusion genes). Also zooming out several times gives better view here.

### 6.3.2 - PSL format
If one uses the `--blat-visualization` command line option of the *FusionCatcher* then the BLAT alignment of the supporting reads will be done automatically by the *FusionCatcher* and the results are saved as [PSL](http://genome.ucsc.edu/FAQ/FAQformat.html#format2) files with names that are ending with `_reads.psl` in the:
  * `supporting-reads_gene-fusions_BOWTIE.zip`,
  * `supporting-reads_gene-fusions_BLAT.zip`,
  * `supporting-reads_gene-fusions_STAR.zip`, and
  * `supporting-reads_gene-fusions_BOWTIE2.zip`, and
  * `supporting-reads_gene-fusions_BWA.zip`.
The files with names ending in `_reads.psl` may be used further for visualization of the candidate fusion genes using [UCSC Genome Browser](http://genome.ucsc.edu/), [IGV (Integrative Genome Viewer)](http://www.broadinstitute.org/igv/) or any other viewer/browser which supports the [PSL](http://genome.ucsc.edu/FAQ/FAQformat.html#format2) format.

### 6.3.3 - SAM format

#### 6.3.3.1 - Automatic method
If one uses the `--visualization-sam` command line option of the *FusionCatcher* then the BOWTIE2 alignment of the supporting reads will be done automatically by the *FusionCatcher* and the results are saved as [SAM](http://samtools.github.io/hts-specs/SAMv1.pdf) files with names that are ending with `_reads.sam` in the:
  * `supporting-reads_gene-fusions_BOWTIE.zip`,
  * `supporting-reads_gene-fusions_BLAT.zip`,
  * `supporting-reads_gene-fusions_STAR.zip`,
  * `supporting-reads_gene-fusions_BOWTIE2.zip`, and
  * `supporting-reads_gene-fusions_BWA.zip`.
The files with names ending in `_reads.sam` (please note, that they still needed to be converted to BAM, coordiante sorted and indexed first) may be used further for visualization of the candidate fusion genes using [UCSC Genome Browser](http://genome.ucsc.edu/), [IGV (Integrative Genome Viewer)](http://www.broadinstitute.org/igv/) or any other viewer/browser which supports the [SAM](http://samtools.github.io/hts-specs/SAMv1.pdf) format.

#### 6.3.3.2 - Manual method
Here is an rough example of manually aligning the supporting reads (that is named as `supporting_reads.fq` in the below example; the FASTQ files needed here are the files ending in `_reads.fq` from the ZIP archives `supporting-reads_gene-fusions_*.zip` produced by *FusionCatcher*) using different aligners.
  * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) aligner (where `your_choice_of_genome_bowtie2_index` may be for human, for example [this](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip))
   * alignment done ignoring the paired-end information (i.e. like single reads):
    ```
bowtie2 \
--local \
-k 10 \
-x your_choice_of_genome_bowtie2_index \
-U supporting_reads.fq \
-S fusion_genes.sam

samtools view -bS fusion_genes.sam | samtools sort - fusion_genes.sorted

samtools index fusion_genes.sorted.bam
    ```
   * alignment done taking into account the paired-end information:
    ```
cat supporting_reads.fq | \
paste - - - - - - - - | \
awk '{print $1"\n"$2"\n"$3"\n"$4 > "r1.fq"; print $5"\n"$6"\n"$7"\n"$8 > "r2.fq"}'

bowtie2 \
--local \
-k 10 \
-x your_choice_of_genome_bowtie2_index \
-1 r1.fq \
-2 r2.fq \
-S fusion_genes.sam

samtools view -bS fusion_genes.sam | samtools sort - fusion_genes.sorted

samtools index fusion_genes.sorted.bam
    ```
  * [STAR](http://github.com/alexdobin/STAR) aligner (where `your_choice_of_genome_star_index` should be built according to the [STAR Manual](http://github.com/alexdobin/STAR/tree/master/doc))
   * alignment done ignoring the paired-end information (i.e. like single reads):
    ```
STAR \
--genomeDir your_choice_of_genome_star_index \
--alignSJoverhangMin 9 \
--chimSegmentMin 17 \
--readFilesIn supporting_reads.fq \
--outFileNamePrefix .

samtools view -bS fusion_genes.sam | samtools sort - fusion_genes.sorted

samtools index fusion_genes.sorted.bam
    ```
   * alignment done taking into account the paired-end information:
    ```
cat supporting_reads.fq | \
paste - - - - - - - - | \
awk '{print $1"\n"$2"\n"$3"\n"$4 > "r1.fq"; print $5"\n"$6"\n"$7"\n"$8 > "r2.fq"}'

STAR \
--genomeDir /your_choice_of_genome_star_index/ \
--alignSJoverhangMin 9 \
--chimSegmentMin 17 \
--readFilesIn r1.fq r2.fq\
--outFileNamePrefix .

samtools view -bS Aligned.out.sam | samtools sort - fusion_genes.sorted

samtools index fusion_genes.sorted.bam
    ```

Further, the files `fusion_genes.sorted.bam` and `fusion_genes.sorted.bam.bai` may be used with your favourite NGS visualizer!

### 6.3.4 - Chimera `R/BioConductor` package
For visualization of fusion genes found by *FusionCatcher* one may use also the `R/BioConductor` package [Chimera](http://www.bioconductor.org/packages/release/bioc/html/chimera.html), which supports *FusionCatcher*.

## 6.4 - Examples

### 6.4.1 - Example 1

Here, is an example of how *FusionCatcher* can be used to search for fusion genes in human RNA-seq sample where:
  1. any distance at chromosomal level between the candidate fusion genes is acceptable, **and**
  1. the candidate fusion genes are allowed to be readthroughs (i.e. the genes forming a fusion gene maybe adjacent on the chromosome)
  1. the candidate fusion genes are not allowed to be less the 1000 bp apart on the same strand
  1. use two methods to find the fusion genes (i.e. use BOWTIE, BLAT, STAR, and BOWTIE2 aligners for mapping the reads and this allows to find the fusion genes even in the case that the annotation from Ensembl database is not entirely correct, like for example find a fusion junction even if it is in the middle of a exon or intron)
```
fusioncatcher \
-d /some/human/data/directory/ \
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/
```

### 6.4.2 - Example 2

Here, is an example of how *FusionCatcher* can be used to search for fusion genes in human RNA-seq sample where:
  1. any distance at chromosomal level between the candidate fusion genes is acceptable, **and**
  1. the candidate fusion genes are **not** allowed to be readthroughs (i.e. there is still at least one known gene situated one the same strand in between the genes which form the candidate fusion gene)
  1. the candidate fusion genes are not allowed to be less the 1000 bp apart on the same strand
  1. use only **one** method to find the fusion genes (i.e. use only BOWTIE aligner for mapping the reads and this allows to find the fusion genes only in the case that the annotation from Ensembl database is correct, like for example find a fusion junction only if it matches perfectly the known exon borders)
```
fusioncatcher \
-d /some/human/data/directory/ \
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--skip-readthroughs \
--skip-blat
```




---


# 7 - ALIGNERS

## 7.1 - Bowtie

By default, *FusionCatcher* its the Bowtie aligner for finding candidate fusion genes. This approach relies heavily on good is the annotation data for the given organism in the Ensembl database. If, for example, a gene is not annotated well and has several exons which are not annotated in the Ensembl database and if one of these exons is the one involved in the fusion point then this fusion gene will not be found by using only the Bowtie aligner. In order to find also the fusion genes where the the junction point is in the middle of exons or introns, `*FusionCatcher*` is using by default the BLAT, and STAR aligners in addition to Bowtie aligner. The command line options '`--skip-blat`','`--skip-star`', '`--skip-bowtie2`', or '`--skip-bwa`' should be used in order to specify what aligners should not be used. The command line option '`--aligners`' specifies which aligners should be used by default. For example, '`--aligners=blat,star,bowtie2,bwa`' forces *FusionCatcher* too use all aligners for finding fusion genes

## 7.2 - Bowtie and Blat

The use of Bowtie and Blat aligners is the **default** approach of *FusionCatcher* for finding fusion genes.

In order not to use this approach the command line option '`--skip-blat`' should be added (or remove the string `blat` from line `aligners` from file `fusioncatcher/etc/configuration.cfg`), as following:

```
fusioncatcher \
-d /some/human/data/directory/ \ 
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--skip-blat
```

Please, read the license of Blat aligner before using this approach in order to see if you may use Blat! *FusionCatcher* will use Blat aligner when using this approach!

## 7.3 - Bowtie and STAR

The use of Bowtie and STAR aligners is the **default** approach of *FusionCatcher* for finding fusion genes.

In order not to use this approach the command line option '`--skip-star`' should be added, as following:

```
fusioncatcher \
-d /some/human/data/directory/ \ 
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--skip-star
```

## 7.4 - Bowtie and Bowtie2

The use of Bowtie and Bowtie2 aligners is **not** the **default** approach of *FusionCatcher* for finding fusion genes.

In order not to use this approach the command line option '`--skip-bowtie2`' should be added, as following:

```
fusioncatcher \
-d /some/human/data/directory/ \ 
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--skip-bowtie2
```

In order to use this approach the command line option '`--aligners`' should contain the string '`bowtie2`', like for example

```
fusioncatcher \
-d /some/human/data/directory/ \ 
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--aligners blat,star,bowtie2
```

## 7.5 - Bowtie and BWA

The use of Bowtie and BWA aligners is **not** the **default** approach of *FusionCatcher* for finding fusion genes.

In order not to use this approach the command line option '`--skip-bowtie2`' should be added, as following:

```
fusioncatcher \
-d /some/human/data/directory/ \ 
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--skip-bowtie2
```

In order to use this approach the command line option '`--aligners`' should contain the string '`bwa`', like for example

```
fusioncatcher \
-d /some/human/data/directory/ \ 
-i /some/input/directory/containing/fastq/files/ \
-o /some/output/directory/ \
--aligners blat,star,bwa
```

---
# 8 - Command line options

## fusioncatcher
It searchers for fusion genes and/or translocations in RNA-seq data (paired-end reads FASTQ files produced by Illumina next-generation sequencing platforms like Illumina Solexa and Illumina `HiSeq`) in diseased samples. Its command line is:
```
fusioncatcher [options]
```
and the command line options are:
```
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILENAME, --input=INPUT_FILENAME
                        The input file(s) or directory. The files should be in
                        FASTQ or SRA format and may be or not compressed using
                        gzip or zip. A list of files can be specified by given
                        the filenames separated by comma. If a directory is
                        given then it will analyze all the files found with
                        the following extensions: .sra, .fastq, .fastq.zip,
                        .fastq.gz, .fastq.bz2, fastq.xz, .fq, .fq.zip, .fq.gz,
                        .fq.bz2, fz.xz, .txt, .txt.zip, .txt.gz, .txt.bz2 .
  --batch               If this is used then batch mode is used and the input
                        specified using '--input' or '-i' is: (i) a tab-
                        separated text file containing a each line such that
                        there is one sample per line and first column are the
                        FASTQ files' full pathnames/URLs, separated by commas,
                        corresponding to the sample and an optional second
                        column containing the name for the sample, or (ii) a
                        input directory which contains a several
                        subdirectories such that each subdirectory corresponds
                        to only one sample and it contains all the FASTQ files
                        corresponding to that sample. This is useful when
                        several samples needs to be analyzed.
  --single-end          If this is used then it is assumed that all the input
                        reads are single-end reads which must be longer than
                        130 bp. Be default it is assumed that all input reads
                        come from a paired-end reads.
  -I NORMAL_MATCHED_FILENAME, --normal=NORMAL_MATCHED_FILENAME
                        The input file(s) or directory containing the healthy
                        normal-matched data. They should be given in the same
                        format as for '--input'. In case that this option is
                        used then the files/directory given to '--input' is
                        considered to be from the sample of a patient with
                        disease. This is optional.
  -o OUTPUT_DIRECTORY, --output=OUTPUT_DIRECTORY
                        The output directory where all the output files
                        containing information about the found candidate
                        fusiongenes are written. Default is 'none'.
  -d DATA_DIRECTORY, --data=DATA_DIRECTORY
                        The data directory where all the annotations files
                        from Ensembl database are placed, e.g. 'data/'. This
                        directory should be built using 'fusioncatcher-build'.
                        If it is not used then it is read from configuration
                        file specified with '--config' from 'data = ...' line.
  -T TMP_DIRECTORY, --tmp=TMP_DIRECTORY
                        The temporary directory where all the outputs files
                        and directories will be written. Default is directory
                        'tmp' in the output directory specified with '--
                        output'.
  -p PROCESSES, --threads=PROCESSES
                        Number or processes/threads to be used for running
                        SORT, Bowtie, BLAT, STAR, BOWTIE2 and other
                        tools/programs. If it is 0 (as it is by default) then
                        the number of processes/threads will be read first
                        from 'fusioncatcher/etc/configuration.cfg' file. If
                        even there it is still set to 0 then 'min(number-of-
                        CPUs-found,16)' processes will be used. Setting number
                        of threads in 'fusioncatcher/etc/configuration.cfg'
                        might be usefull in situations where one server is
                        shared between several users and in order to limit
                        FusionCatcher using all the CPUs/resources.Default is
                        '0'.
  --config=CONFIGURATION_FILENAME
                        Configuration file containing the paths to external
                        tools (e.g. Bowtie, Blat, fastq-dump.) in case that
                        they are not specified in PATH variable! Default is '/
                        apps/fusioncatcher/etc/configuration.cfg,/apps/fusionc
                        atcher/bin/configuration.cfg'.
  -z, --skip-update-check
                        Skips the automatic routine that contacts the
                        FusionCatcher server to check for a more recent
                        version. Default is 'False'.
  -V, --keep-viruses-alignments
                        If it is set then the SAM alignments files of reads
                        mapping on viruses genomes are saved in the output
                        directory for later inspection by the user. Default is
                        'False'.
  -U, --keep-unmapped-reads
                        If it is set then the FASTQ files, containing the
                        unmapped reads (i.e. reads which do not map on genome
                        and transcriptome), are saved in the output directory
                        for later inspection by the user. Default is 'False'.
  --aligners=ALIGNERS   The aligners to be used on Bowtie aligner. By default
                        always BOWTIE aligner is used and it cannot be
                        disabled. The choices are:
                        ['blat','star','bowtie2','bwa']. Any combination of
                        these is accepted if the aligners' names are comma
                        separated. For example, if one wants to used all four
                        aligners then 'blat,star,bowtie2,bwa' should be given.
                        The command line options '--skip-blat', '--skip-star',
                        and '--skip-bowtie2' have priority over this option.
                        If the first element in the list is the configuration
                        file (that is '.cfg' file) of FusionCatcher then the
                        aligners specified in the list of aligners specified
                        in the configuration file will be used (and the rest
                        of aligner specified here will be ignored). In case
                        that the configuration file is not found then the
                        following aligners from the list will be used. Default
                        is
                        '/apps/fusioncatcher/etc/configuration.cfg,blat,star'.
  --skip-blat           If it is set then the pipeline will NOT use the BLAT
                        aligner and all options and methods which make use of
                        BLAT will be disabled. BLAT aligner is used by
                        default. Please, note that BLAT license does not allow
                        BLAT to be used for commercial activities. Fore more
                        information regarding BLAT please see its license:
                        <http://users.soe.ucsc.edu/~kent/src/>. Default is
                        'False'.
  --skip-star           If it is set then the pipeline will NOT use the STAR
                        aligner and all options and methods which make use of
                        STAR will be disabled. STAR aligner is used by
                        default. Default is 'False'.
  --sort-buffer-size=SORT_BUFFER_SIZE
                        It specifies the buffer size for command SORT. Default
                        is '80%' if less than 32GB installed RAM else is set 
                        to 26 GB.
  --start=START_STEP    It re-starts executing the workflow/pipeline from the
                        given step number. This can be used when the pipeline
                        has crashed/stopped and one wants to re-run it from
                        from the step where it stopped without re-running from
                        the beginning the entire pipeline. 0 is for restarting
                        automatically and 1 is the first step. Default is '0'.


```

## fusioncatcher-build
It downloads the necessary data for a given organism from the Ensembl database and it builds the necessary files/indexes which are needed to running *FusionCatcher*. Its command line is:
```
fusioncatcher-build [options]
```
and the command line options are:
```
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -o OUTPUT_DIRECTORY, --output=OUTPUT_DIRECTORY
                        The output directory where all the outputs files  and
                        directories will be written.
  -c CONFIGURATION_FILENAME, --config=CONFIGURATION_FILENAME
                        Configuration file containing the paths to external
                        tools (e.g. Bowtie, etc.) in case that they are not in
                        PATH! Default is '/apps/fusioncatcher/bin/../etc/confi
                        guration.cfg,/apps/fusioncatcher/bin/configuration.cfg
                        '.
  -g ORGANISM, --organism=ORGANISM
                        Organism for which the data is downloaded from Ensembl
                        database and built, for example: 'homo_sapiens',
                        'mus_musculus', 'rattus_norvegicus',
                        'canis_familiaris', etc. Default is 'homo_sapiens'.
  -w WEB_ENSEMBL, --web=WEB_ENSEMBL
                        Ensembl database web site from where the data is
                        downloaded.  e.g. 'www.ensembl.org',
                        'uswest.ensembl.org', 'useast.ensembl.org',
                        'asia.ensembl.org', etc. Default is 'www.ensembl.org'.
  -e FTP_ENSEMBL, --ftp-ensembl=FTP_ENSEMBL
                        Ensembl database FTP site from where the data is
                        downloaded. Default is 'ftp.ensembl.org'.
  --ftp-ensembl-path=FTP_ENSEMBL_PATH
                        The path for Ensembl database FTP site from where the
                        data is downloaded.
  -x FTP_UCSC, --ftp-ucsc=FTP_UCSC
                        UCSC database FTP site from where the data is
                        downloaded. Default is 'hgdownload.cse.ucsc.edu'.
  -n FTP_NCBI, --ftp-ncbi=FTP_NCBI
                        NCBI database FTP site from where the data is
                        downloaded. Default is 'ftp.ncbi.nlm.nih.gov'.
  --skip-blat           If it is set then the pipeline will NOT use the BLAT
                        aligner and all options and methods which make use of
                        BLAT will be disabled. BLAT aligner is used by
                        default. Please, note that BLAT license does not allow
                        BLAT to be used for commercial activities. Fore more
                        information regarding BLAT please see its license:
                        <http://users.soe.ucsc.edu/~kent/src/>. Default is
                        'False'.
  --enlarge-genes       If it is set then the genes are enlarged (i.e. their
                        introns include also in the transcriptome). Default is
                        'False'.
  -p PROCESSES, --threads=PROCESSES
                        Number or processes/threads to be used. Default is
                        '0'.
  --skip-database=SKIP_DATABASE
                        If it is set then the pipeline will skip the specified
                        database(s). The choices are ['cosmic','conjoing','chi
                        merdb2','ticdb','cgp','cacg']. If several databases
                        should be skipped, then their names shall be separated
                        by comma. Default is ''.
  -s START_STEP, --start=START_STEP
                        It starts executing the workflow from the given step
                        number. This can be used when the pipeline has
                        crashed/stopped and one wants to re-run it from from
                        the step where it stopped without re-running from the
                        beginning the entire pipeline. 0 is for restarting
                        automatically and 1 is the first step. This is
                        intended to be used for debugging. Default is '0'.
  -l HASH, --hash=HASH  Hash to be used for computing checksum. The choices
                        are ['no','crc32','md5','adler32','sha512','sha256'].
                        If it is set up to 'no' then no checksum is used and
                        the entire pipeline is executed as a normal shell
                        script. For more information see 'hash_library' in
                        'workflow.py'. This is intended to be used for
                        debugging. Default is 'no'.
  -k, --keep            Preserve intermediate files produced during the run.
                        By default, they are NOT deleted upon exit. This is
                        intended to be used for debugging. Default value is
                        'False'.
  -u CHECKSUMS_FILENAME, --checksums=CHECKSUMS_FILENAME
                        The name of the checksums file. This is intended to be
                        used for debugging. Default value is 'checksums.txt'.

```


---

# 9 - Methods

The main goal of *FusionCatcher* is to find **somatic** (and/or pathogenic) fusion genes in RNA-seq data.

*FusionCatcher* is doing its own quality filtering/trimming of reads. This is needed because most a very important factor for finding fusion genes in RNA-seq experiment is the length of RNA fragments.  **Ideally** the RNA fragment size for finding fusion genes should be over 300 bp.  Most of the RNA-seq experiments are designed for doing differentially expression analyses and not for finding fusion genes and therefore the RNA fragment size many times is less than 300bp and the trimming and quality filtering should be done in such a way that it does not decrease even more the RNA fragment size.

*FusionCatcher* is able to finding fusion genes even in cases where the fusion junction is within known exon or within known intron (for example in the middle of an intron) due to the use of BLAT aligner. The minimum condition for *FusionCatcher* to find a fusion gene is that both genes involved in the fusion are annotated in Ensembl database (even if their gene structure is "wrong").

*FusionCatcher* is spending most of computational analysis on the most promising fusion genes candidate and tries as early as possible to filter out the candidate fusion genes which do not look promising, like for example:
  * candidate fusion gene is composed of a gene and its pseudogene, or
  * candidate fusion gene is composed of a gene and its paralog gene, or
  * candidate fusion gene is composed of a gene and a miRNA gene (but a gene which contains miRNA genes are not skipped), or
  * candidate fusion gene is composed of two genes which have a very sequence similarity (i.e. *FusionCatcher* is computing its homology score), or
  * candidate fusion gene is known to be found in samples from healthy persons (using the 16 organs RNA-seq data from the Illumina BodyMap2), or
  * candidate fusion gene is in one of the known databases of fusion genes found in healthy persons, i.e. ChimerDB2, CACG, and ConjoinG.

*FusionCatcher* is using by default three aligners for mapping the reads. The aligners are Bowtie, BLAT, and STAR. STAR is used here only and only for "splitting" the reads while aligning them.


---

# 10 - Comparisons to other tools
When performing comparisons where *FusionCatcher* is compared with other gene fusions finder we **always recommend strongly to use the default/recommended parameters** for *FusionCatcher* and also to use the raw FASTQ files which came directly from the Illumina sequencer.


The performance of *FusionCatcher* is decreased drastically, when using other parameters than the default/recommended ones! Especially **do not change** the defaults for: `--5keep`, `--anchor-fusion`, `--reads-fusion`, `--pairs-fusion`, `--pairs-fusion2`! The default parameters should work just fine for input reads which have the size range between 35 bp to 250 bp.


Also, when comparing the fusion genes found by *FusionCatcher* with fusion genes found by other tools one needs to keep in mind that *FusionCatcher* is a **SOMATIC** fusion gene finder and **NOT** a (general) fusion gene finder. This means that if a fusion gene is already known to exist in healthy individuals (from public literature or from our internal RNA-seq database of healthy sample) then that fusion gene will be skipped by *FusionCatcher* and it will not be reported at all! An example is the well known fusion gene TTTY15-USP9Y which is known to be found in healthy individuals (see [here](http://www.sciencedirect.com/science/article/pii/S0002944015001996)) and which *FusionCatcher* will skip it and will not report it on purpose because **it is not a somatic fusion gene**!

Also, when one is running *FusionCatcher* on some synthetic/simulated RNA-seq datasets which contain a set of random/ad-hoc fusion genes which are created randomly and without any biological support (for example, that fusion gene has never been reported in the literature to exist in a diseased patient), there most likely *FusionCatcher* will detect that these **random/ad-hoc** fusion genes are not fitting the already known biological knowledge (e.g. ad-hoc/random fusion gene might have been reported already to exist in healthy patients, or ad-hoc/random fusion is between a gene its paralog/homolog/pseudogene) and will skip them and will not report them even if it finds them. Therefore we strongly recommend not to run *FusionCatcher* on synthetic/simulated RNA-seq dataset which are known to contain fusion genes which are **not** somatic fusion genes. Also, we strongly recommend not to run *FusionCatcher* on downsampled input datasets, like for example, choosing randomly 30 million reads from an original datasets with 60 million reads. *FusionCatcher* has been specifically built for analyzing real input RNA-seq datasets which come directly from the sequencing machine.



---

# 11 - License
*FusionCatcher*'s code is released under [GNU GPL version 3 license](http://www.gnu.org/copyleft/gpl.html). *FusionCatcher* is using third-party tools and databases. The user is responsible to obtain licenses for the third-party tools and databases which are used by *FusionCatcher*.

**Most** (but not all) of the third-party tools and databases used by *FusionCatcher* are (i) free to use, or (ii) are released under GPL/MIT-type licenses. The most notable exception here of which we are aware is BLAT's aligner license, which requires one to buy a license when BLAT is used in commercial environment (please, see for more [here](http://www.kentinformatics.com/contact-us.html)). In case that one does not wish to use BLAT aligner then it is still possible to use *FusionCatcher* for finding fusion genes, by telling *FusionCatcher* not to use BLAT aligner but instead to use the BOWTIE2 aligner (BLAT is used by default and BOWTIE2 is not used by default), as following:

```
/apps/fusioncatcher/bin/fusioncatcher \
--aligners star,bowtie2
```



---

# 12 - Citing
If you use *FusionCatcher*, please cite:

D. Nicorici, M. Satalan, H. Edgren, S. Kangaspeska, A. Murumagi, O. Kallioniemi, S. Virtanen, O. Kilkku, **FusionCatcher – a tool for finding somatic fusion genes in paired-end RNA-sequencing data**, bioRxiv, Nov. 2014, [DOI:10.1101/011650](http://biorxiv.org/content/early/2014/11/19/011650)


---

# 13 - Reporting Bugs

Please, when reporting bugs include also the following files:
  * "fusioncatcher.log" (this contains just a list of the commands executed by *FusionCatcher*), and
  * "info.txt" (this contains info: regarding the version of *FusionCatcher* and tools used, statistics about input FASTQ files, counts of found reads, etc.)
which were generated by *FusionCatcher* during the run.

**NOTE**: Giving only step number where the error has appeared is not enough because the step numbers depend on the input type (e.g. raw FASTQ file, ZIP compressed FASTQ file, SRA file, etc.) and command line options used to run *FusionCatcher* (e.g. some command line option skip some steps). The step numbers are used by *FusionCatcher* for being able to re-start from the last step which was executed successfully last time in case that last time the run ended prematurely due to reasons which didn't depend on *FusionCatcher* (e.g. server crashed).


---

# NOTES
  * <font color='red'>The performance of <b>FusionCatcher</b> is decreased drastically, when using other parameters than the default/recommended ones! Especially <b>do not change</b> the defaults for: <code>--5keep, --anchor-fusion, --reads-fusion, --pairs-fusion, --pairs-fusion2</code>! The default parameters should work just fine for input reads which have the size range between 35 bp to 250 bp.</font>
  * * <font color='red'>The performance of <b>FusionCatcher</b> is decreased drastically, when running <b>FusionCatcher</b> on a subset of the reads. It is not recommended to run <b>FusionCatcher</b> on 20 million paired-reads sampled from a sample. All the reads from the sample shall be given as input to <b>FusionCatcher</b></font>
  * `fusioncatcher-build` takes several hours to run and it depends on the local internet connection speed. It needs to be run only once!
  * *FusionCatcher* can be run many times using the same data produced by the `fusioncatcher-build`;
  * Ensembl version 84 was found to work fine with *FusionCatcher* as May 2016;
  * *FusionCatcher* and `fusioncatcher-build` restart automatically from the point where have been interrupted at the previous run.
  * *FusionCatcher* by default is focusing on finding fusion genes specific to diseased/tumor/cancer samples. That means that *FusionCatcher* will skip the fusion genes which are already known to exist in healthy samples. If one wishes to find fusion genes in healthy samples then we suggest other fusion finders to be used.
  * *FusionCatcher* is able to find fusion genes **also** without using BLAT aligner but in this case we recommend to user BOWTIE2 aligner (which is not used by default) also in order to compensate!



