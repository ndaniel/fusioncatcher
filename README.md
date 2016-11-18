
FusionCatcher
=============

Finder of somatic fusion-genes in RNA-seq data.


Download / Install / Update / Upgrade [FusionCatcher](http://github.com/ndaniel/fusioncatcher)
----------------------------------------------------------------------------------------------

Use this one-line command:

```bash
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download
```

For more installing options, see:

```bash
bootstrap.py --help
```

Description
-----------
FusionCatcher searches for novel/known somatic fusion genes, translocations, and
chimeras in RNA-seq data (paired-end or single-end reads from Illumina NGS platforms 
like Solexa/HiSeq/NextSeq/MiSeq/MiniSeq) from diseased/cancer samples.

The aims of FusionCatcher are:
 * very good detection rate for finding candidate somatic fusion
   genes (see somatic mutations; using a matched normal sample is
   optional; several databases of known fusion genes found in healthy
   samples are used as a list of known false positives; biological
   knowledge is used, like for example gene fusion between a gene and
   its pseudogene is filtered out),
 * very good RT-PCR validation rate of found candidate somatic fusion
   genes (this is very important for us),
 * very easy to use (i.e. no a priori knowledge of bioinformatic
   databases and bioinformatics is needed in order to run FusionCatcher BUT
   Linux/Unix knowledge is needed; it allows a very high level of control
   for expert users),
 * to be as automatic as possible (i.e. choose automatically the best
   parameters for the given input data, e.g. finding automatically the 
   adapters, quality trimming of reads, building the exon-exon junctions 
   automatically based on the length of the reads given as input, etc. 
   whilst giving full control to expert users) whilst providing the best 
   possible detection rate for finding somatic fusion genes (with very 
   low rate of false positives and very good sensitivity).


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
Old releases and the latest release of FusionCatcher are on https://sourceforge.net/projects/fusioncatcher/files/

Citing
------
D. Nicorici, M. Satalan, H. Edgren, S. Kangaspeska, A. Murumagi, O. Kallioniemi,
S. Virtanen, O. Kilkku, FusionCatcher – a tool for finding somatic fusion genes
in paired-end RNA-sequencing data, bioRxiv, Nov. 2014, 
[DOI:10.1101/011650](http://dx.doi.org/10.1101/011650)

