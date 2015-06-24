#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==============================================================================
FusionCatcher
==============================================================================
FusionCatcher searches for novel somatic fusion genes in RNA-seq paired-end
reads data produced by the Illumina Solexa platforms (for example: Illumina
HiSeq 2000, GAIIx, GAII).

Required dependencies:
- Linux 64-bit (e.g. Ubuntu 12.04/14.04 or newer)
- Python version 2.7.x (>=2.6.0 and < 3.0 is fine)
- BioPython version 1.65 (>=1.50 is fine)
- Bowtie 64-bit version 1.1.1  <http://bowtie-bio.sourceforge.net/index.shtml>
- Bowtie2 64-bit version 2.2.5  <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>
- STAR version 2.4.1c <http://github.com/alexdobin/STAR>
- BWA version 0.7.12  <http://sourceforge.net/projects/bio-bwa/>
- SEQTK version 1.0-r68e-dirty  <http://github.com/ndaniel/seqtk>

Optional dependencies:
- strongly recommended (expected by default to be installed but their use can
  be disabled by using the command line option '--skip-blat'):
   * BLAT version 0.35 <http://users.soe.ucsc.edu/~kent/src/>. Executables are
     available at <http://hgdownload.cse.ucsc.edu/admin/exe/>. Please, check the
     license to see if it allows you to run/use it!
   * faToTwoBit <http://users.soe.ucsc.edu/~kent/src/>. Executables are
     available at <http://hgdownload.cse.ucsc.edu/admin/exe/>. Please, check the
     license to see if it allows you to run/use it!

- nice to have:
   * Velvet (de novo assembler) version 1.2.10
     <http://www.ebi.ac.uk/~zerbino/velvet/>
   * fastq-dump version 2.3.2 from NCBI SRA Toolkit
     <http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>
   * Python library openpyxl <http://pypi.python.org/pypi/openpyxl>
   * Python library xlrd <http://pypi.python.org/pypi/xlrd>
   * SAMTools <http://www.htslib.org/>

Genomic Databases (used automatically by FusionCatcher):
- ENSEMBL database <http://www.ensembl.org>
- COSMIC database <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>
- TICdb database <http://www.unav.es/genetica/TICdb/>
- IMGT/HLA Database <http://www.ebi.ac.uk/ipd/imgt/hla/>
- ChimerDB 2.0 database (literature-based annotation)
  <http://ercsb.ewha.ac.kr/FusionGene/index.jsp>
- ConjoinG database <http://metasystems.riken.jp/conjoing/>
- Cancer Genome Project (CGP) translocations database
  <http://www.sanger.ac.uk/genetics/CGP/Census/>
- CACG database <http://cgc.kribb.re.kr/map/>
- UCSC database <http://hgdownload.cse.ucsc.edu/downloads.html#human>
- RefSeq database (thru UCSC database)
- Viruses/bacteria/phages NCBI genomes database <ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/>
- DGD database <http://dgd.genouest.org/>


Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

This file is part of FusionCatcher.

FusionCatcher is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FusionCatcher is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FusionCatcher (see file 'COPYING.txt').  If not, see
<http://www.gnu.org/licenses/>.

By default, FusionCatcher is running BLAT aligner
<http://users.soe.ucsc.edu/~kent/src/> but it offers also the option to disable
all its scripts which make use of BLAT aligner if you choose explicitly to do so.
BLAT's license does not allow to be used for commercial activities. If BLAT
license does not allow to be used in your case then you may still use
FusionCatcher by forcing not use the BLAT aligner by specifying the option
'--skip-blat'. Fore more information regarding BLAT please see its license.

Please, note that FusionCatcher does not require BLAT in order to find
candidate fusion genes!
"""

import sys
import os
import struct
import optparse
import multiprocessing
import subprocess
import shutil
import socket
import locale
import math
import configuration


# for sort in linux
locale.setlocale(locale.LC_ALL, 'C')

def expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))


def islink(alink = None):
    """
    Wrapper for: os.path.islink()
    """
    f = False
    if alink:
        alink = alink[:-1] if alink.endswith(os.sep) else alink
        if os.path.islink(alink):
            f = True
    return f

# get the path of this script
pipeline_path = os.path.dirname(expand(sys.argv[0]))

def outdir(*more_paths):
    global out_dir
    return os.path.join(out_dir,*more_paths)

def datadir(*more_paths):
    global data_dir
    return os.path.join(data_dir,*more_paths)

def tmpdir(*more_paths):
    global tmp_dir
    return os.path.join(tmp_dir,*more_paths)

# make sure that a directory ends with path separator such that workflow can
# recognize it as directory
def adir(a_dir):
    if (not a_dir.endswith('\\')) and (not a_dir.endswith('/')):
        a_dir = a_dir + os.sep
    return a_dir

#
#
# test if a command line option has been passed
def is_optparse_provided(parser, dest):
    r = False
    for opt in parser._get_all_options():
        if opt.dest == dest:
            if opt._long_opts and opt._long_opts[0] in sys.argv[1:]:
                r = True
                break
            if opt._short_opts and opt._short_opts[0] in sys.argv[1:]:
                r = True
                break
    return r
#    if any (opt.dest == dest and (opt._long_opts[0] in sys.argv[1:] or (False if (not opt._short_opts) else opt._short_opts[0] in sys.argv[1:])) for opt in parser._get_all_options()):
#        return True
#    return False


#
#
#
def empty(a_file):
    f = True
    if (os.path.isfile(a_file) or islink(a_file)):
        s = os.path.getsize(a_file)
        if s < 1000:
            d = [line for line in file(a_file,'r').readlines() if line.rstrip('\r\n')]
            if d:
                f = False
        else:
            f = False
    return f


#
#
#
def delete_file(some_file):
    some_file = some_file[:-1] if some_file.endswith(os.sep) else some_file
    if os.path.isfile(some_file) or islink(some_file):
        os.remove(some_file)

#
#
#
def info(ajob, fromfile, tofile , top = "\n\n\n", bottom = "\n\n\n" , temp_path = 'no'):
    if ajob.run():
        aux = open(tofile,'a')
        top = str(top).splitlines() if type(top).__name__ == 'str' else top
        bottom = str(bottom).splitlines() if type(bottom).__name__ == 'str' else bottom
        ajob.write("APPENDING to file: '%s'.\n"% (tofile,))
        for line in top:
            t = line.rstrip('\r\n')+'\n'
            aux.write(t)
            ajob.write(">%s" % (t,))
        if fromfile:
            ajob.write(">from file: '%s'\n"% (fromfile,))
            for line in file(fromfile,'r').readlines():
                t = line.rstrip('\r\n')+'\n'
                aux.write(t)
        for line in bottom:
            t = line.rstrip('\r\n')+'\n'
            aux.write(t)
            ajob.write(">%s" % (t,))
        aux.close()

    if fromfile:
        ajob.clean(fromfile,temp_path = temp_path)

#
# command line parsing
#

class MyOptionParser(optparse.OptionParser):
    def format_epilog(self, formatter):
        return self.epilog


def is_known_extension(something):
    kx = ['fastq.gz','.fq.gz',
          '.fastq.bz2','.fq.bz2',
          '.fastq.zip','.fq.zip',
          '.fastq.xz','.fq.xz',
          '.fastq', '.fq',
          '.sra',
          '.bam','.sam']
    sign = False
    for ekx in kx:
        if something.lower().endswith(ekx):
            sign = True
            break
    # skip readme files
    if (something.lower().startswith('readme') or
        something.lower().startswith('index.') or
        something.lower().startswith('checksum') or
        something.startswith('.') or
        something.lower().startswith('md5')):
        sign = False
    return sign


usage = "%prog [options]"

epilog = ("\n" +
         "Author: Daniel Nicorici \n" +
         "Email: Daniel.Nicorici@gmail.com \n" +
         "Copyright (c) 2009-2015, Daniel Nicorici \n " +
         "\n")

description = ("FusionCatcher searches for novel gene fusions in RNA-seq \n"+
               "paired-end reads data produced by the Illumina sequencing \n"+
               "platforms (like for example: Illumina HiSeq 2500, \n"+
               "Illumina HiSeq 2000, Illumina HiSeq X, Illumina NextSeq 500, \n"+
               "Illumina GAIIx, Illumina GAII, Illumina MiSeq). \n")

version = "%prog 0.99.4c beta"


if __name__ == "__main__":

    parser = MyOptionParser(
                usage       = usage,
                epilog      = epilog,
                description = description,
                version     = version
             )

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = "The input file(s) or directory. The files should be "+
                             "in FASTQ or SRA format and may be or not compressed "+
                             "using gzip or zip. "+
                             "A list of files can be specified by given the "+
                             "filenames separated by comma. If a directory is given "+
                             "then it will analyze all the files found with the "+
                             "following extensions: .sra, "+
                             ".fastq, .fastq.zip, .fastq.gz, .fastq.bz2, fastq.xz, "+
                             ".fq, .fq.zip, .fq.gz, .fq.bz2, fz.xz, "+
                             ".txt, .txt.zip, .txt.gz, .txt.bz2 ."
                      )

    parser.add_option("--batch",
                      action = "store_true",
                      dest = "batch_mode",
                      default = False,
                      help = "If this is used then batch mode is used "+
                             "and the input specified using '--input' or '-i' is: "+
                             "(i) a tab-separated text file containing a each line such "+
                             "that there is one sample per line and first column are the "+
                             "FASTQ files' full pathnames/URLs, separated by commas, corresponding to the "+
                             "sample and an optional second column containing the name for the sample, or "+
                             "(ii) a input directory which contains a several subdirectories such that each "+
                             "subdirectory corresponds to only one sample and it contains all the FASTQ files "+
                             "corresponding to that sample. This is useful when several samples needs to be analyzed."
                      )

    parser.add_option("--single-end",
                      action = "store_true",
                      dest = "single_end",
                      default = False,
                      help = "If this is used then it is assumed that all the input reads are single-end reads "+
                             "which must be longer than 130 bp. "+
                             "Be default it is assumed that all input reads come from a paired-end reads."
                      )
                      


    parser.add_option("--normal","-I",
                      action = "store",
                      type = "string",
                      dest = "normal_matched_filename",
                      help = "The input file(s) or directory containing the "+
                             "healthy normal-matched data. They should be given in the same "+
                             "format as for '--input'. In case that this option is used "+
                             "then the files/directory given to '--input' is considered "+
                             "to be from the sample of a patient with disease. This is optional."
                      )

    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_directory",
                      help = "The output directory where all the output files "+
                             "containing information about the found candidate fusion"+
                             "genes are written. Default is '%default'.")

    parser.add_option("--data","-d",
                      action = "store",
                      type = "string",
                      dest = "data_directory",
                      help = "The data directory where all the annotations files "+
                             "from Ensembl database are placed, e.g. 'data/'. "+
                             "This directory should be built using 'fusioncatcher-build'. "+
                             "If it is not used then it is read from configuration file "+
                             "specified with '--config' from 'data = ...' line.")

    parser.add_option("--tmp","-T",
                      action = "store",
                      type = "string",
                      dest = "tmp_directory",
                      default = "tmp",
                      help = "The temporary directory where all the outputs files "+
                             "and directories will be written. Default is directory "+
                             "'%default' in the output directory specified with '--output'. ")

    parser.add_option("--threads","-p",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = multiprocessing.cpu_count(),
                      help = "Number or processes to be used for running Bowtie, "+
                             "BLAT, STAR, BOWTIE2 and other tools/programs. "+
                             "By default the number of processes is chosen to "+
                             "be equal to the number of found CPUs. "+
                             "Default is '%default'. ")

    parser.add_option("--config",
                      action = "store",
                      type = "string",
                      dest = "configuration_filename",
                      default = os.path.abspath(os.path.join(pipeline_path,"..","etc","configuration.cfg"))+','+os.path.abspath(os.path.join(pipeline_path,"configuration.cfg")),
                      help = "Configuration file containing the paths to external "+
                             "tools (e.g. Bowtie, Blat, fastq-dump.) in case that "+
                             "they are not specified in PATH variable! "+
                             "Default is '%default'.")

    parser.add_option("--skip-update-check","-z",
                      action = "store_true",
                      dest = "skip_update_check",
                      default = False,
                      help = "Skips the automatic routine that contacts the "+
                             "FusionCatcher server to check for a more recent version. "+
                             "Default is '%default'. ")


    group = optparse.OptionGroup(parser, "Trimming options")

    group.add_option("--5keep","-l",
                      action = "store",
                      type = "int",
                      dest = "trim_3end_keep",
                      default = 60, # 60 # 68
                      help = "This may be seen as seed length. For Bowtie aligner the reads "+
                             "longer than '%default' will be trimmed from "+
                             "their 3-end such that to become exactly '%default' bp long. "+
                             "Reads shorter than '%default' will not be trimmed. " +
                             "The trimming priorities are '--5end','--3end','--5keep'. " +
                             "if several trimming options are used simultaneously. "+
                             "The trimming is done by default only to the reads used "+
                             "for BOWTIE aligner but not for BLAT/STAR/BOWTIE2/BWA aligners. In order "+
                             "to apply the trimming also the reads used by BLAT/STAR/BOWTIE2/BWA aligners "+
                             "option '--trim-psl' should be used! The trimming of reads for "+
                             "BLAT/STAR/BOWTIE2/BWA aligners is done using the option '--trim-psl-5keep'. "+
                             "Default is '%default'.")

    group.add_option("--5end","-5",
                      action = "store",
                      type = "int",
                      dest = "trim_5end",
                      default = 0,
                      help = "It trims all the reads from their 5-end with the "+
                             "given size. "+
                             "The trimming priorities are '--5end','--3end','--5keep'. " +
                             "if several trimming options are used simultaneously. "+
                             "The trimming is done by default only to the reads used "+
                             "for BOWTIE aligner but not for BLAT aligner. In order "+
                             "to apply the trimming also the reads used by BLAT/STAR/BOWTIE2/BWA aligners "+
                             "option '--trim-psl' or '--trim-psl-5end' should be used! "+
                             "Default is '%default'.")

    group.add_option("--3end","-3",
                      action = "store",
                      type = "int",
                      dest = "trim_3end",
                      default = 0,
                      help = "It trims all the reads from their 3-end with the "+
                             "given size. "+
                             "The trimming priorities are '--5end','--3end','--5keep'. " +
                             "if several trimming options are used simultaneously. "+
                             "The trimming is done by default only to the reads used "+
                             "for BOWTIE aligner but not for BLAT aligner. In order "+
                             "to apply the trimming also the reads used by BLAT/STAR/BOWTIE2/BWA aligners "+
                             "option '--trim-psl' should be used! "+
                             "Default is '%default'.")

    group.add_option("--trim-psl",
                      action = "store_true",
                      dest = "trim_psl",
                      default = False,
                      help = "If it is specified then also the reads given as input "+
                             "to BLAT/STAR/BOWTIE2/BWA aligners are trimmed using the parameters given "+
                             "by command line arguments '--5keep', '--5end', and '--3end'. "+
                             "By default the trimming options "+
                             "'--5keep', '--5end', '--3end' are trimming the reads only for "+
                             "for the BOWTIE method but not when BLAT/STAR/BOWTIE2/BWA are used. "+
                             "Default is '%default'.")


    group.add_option("--trim-psl-5keep","-x",
                      action = "store",
                      type = "int",
                      dest = "trim_psl_3end_keep",
                      default = 82, # 80
                      help = "This may be seen as seed length. All reads given as input "+
                             "to BLAT/STAR/BOWTIE2/BWA aligners and which "+
                             "are longer than '%default' will be trimmed from "+
                             "their 3-end such that to become exactly '%default' bp long. "+
                             "The reads given as input to Bowtie are not trimmed using this "+
                             "option. It should be set to 0 if no trimming should be done "+
                             "for BLAT/STAR/BOWTIE2/BWA. "+
                             "Default is '%default'.")

    group.add_option("--trim-psl-5end",
                      action = "store_true",
                      dest = "trim_psl_5end",
                      default = False,
                      help = "If it is specified then also the reads given as input "+
                             "to BLAT/STAR/BOWTIE2/BWA aligners are trimmed using the parameters given "+
                             "by command line argument '--5end'. "+
                             "By default the trimming options "+
                             "'--5keep', '--5end', '--3end' are trimming the reads only for "+
                             "for the BOWTIE method but not when BLAT/STAR/BOWTIE2/BWA are used. "+
                             "Default is '%default'.")

    group.add_option("--trim-quality","-Q",
                      action = "store",
                      dest = "trim_quality",
                      type = "int",
                      default = 5,
                      help = "The input reads will be trimmed from their 3'end "+
                             "when the quality scores are below the given threshold, e.g. 5 for Q5. "+
                             "Default is '%default'.")

    group.add_option("--trim-wiggle",
                      action = "store",
                      dest = "trim_wiggle",
                      type = "int",
                      default = 0, # it was 2
                      help = "The input reads will be trimmed during the alignment "+
                             "from their 5' and 3' ends for filtering only purposes. "+
                             "Default is '%default'.")

    group.add_option("--trimfq",
                      action = "store",
                      dest = "trimfq",
                      type = "float",
                      default = 1.00,
                      help = "If this is set less than 1.00 the quality then the quality "+
                             "trimming will be done using Phred algorithm in addition to "+
                             "quality filtering which is already done by default. "+
                             "For this the 'seqtk trimfq' tool is used and also the input "+
                             "reads should have quality score in Sanger format. A recommended value "+
                             "here for quality trimming is 0.05 (which is the default value of 'seqtk trimfq') or 0.10."
                    )


    parser.add_option_group(group)


    group = optparse.OptionGroup(parser, "Search fusion genes options")

    group.add_option("--mismatches","-m",
                      action = "store",
                      type = "int",
                      dest = "mismatches",
                      default = 2,
                      help = "Maximum number of mismatches to be allowed for "+
                             "mapping reads using Bowtie aligner. "+
                             "Minimum accepted value is zero and maximum is 3. "+
                             "Default is '%default'.")

    group.add_option("--mismatches-psl",
                      action = "store",
                      type = "int",
                      dest = "mismatches_psl",
                      default = 2,
                      help = "Maximum number of mismatches to be allowed for "+
                             "mapping reads using BLAT/STAR/BOWTIE2/BWA aligner. "+
                             "Default is '%default'.")

    group.add_option("--mismatches-ambiguous",
                      action = "store",
                      type = "int",
                      dest = "ambiguous_mismatches",
                      default = 2,
                      help = "Maximum number of mapping mismatches for which the "
                             "same reads are considered mapping ambiguously. "+
                             "Default is '%default'.")

    group.add_option("--mismatches-filtering",
                      action = "store",
                      type = "int",
                      dest = "filter_mismatches",
                      default = 2,
                      help = "Maximum number of mapping mismatches used for filtering the reads. "+
                             "Default is '%default'.")

    group.add_option("--pairs-fusion","-s",
                      action = "store",
                      dest = "spanning_pairs",
                      default = "3,3,3,3,3",
                      help = "The minimum number of paired-end reads which "+
                             "support a candidate fusion gene and which will be "+
                             "considered for further analysis.  "+
                             "It is given separated by commas for each of "+
                             "the aligners: BOWTIE, BLAT, STAR, BOWTIE2, BWA (in this order). " +
                             "Default is '%default'.")

    group.add_option("--top-pairs-fusion",
                      action = "store",
                      type = "int",
                      dest = "spanning_pairs_count",
                      default = 10000,
                      help = "If the '--pairs-fusion' selects more than N preliminary "+
                             "candidate fusion genes then only the first N will be "+
                             "considered for further analysis. N is set here. "+
                             "Default is '%default'.")

    group.add_option("--reads-fusion","-r",
                      action = "store",
                      dest = "spanning_reads",
                      default = "2,2,2,2,2",
                      help = "The minimum number of reads which "+
                             "support a candidate fusion gene that is the minimum "+
                             "number of reads which overlap over the fusion "+
                             "junction. It is given separated by commas for each of "+
                             "the aligners: BOWTIE, BLAT, STAR, BOWTIE2, BWA (in this order). " +
                             "Default is '%default'.")


    group.add_option("--anchor-fusion","-a",
                      action = "store",
                      dest = "length_anchor",
                      default = "17,17,17,17,17", # default 14; 17; 18?
                      help = "The minimum length which a read should overlap over (or anchor/overhang for) "+
                             "fusion junction of a candidate fusion gene in order to be considered for " +
                             "further analysis. Minimum accepted value is 10 and it should not exceed half "+
                             "of the length of the longest read from the input data. It is given separated "+
                             "by commas for each of the aligners: BOWTIE, BLAT, STAR, BOWTIE2, BWA (in this order). " +
                             "Default is '%default'.")

    group.add_option("--anchor-fusion2","-W",
                      action = "store",
                      type = "int",
                      dest = "length_anchor2",
                      default = 47, # default 22?
                      help = "If the anchor/overhang which supports the fusion is longer (or equal) than "+
                             "this value than the required number of reads supporting the fusion is 1. " +
                             "It basically overrides '--reads-fusion*' for anchors longer (or equal) than"+
                             "the value specified here. It always should be larger than the value "+
                             "specified by '--reads-fusion*'. "+
                             "Default is '%default'.")

    mydefault = sorted([
            "paralogs",
            "pair_pseudo_genes",
            "similar_reads",
            "ambiguous",
            "similar_symbols",
            "ensembl_fully_overlapping",
            "ensembl_same_strand_overlapping",
#            'ucsc_fully_overlapping',
#            'ucsc_same_strand_overlapping',
            'refseq_fully_overlapping',
            'refseq_same_strand_overlapping',
            "distance1000bp",
            "rrna",
            "trna",
            "mt",
            "mirna",
            "yrna",
            "7skrna",
            "snorna",
            "snrna",
            "rp11_gene",
            "cta_gene",
            "ctb_gene",
            "ctc_gene",
            "ctd_gene",
            "rp_gene",
            "banned",
            "healthy",
            "conjoing",
            "metazoa",
            "bodymap2",
            "non-tumor_cell-lines",
            "removed"])
    all_choices = sorted([
            'paralogs',
            'adjacent',
            'ambiguous',
            'distance1000bp',
            'chimerdb2',
            'cacg',
            'duplicates',
            'bodymap2',
            'metazoa',
            'similar_reads',
            'similar_symbols',
            'short_distance',
            'yrna',
            '7skrna',
            'rrna',
            'trna',
            'mt',
            'lincrna',
            'mirna',
            'pseudogene',
            'snorna',
            'snrna',
            'pair_pseudo_genes',
            'rp_gene' ,
            'rp11_gene',
            'ensembl_fully_overlapping',
            'ensembl_partially_overlapping',
            'ensembl_same_strand_overlapping',
            'ribosomal_protein',
            'cta_gene',
            'ctb_gene',
            'ctc_gene',
            'ctd_gene',
            'conjoing',
            'healthy',
            'ucsc_fully_overlapping',
            'ucsc_partially_overlapping',
            'ucsc_same_strand_overlapping',
            'refseq_fully_overlapping',
            'refseq_partially_overlapping',
            'refseq_same_strand_overlapping',
            'gencode_fully_overlapping',
            'gencode_partially_overlapping',
            'gencode_same_strand_overlapping',
            'distance10kbp',
            'distance100kbp',
            'banned',
            'non-tumor_cell-lines'
            'removed'])
    group.add_option("--filter-fusion","-b",
                      action = "store",
                      type = "string",
                      dest = "biotypes",
                      default = ','.join(sorted(mydefault)),
                      help = "Candidate gene fusions to be skipped from further "+
                             "analysis in case that one of "+
                             "partner gene or both genes (which form a fusion) "+
                             "are specified here. "+
                             "All possible values are: ["+', '.join(sorted(all_choices))+"]. "+
                             "'short_distance' is used for labeling the "+
                             "candidate fusion genes which do meet the criteria "+
                             "specified with '--min-dist-fusion'. "+
                             "Several can be chosen but in this case they " +
                             "should comma separated. "+
                             "Default is '%default'.")

    group.add_option("--filter-fusion-add","-B",
                      action = "store",
                      type = "string",
                      dest = "biotypes_more",
                      help = "Any label of fusion genes specified here will be "+
                             "appended to the list given to '--filter-fusion'. "+
                             "This is just an easy way to add more to '--filter-fusion'. "+
                             "For more read the description of '--filter-fusion'. "+
                             "Default is '%default'.")

    group.add_option("--dist-fusion","-D",
                      action = "store",
                      type = "int",
                      dest = "min_dist",
                      default = 200000,
                      help = "The candidate fusion genes where the distance "+
                             "between the genes is below this threshold will be marked "+
                             "using the label 'custom_distance' "+
                             "Default is '%default'.")

    group.add_option("--all-reads-fusion","-A",
                      action = "store_true",
                      dest = "all_reads_junction",
                      default = False,
                      help = "If it is specified then all reads (reads which form "+
                             "a pair and single reads which do not have a mate "+
                             "read because their mate has been removed due to "+
                             "different reasons, like for example low quality), "+
                             "will be used for finding the fusion point, which "+
                             "is the exon-exon junction. If not specified then only "+
                             "reads which form a pair will be used for "+
                             "finding the exon-exon junction (one read maps on one "+
                             "of the transcripts of the gene involved in the fusion "+
                             "and its mate will map on the exon-exon junction). "+
                             "Default is '%default'.")

    group.add_option("--homolog-fusion","-H",
                      action = "store",
                      type = "float",
                      dest = "homolog",
                      default =  float(1)/float(8*(10**4)),#float(1)/float(2*(10**5)), # float(1)/float(8*(10**4)),# float(1)/float(5*(10**4)),
                      help = "The minimum number of reads (as percentage [0..1]) "+
                             "which map simultaneously "+
                             "onto two genes in order to be considered homologous. "+
                             "If set to 0 then no homology analysis is done. "+
                             "This information is used for filtering out candidate "+
                             "fusion genes which are homologous. "+
                             "Default is '%default'.")

    group.add_option("--visualization-psl",
                      action = "store_true",
                      dest = "psl_visualization",
                      default = False,
                      help = "If it is set then the pipeline will use the BLAT "+
                             "aligner for aligning the reads which support the "+
                             "newly found candidate fusion genes.  Please, note "+
                             "that BLAT license does not allow BLAT to be used for "+
                             "commercial activities.  Fore more information "+
                             "regarding BLAT please see its license: "+
                             "<http://users.soe.ucsc.edu/~kent/src/>.  Also please, note "+
                             "that this option is not actively supported anymore and "+
                             "in the future will be deprecated.  If one still wants "+
                             "to use it, then one should run this 'faToTwoBit genome.fa genome.2bit -noMask') "+
                             "in 'fusioncatcher/data/current/'.  Instead it is recommended to use "+
                             "'--visualization-sam'. "+
                             "Default is '%default'.")

    group.add_option("--visualization-sam",
                      action = "store_true",
                      dest = "sam_visualization",
                      default = False,
                      help = "If it is set then the pipeline will use the BOWTIE2 "+
                             "aligner for aligning the reads which support the "+
                             "newly found candidate fusion genes. "+
                             "Default is '%default'.")

    group.add_option("--assembly","-M",
                     action = "store_true",
                     dest = "assembly",
                     default = False,
                     help = "If used then the reads found to support the newly "+
                            "found candidate fusion genes are assembled using "+
                            "VELVET <http://www.ebi.ac.uk/~zerbino/velvet/>. "+
                            "Default is '%default'.")


    parser.add_option_group(group)


    group = optparse.OptionGroup(parser, "Bioinformatics sonication (for input reads longer than 130bp only)")


    group.add_option("--sonication",
                      action = "store",
                      type = "int",
                      dest = "sonication",
                      default =  130,
                      help = "In case that the input reads are longer than the threshold set here "+
                             "then they will be broken up bioinformatically in smaller reads. "+
                             "If this is set to 0 then no break up will be done. "+
                             "Default is '%default'.")

    group.add_option("--bridges",
                      action = "store",
                      type = "int",
                      dest = "bridges",
                      default = 0,
                      help = "Number of encompasses paired-reads to be generated for each input long read. "+
                             "If it is set to 0 then the number will chosen automatically based on "+
                             "the length of input reads, i.e. ceil(length_read/160). " +
                             "Default is '%default'.")

    parser.add_option_group(group)



    group = optparse.OptionGroup(parser, "Reads filtering/processing options")

    group.add_option("--skip-filter-mt",
                      action = "store_true",
                      dest = "skip_mitochondrion_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads "+
                             "which map on the mitochondrion. "+
                             "Default is '%default'.")

    group.add_option("--skip-filter-phix174",
                      action = "store_true",
                      dest = "skip_phix174_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads "+
                             "which map on the 'Enterobacteria phage phiX174' genome. "+
                             "Default is '%default'.")

    group.add_option("--skip-filter-hla",
                      action = "store_true",
                      dest = "skip_hla_filtering",
                      default = False,
                      help = "If specified then it filters out the reads "+
                             "which map on the HLA (histocompatibility) genes. "+
                             "Default is '%default'.")

    group.add_option("--skip-filter-vir",
                      action = "store_true",
                      dest = "skip_viruses_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads "+
                             "which map on known genomes of viruses. "+
                             "Default is '%default'.")

    group.add_option("--skip-filter-str",
                      action = "store_true",
                      dest = "skip_str_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads "+
                             "which contain STR (short tandem repeats). "+
                             "Default is '%default'.")

    group.add_option("--skip-filter-b",
                      action = "store_true",
                      dest = "skip_b_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads with "+
                             "B quality scores (i.e. low quality) which are a special "+
                             "indicator in "+
                             "Fastq Illumina files. Default is '%default'.")

    group.add_option("--skip-filter-ambiguous",
                      action = "store_true",
                      dest = "skip_ambiguous_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads which "+
                             "maps ambiguously (i.e. same read map simultaneously on two "+
                             "locuses on genome/transcriptome within 0-3 mismatches. "
                             "Default is '%default'.")

    group.add_option("--skip-filter-adapter",
                      action = "store_true",
                      dest = "skip_adapter_filtering",
                      default = False,
                      help = "If specified then it skips filtering out the reads which "+
                             "contains the adapters. "+
                             "Default is '%default'.")

    group.add_option("--skip-filter-psl",
                      action = "store_true",
                      dest = "skip_prefiltering_psl",
                      default = False,
                      help = "If it is set then the pipeline will not prefilter "+
                             "the short reads which will be used for doing BLAT/STAR/BOWTIE2/BWA alignment. "+
                             "By default, the short reads are prefiltered before "+
                             "being aligned using BLAT/STAR/BOWTIE2 in order to speed up the BLAT/STAR/BOWTIE2/BWA "+
                             "alignment which is time and computationally demanding. "+
                             "The disadvantage of doing prefiltering is that the sensitivity "+
                             "of BLAT/STAR/BOWTIE2/BWA alignment is somewhat lowered. "+
                             "Default is '%default'.")

    group.add_option("--skip-interleave",
                      action = "store_true",
                      dest = "skip_interleave_processing",
                      default = False,
                      help = "If specified then it skips interleaving the short reads "+
                             "from the input FASTQ files. The program tries automatically "+
                             "to pair the forward and reverse short reads based on file "+
                             "names. In case that the pair is done wronlgy then this "+
                             "argument should be used to remedy the problem. "+
                             "Default is '%default'.")

    group.add_option("--skip-known-fusions",
                      action = "store_true",
                      dest = "skip_known_fusions",
                      default = False,
                      help = "If it is set then the pipeline will not use its own database "+
                             "and COSMIC database of already known fusion genes! "+
                             "Here skipping means that the known fusion genes will "+
                             "treated as any other candidate fusion genes "+
                             "and if there is enough evidence will be shown in the "+
                             "final list of found fusion genes. By default, the known "+
                             "fusion genes are treated preferentially and are pushed "+
                             "directly to the very final step of finding the junction "+
                             "point. " +
                             "Default is '%default'.")

    group.add_option("--skip-adjacent",
                      action = "store_true",
                      dest = "skip_adjacent",
                      default = False,
                      help = "If it is set then the pipeline will not seach for "+
                             "candidate fusion genes where the genes are adjacent! "+
                             "By default the candidate fusion genes which have "+
                             "genes that are adjacent are analysed also but in many cases "+
                             "they are just annotation errors in the Ensembl database "+
                             "and maybe they are not real fusion genes. "+
                             "Default is '%default'.")

    group.add_option("--skip-banned-fusions",
                      action = "store_true",
                      dest = "skip_banned_fusions",
                      default = False,
                      help = "If it is set then the list of known banned fusion "+
                             "genes (which are found in healthy samples) is not used. "+
                             "Default is '%default'.")

    group.add_option("--keep-viruses-alignments","-V",
                      action = "store_true",
                      dest = "keep_viruses",
                      default = False,
                      help = "If it is set then the SAM alignments files of reads mapping on "+
                             "viruses genomes are saved in the output directory "+
                             "for later inspection by the user. "+
                             "Default is '%default'.")

    group.add_option("--keep-unmapped-reads","-U",
                      action = "store_true",
                      dest = "keep_unmapped_reads",
                      default = False,
                      help = "If it is set then the FASTQ files, containing "+
                             "the unmapped reads (i.e. reads which do not map "+
                             "on genome and transcriptome), are saved in the output directory "+
                             "for later inspection by the user. "+
                             "Default is '%default'.")

    group.add_option("--skip-compress-ids",
                      action = "store_true",
                      dest = "skip_compress_ids",
                      default = False,
                      help = "If it is set then the reads ids will not be compressed "+
                             "(i.d. renamed) using lossy compression and "+
                             "the original reads ids will be kept thru the whole "+
                             "run of FusionCatcher. Be default the reads ids will be "+
                             "compressed using lossy compression. "+
                             "Default is '%default'.")

    group.add_option("--skip-automatic-scaling",
                      action = "store_true",
                      dest = "skip_automatic_scaling",
                      default = False,
                      help = "If it is set then the thresholds for anchor length, "+
                             "spanning reads, and spanning pairs will not be adjusted "+
                             "automatically according to the input reads. "+
                             "Default is '%default'.")

    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, "Usage of aligners")


    group.add_option("--aligners",
                      action = "store",
                      type = "string",
                      dest = "aligners",
                      default = "blat,star",
                      help = "The aligners to be used on Bowtie aligner. "+
                             "By default always BOWTIE aligner is used and it "+
                             "cannot be disabled. The choices are: "+
                             "['blat','star','bowtie2','bwa']. Any combination of these is "+
                             "accepted if the aligners' names are comma separated. "+
                             "For example, if one wants to used all four aligners "+
                             "then 'blat,star,bowtie2,bwa' should be given. "+
                             "The command line options '--skip-blat', '--skip-star', "+
                             "and '--skip-bowtie2' have priority over this option. "
                             "Default is '%default'.")

    group.add_option("--skip-blat",
                      action = "store_true",
                      dest = "skip_blat",
                      default = False,
                      help = "If it is set then the pipeline will NOT use the BLAT "+
                             "aligner and all options and methods which make use of "+
                             "BLAT will be disabled. "+
                             "BLAT aligner is used by default. Please, note "+
                             "that BLAT license does not allow BLAT to be used for "+
                             "commercial activities. Fore more information "+
                             "regarding BLAT please see its license: "+
                             "<http://users.soe.ucsc.edu/~kent/src/>. "+
                             "Default is '%default'.")

    group.add_option("--skip-star",
                      action = "store_true",
                      dest = "skip_star",
                      default = False,
                      help = "If it is set then the pipeline will NOT use the STAR "+
                             "aligner and all options and methods which make use of "+
                             "STAR will be disabled. "+
                             "STAR aligner is used by default. " +
                             "Default is '%default'.")

    group.add_option("--skip-star-bowtie",
                      action = "store_true",
                      dest = "skip_star_bowtie",
                      default = False,
                      help = "If it is set then the pipeline will NOT use the BOWTIE "+
                             "aligner within the the usage of STAR aligner. By default it is tried to rescue "+
                             "STAR's partially mapped reads by runing "+
                             "again BOWTIE aligner and stich together the the partialy mapped reads. " +
                             "Default is '%default'.")

    group.add_option("--skip-bowtie2",
                      action = "store_true",
                      dest = "skip_bowtie2",
                      default = False,
                      help = "If it is set then the pipeline will NOT use the BOWTIE2 "+
                             "aligner and all options and methods which make use of "+
                             "BOWTIE2 will be disabled. "+
                             "BOWTIE2 aligner is used by default. " +
                             "Default is '%default'.")

    group.add_option("--skip-bwa",
                      action = "store_true",
                      dest = "skip_bwa",
                      default = False,
                      help = "If it is set then the pipeline will NOT use the BWA "+
                             "aligner and all options and methods which make use of "+
                             "BWA will be disabled. "+
                             "BWA aligner is used by default. " +
                             "Default is '%default'.")


    group.add_option("--skip-conversion-grch37",
                      action = "store_true",
                      dest = "skip_conversion_grch37",
                      default = False,
                      help = "If it is set then the fusion coordinates for human "+
                             "genome version GRCh38 (which is default) will NOT be "+
                             "reported also using version GRCh37. " +
                             "Default is '%default'.")

    group.add_option("--limit-blat",
                      action = "store",
                      type = "int",
                      dest = "limit_blat",
                      default = 3 * (2**30),
                      help = "The maximum limit of the genome's size which BLAT aligner "+
                             "is able to handle.  If the genome is larger than this limit "+
                             "then it will be split automatically in smaller pieces such that "+
                             "the aligner can handle them without an error. "+
                             "Default is '%default'.")

    group.add_option("--limit-bowtie",
                      action = "store",
                      type = "int",
                      dest = "limit_bowtie",
                      default = 2**32 - 100000,
                      help = "The maximum limit of the genome's size which BOWTIE aligner "+
                             "is able to handle.  If the genome is larger than this limit "+
                             "then it will be split automatically in smaller pieces such that "+
                             "the aligner can handle them without an error. "+
                             "Default is '%default'.")

    group.add_option("--limit-bowtie2",
                      action = "store",
                      type = "int",
                      dest = "limit_bowtie2",
                      default = 30*(10**6), # 300*(10**6)
                      help = "The maximum limit of the genome's size which BOWTIE2 aligner "+
                             "is able to handle.  If the genome is larger than this limit "+
                             "then it will be split automatically in smaller pieces such that "+
                             "the aligner can handle them without an error. "+
                             "Default is '%default'.")

    group.add_option("--limit-bwa",
                      action = "store",
                      type = "int",
                      dest = "limit_bwa",
                      default = 1,
                      help = "The maximum limit of the genome's size which MWA aligner "+
                             "is able to handle.  If the genome is larger than this limit "+
                             "then it will be split automatically in smaller pieces such that "+
                             "the aligner can handle them without an error. "+
                             "Default is '%default'.")

    group.add_option("--limit-star",
                      action = "store",
                      type = "int",
                      dest = "limit_star",
                      default = 500*(10**6),#1,#30*(10**6), # 3 * (2**30)
                      help = "The maximum limit of the genome's size which STAR aligner "+
                             "is able to handle.  If the genome is larger than this limit "+
                             "then it will be split automatically in smaller pieces such that "+
                             "the aligner can handle them without an error. "+
                             "Default is '%default'.")





    parser.add_option_group(group)


    group = optparse.OptionGroup(parser, "Search sensitivity of fusions (they are mutually exclusive)")

    group.add_option("--sensitive",
                      action = "store_true",
                      dest = "sensitive",
                      default = False,
                      help = "This will increase the sensitivity of finding fusion genes "+
                             "at the expense of finding slightly more false positives. "+
                             "Default is '%default'.")

    group.add_option("--mildly-sensitive",
                      action = "store_true",
                      dest = "mildly_sensitive",
                      default = False,
                      help = "This will increase the sensitivity of finding fusion genes "+
                             "at the expense of finding slightly more false positives. "+
                             "Default is '%default'.")

    group.add_option("--highly-sensitive",
                      action = "store_true",
                      dest = "highly_sensitive",
                      default = False,
                      help = "This will increase the sensitivity of finding fusion genes "+
                             "at the expense of finding slightly more false positives. "+
                             "Default is '%default'.")

    group.add_option("--paranoid-sensitive",
                      action = "store_true",
                      dest = "paranoid_sensitive",
                      default = False,
                      help = "This will increase the sensitivity of finding fusion genes "+
                             "at maximum at the expense of very high rate of false positives. "+
                             "It is not recommended to be used for finding novel fusion genes. "+
                             "It may be used in cases when one searches for already known fusion "+
                             "genes which were missed in previous runs with default sensitivity. "+
                             "Default is '%default'.")

    parser.add_option_group(group)



    group = optparse.OptionGroup(parser, "Custom labels for fusions and (pre)filtering of fusions")

    group.add_option("--label-title",
                      action = "store",
                      type = "string",
                      dest = "label_title",
                      help = "The label title(s) which will be used to mark the "+
                             "candidate fusion genes given as input to '--label-file'. "+
                             "If several are given then they should be comma separated. "+
                             "If '--label-title' is used then also '--label-file' should be specified.")

    group.add_option("--label-file",
                      action = "store",
                      type = "string",
                      dest = "label_file",
                      help = "File(s) containing pairs of Ensembl gene ids (separated by tab, "+
                             "i.e. first two columns) "+
                             "candidate fusion genes given as input to '--label-file'. "+
                             "If several are given then they should be comma separated. "+
                             "If '--label-file' is used then also '--label-title' should be specified.")

    group.add_option("--label-threshold",
                      action = "store",
                      type = "string",
                      dest = "label_threshold",
                      help = "The thresholds which might be given as an optional column 3 "+
                             "in '--label-file'. All the pairs given in '--label-file' "+
                             "which have the number of column number 3 strictly less than "+
                             "this given threshold will be ignored/skipped. "+
                             "If several are given then they should be comma separated.")

    group.add_option("--focus",
                      action = "store",
                      type = "string",
                      dest = "focus_fusions",
                      help = "It contains a tab separated file text containd two columns "+
                             "with Ensembl gene ids for candidate fusion genes which will "+
                             "be forced to pass the preliminary filtering. This lines should "+
                             "be sorted and also the columns should be sorted.")

    group.add_option("--reads-preliminary-fusions",
                      action = "store_true",
                      dest = "reads_preliminary_fusions",
                      default = False,
                      help = "The sequences of all reads which support the preliminary candidate fusion genes "+
                             "are extracted. "+
                             "Default is '%default'.")

    parser.add_option_group(group)

    #
    # debugging arguments
    #

    group = optparse.OptionGroup(parser, "SeqTK subseq tool")

    group.add_option("--extract-buffer-size",
                      action = "store",
                      type = "int",
                      dest = "extract_buffer_size",
                      default = 1*(10**9), # default = 2*(10**9)
                      help = "The size of memory buffer used by Python script (for extracting reads from a FASTQ file based on a list of reads ids). "+
                             "This depends more on the amount of memory which Python environment is able to handle and less than the free actual free RAM memory on the computer where this is run. "+
                             "It might be that the default value is too high and needs to be lowered, e.g. '500000000' "+
                             "be required to be lowered. This is meant to be used together with '--split-seqtk-subseq 0'. " +
                             "Default is '%default'.")

    group.add_option("--split-seqtk-subseq",
                      action = "store",
                      type = "int",
                      dest = "split_seqtk_subseq",
                      default = 1,
                      help = "The input file (i.e. file containing read ids) of 'SEQTK SUBSEQ' will "+
                             "be splitted in a number of parts specified here. If it is 1 then no spliting is done. "+
                             "If it is set to 0 then 'SEQTK SUBSEQ' will not be used and instead an alternative Python script is used. "+
                             "Setting this to 0 or 2 or larger values is meant to be used in cases when "+
                             "'SEQTK SUBSEQ' fails due to not enough memory. "+
                             "Default is '%default'.")

    parser.add_option_group(group)

    #
    # debugging arguments
    #

    group = optparse.OptionGroup(parser, "Debug Options")

    group.add_option("--start",
                      action = "store",
                      type = "int",
                      dest = "start_step",
                      default = 0,
                      help = "It starts executing the workflow from the given "+
                             "step number. This can be used when the pipeline "+
                             "has crashed/stopped and one wants to re-run it from "+
                             "from the step where it stopped without re-running " +
                             "from the beginning the entire pipeline. "+
                             "0 is for restarting automatically and 1 is the first step. "+
                             "Default is '%default'.")

    choices = ('no','crc32','md5','adler32','sha512','sha256')
    group.add_option("--hash",
                      action = "store",
                      type = "choice",
                      choices = choices,
                      dest = "hash",
                      default = "no",
                      help = "Hash to be used for computing checksum. The choices "+
                             "are ['"+"','".join(choices)+"']. "+
                             "If it is set up to 'no' then no checksum is used and "+
                             "the entire pipeline is executed as a normal shell "+
                             "script. For more information see 'hash_library' in "+
                             "'workflow.py'. "+
                             "Default is '%default'.")

    group.add_option("--keep",
                      action = "store_true",
                      dest = "keep_temporary_files",
                      default = False,
                      help = "Preserve intermediate files produced during the run. "+
                             "By default, they are deleted upon exit. "+
                             "Default value is '%default'.")

    group.add_option("--checksums",
                      action = "store",
                      type = "string",
                      dest = "checksums_filename",
                      default = 'checksums.txt',
                      help = "The name of the checksums file. "+
                             "Default value is '%default'. ")

    group.add_option("--bowtie-chunkmbs",
                      action = "store",
                      type = "int",
                      dest = "chunkmbs",
                      default = 128, # 128
                      help = "The value to be passed to the '--chunkmbs' command line option of Bowtie aligner. "+
                             "Default is '%default'.")

    group.add_option("--long",
                      action = "store_true",
                      dest = "long_report",
                      default = False,
                      help = "A slightly longer report for fusion genes will be generated (i.e. fusions genes will be given per each aligner used). "+
                             "Default value is '%default'.")


    parser.add_option_group(group)

    ################################################################################
    #
    # MAIN
    #
    ################################################################################

    #command line parsing
    (options, args) = parser.parse_args()

    #
    # validate options
    #
    if (
        (not options.input_filename) or
        (not options.output_directory)
        ):
        parser.print_help()
        print "EXAMPLE:"
        print "========"
        print ""
        print "fusioncatcher  -d /some/data/directory/  -i /some/input/directory/containing/fastq/files/  -o /some/output/directory/"
        print ""
        print ""
        print "where /some/data/directory/ contains data which was built previously by running:"
        print ""
        print ""
        print "fusioncatcher-build -g homo_sapiens -o /some/data/directory/"
        print ""
        print "or it has been downloaded (see: FusionCatcher's manual, section 'Downloading/building organism's data')."
        print ""
        print "NOTE:"
        print "'fusioncatcher-build' needs to be run only once (for each organism"
        print "or when the Ensembl database is updated) and 'fusioncatcher'"
        print "will reuse the '/some/data/directory/'."
        print ""
        print ""
        print >>sys.stderr, "ERROR: input/output directory is not specified!"
        print >>sys.stderr, ""
        sys.exit(1)

    print "Checking Python version..."
    pythonversion = sys.version_info
    if pythonversion >= (2,6) and pythonversion < (3,0) and struct.calcsize("P") * 8 >= 64:
        print "  * Compatible Python version found!"
    else:
        print >>sys.stderr, "  * ERROR: Found Python version: %s.%s !\n" % (pythonversion[0],pythonversion[1])
        print >>sys.stderr, "           The Python should be 64-bit and the version should be >=2.6.0 and < 3.0 !"
        sys.exit(1)

    print "Checking size of installed RAM memory ..."
    total_memory = 0
    try:
        total_memory = int(os.popen("free -m").readlines()[1].split()[1])
    except:
        pass
    if total_memory != 0:
        if total_memory < 23000:
            print >>sys.stderr, "  * ERROR: %d MB of RAM memory found (minimum of 24 GB of RAM memory is needed)!" % (total_memory,)
            sys.exit(1)
        else:
            print "  * %d MB of RAM memory found!" % (total_memory,)
    else:
        print >>sys.stderr, "  * Warning: Not able to detect the size of installed RAM memory!"

    if options.trim_3end > 0 and options.trim_3end_5end_keep > 0:
        parser.error("ERROR: Arguments '--5keep'and '--3end' are mutually exclusive!")
        sys.exit(1)

    if options.tmp_directory == "tmp":
        options.tmp_directory = os.path.join(expand(options.output_directory),"tmp")


    if expand(options.input_filename) == expand(options.output_directory):
        parser.error("ERROR: Input and output paths should be different!")
        sys.exit(1)

    if expand(options.input_filename) == expand(options.tmp_directory):
        parser.error("ERROR: Input and temporary paths should be different!")
        sys.exit(1)

    if expand(options.output_directory) == expand(options.tmp_directory):
        parser.error("ERROR: Output and temporary paths should be different!")
        sys.exit(1)

    x1 = expand(options.output_directory)
    x2 = options.output_directory
    if  x1.find(",") != -1 or x2.find(",") != -1:
        parser.error("ERROR: Output path contains comma(s)!")
        sys.exit(1)

    x1 = expand(options.tmp_directory)
    x2 = options.tmp_directory
    if  x1.find(",") != -1 or x2.find(",") != -1:
        parser.error("ERROR: Temporary path contains comma(s)!")
        sys.exit(1)



    if options.batch_mode:
        cc = sys.argv[:]
        if cc[0].endswith('fusioncatcher.py'):
            cc[0].replace('fusioncatcher.py','fusioncatcher-batch.py')
        cc = [e for e in cc if e and e!='--batch']
        r = os.system(' '.join(cc))
        if r:
            print >>sys.stderr,"Error while running fusioncatcher-batch.py!"
            sys.exit(1)
        sys.exit(0)

    if options.normal_matched_filename:
        cc = sys.argv[:]
        tempo = adir(expand(options.tmp_directory))
        if not os.path.exists(tempo):
            os.makedirs(tempo)
        tempo_input = os.path.join(tempo,'fusioncatcher-input.log')
        tempo_normal = os.path.join(tempo,'fusioncatcher-normal.log')
        file(tempo_input,'w').write(options.input_filename)
        file(tempo_normal,'w').write(options.normal_matched_filename)
        if cc[0].endswith('fusioncatcher.py'):
            cc[0] = cc[0][:-3]+'-batch.py'
        elif cc[0].endswith('fusioncatcher'):
            cc[0] = cc[0]+'-batch.py'
        dd = []
        for ik in cc:
            if ik.find("=") == -1:
                dd.append(ik)
            else:
                gi = ik.split("=")
                dd.append(gi[0])
                dd.append(gi[1])
        next_i = False
        next_n = False
        com = []
        for ik in dd:
            if ik == '-i' or ik == '--input':
                next_i = True
                continue
            if next_i:
                next_i = False
                continue
            if ik == '--normal' or ik == '-I':
                next_n = True
                continue
            if next_n:
                next_n = False
                continue
            com.append(ik)

        dn = com + ['--input',tempo_input,'--normal',tempo_normal]
        dn = ' '.join(dn)
        print "--------------------------"
        print dn
        print "--------------------------"
        r = os.system(dn)
        if r:
            print >>sys.stderr,"Error while running fusioncatcher-batch.py!"
            sys.exit(1)
        os.remove(tempo_input)
        os.remove(tempo_normal)
        shutil.rmtree(tempo)
        sys.exit(0)

    #
    # Reading the configuration file: "configuration.cfg"
    #
    config_files = [el for el in options.configuration_filename.split(",") if el and (os.path.isfile(el) or islink(el))]
    configfile = ''
    if config_files:
        configfile = config_files[0] # first one has priority
    confs = configuration.manage(configfile,skip_python=['openpyxl','xlrd'])
    if not options.data_directory:
        p = confs.get("DATA",None)
        if p and (os.path.isdir(p) or islink(p)):
            options.data_directory = p
        else:
            parser.error("ERROR: Argument '--data' needs to be specified as command line (or in 'configuration.cfg' file)!")
            sys.exit(1)
    # check if version of fusioncatcher.py matches the configuration.cfg file
    p = confs.get("FUSIONCATCHER",None)
    if p:
        t = parser.get_version()
        t = t.lower().split(".py")
        if t and len(t) == 2 and t[1].strip() == p.lower():
            pass
        else:
            print >>sys.stderr,"................................................................................"
            print >>sys.stderr,"ERROR: The version of configuration.cfg file does not match the version of the fusioncatcher.py!"
            print >>sys.stderr,"Please, fix this!"
            print >>sys.stderr,"................................................................................"
            sys.exit(1)

    #
    # DIRECTORIES
    #
    data_dir = adir(expand(options.data_directory))
    out_dir = adir(expand(options.output_directory))
    tmp_dir = adir(expand(options.tmp_directory))
    log_file = expand(outdir('fusioncatcher.log'))
    info_file = expand(outdir('info.txt'))


    ################################################################################
    # Managing OPTIONS
    ################################################################################

    if options.skip_mitochondrion_filtering:
        options.biotypes = options.biotypes.replace('mt','').replace(',,',',')

    if options.homolog == 0:
        options.biotypes = options.biotypes.replace('similar_reads','').replace(',,',',')

    if options.skip_adjacent and options.biotypes.find('adjacent') == -1:
        options.biotypes = options.biotypes + ',adjacent'

    if options.skip_banned_fusions:
        options.biotypes = options.biotypes.replace('banned','').replace(',,',',')
        options.biotypes = options.biotypes.replace('healthy','').replace(',,',',')

    if options.skip_blat:
        options.trim_blat = False
        options.psl_visualization = False
        options.skip_prefiltering_blat = False

    if options.keep_viruses:
        options.skip_viruses_filtering = False

    alg = set([el for el in options.aligners.lower().split(',') ])
    if (not options.skip_blat) and ('blat' not in alg):
        options.skip_blat = True
    if (not options.skip_star) and ('star' not in alg):
        options.skip_star = True
    if (not options.skip_bowtie2) and ('bowtie2' not in alg):
        options.skip_bowtie2 = True
    if (not options.skip_bwa) and ('bwa' not in alg):
        options.skip_bwa = True

    # create the output directory
    if (not os.path.isdir(out_dir)) and (not islink(out_dir)):
        if os.path.isfile(out_dir):
            print >>sys.stderr, "ERROR: The output directory is a actually a file! Please, delete it before proceeding further!"
            sys.exit(1)
        else:
            os.makedirs(out_dir)

    # deal with temporary files flag
    temp_flag = 'yes'
    if options.keep_temporary_files or (options.hash and options.hash != 'no'):
        temp_flag = 'no'

    ################################################################################
    # SENSITIVE
    ################################################################################
    sensitive = 0
    if options.sensitive:
        options.spanning_pairs = "2,2,2,2,2"
        options.spanning_reads = "2,2,2,2,2"
        options.length_anchor = "14,17,17,17,17"
        options.length_anchor2 = 40
        sensitive = sensitive + 1

    ################################################################################
    # MILD SENSITIVE
    ################################################################################
    if options.mildly_sensitive:
        options.spanning_pairs = "2,2,2,2,2"
        options.spanning_reads = "2,2,2,2,2"
        options.length_anchor = "13,15,15,15,15"
        options.length_anchor2 = 22
        options.mismatches = 2
        options.mismatches_psl = 4
        sensitive = sensitive + 1

    ################################################################################
    # HIGHLY SENSITIVE
    ################################################################################
    if options.highly_sensitive:
        options.spanning_pairs = "2,2,2,2,2"
        options.spanning_reads = "1,1,1,1,1"
        options.length_anchor = "13,14,14,14,14"
        options.length_anchor2 = 22
        options.mismatches = 2
        options.mismatches_psl = 4
        options.skip_prefiltering_psl = True
        sensitive = sensitive + 1

    ################################################################################
    # PARANOID SENSITIVE
    ################################################################################
    if options.paranoid_sensitive:
        options.spanning_pairs = "2,2,2,2,2"
        options.spanning_reads = "1,1,1,1,1"
        options.length_anchor = "11,11,11,11,11"
        options.length_anchor2 = 22
        options.mismatches = 2
        options.mismatches_psl = 4
        options.spanning_pairs_count = 20000
        options.homolog = 0
        options.skip_prefiltering_psl = True
#        options.all_reads_junction = True
        sensitive = sensitive + 1

    if sensitive > 1:
        parser.error("ERROR: The command line options: '--paranoid-sensitive','--sensitive','--midly-sensitive','--highly-sensitive' are mutually exclusive!")
        sys.exit(1)


    ################################################################################

    spanning_pairs = options.spanning_pairs.split(',')
    if len(spanning_pairs) != 5:
        print >>sys.stderr, "ERROR: the command option SPANNING_PAIRS has been given incorrectly! Expecting 4 commas!"
        sys.exit(1)
    spanning_pairs_bowtie = int(spanning_pairs[0])
    spanning_pairs_blat = int(spanning_pairs[1])
    spanning_pairs_star = int(spanning_pairs[2])
    spanning_pairs_bowtie2 = int(spanning_pairs[3])
    spanning_pairs_bwa = int(spanning_pairs[4])
    spanning_pairs_minimum = min(map(int,spanning_pairs))



    spanning_reads = options.spanning_reads.split(',')
    if len(spanning_reads) != 5:
        print >>sys.stderr, "ERROR: the command option SPANNING_READS has been given incorrectly! Expecting 4 commas!"
        sys.exit(1)
    spanning_reads_bowtie = int(spanning_reads[0])
    spanning_reads_blat = int(spanning_reads[1])
    spanning_reads_star = int(spanning_reads[2])
    spanning_reads_bowtie2 = int(spanning_reads[3])
    spanning_reads_bwa = int(spanning_reads[4])
    spanning_reads_minimum = min(map(int,spanning_reads))



    length_anchor = options.length_anchor.split(',')
    if len(length_anchor) != 5:
        print >>sys.stderr, "ERROR: the command option LENGTH_ANCHOR has been given incorrectly! Expecting 4 commas!"
        sys.exit(1)
    length_anchor_bowtie = int(length_anchor[0])
    length_anchor_blat = int(length_anchor[1])
    length_anchor_star = int(length_anchor[2])
    length_anchor_bowtie2 = int(length_anchor[3])
    length_anchor_bwa = int(length_anchor[4])
    length_anchor_minimum = min(map(int,length_anchor))



    length_anchor2 = options.length_anchor2



    if spanning_reads_minimum < 1:
        parser.error("ERROR: The minimum value of the SPANNING_READS is 1 but the value %s was given!" % (options.spanning_reads,))
        sys.exit(1)

    if length_anchor_minimum < 10:
        parser.error("ERROR: The minimum value of the LENGTH_ANCHOR is 10 but the value %s was given!" % (options.length_anchor,))
        sys.exit(1)

    if length_anchor2 <= length_anchor_minimum:
        parser.error("ERROR: --anchor-fusion2 (%d) should be larger than anchor-fusion (%s)!" % (options.length_anchor2,options.length_anchor))
        sys.exit(1)

    if spanning_pairs_bowtie != spanning_pairs_minimum or spanning_pairs_bowtie < 2:
        parser.error("ERROR: The minimum value of the SPANNING_PAIRS should have values larger than 2 but the value %s was given!" % (options.spanning_pairs,))
        sys.exit(1)

    if options.skip_star_bowtie and (not is_optparse_provided(parser,'limit_star')):
        options.limit_star = 3 * (2**30)


    # test that the version of build data matches
    if empty(datadir('version.txt')) or ( not os.path.exists(datadir('transcripts.fa'))) or ( not os.path.exists(datadir('genes.txt'))):
        print >>sys.stderr,"\n\n"
        print >>sys.stderr,"ERROR: The build data not found in '%s'!" % (data_dir,)
        print >>sys.stderr,"\n\n"
        sys.exit(1)

    ################################################################################
    # Contacts the FusionCatcher server to check for a more recent version
    ################################################################################
    old = []
    if not options.skip_update_check:
        try:
            import urllib2
            serverversion  = urllib2.urlopen('http://fusioncatcher.blogsite.org/fusioncatcher-version.txt').read()
            if serverversion:
                serverversion = serverversion.splitlines()
                serverversion = serverversion[0].strip()
                version = parser.get_version()
                if serverversion != version:
                    old = ["",
                           "="*80,
                           "WARNING: This is an OLD version of FusionCatcher! There is a newer",
                           "         version available! Please, update to the newest version!",
                           ""
                           " - Current version:   %s" % (version.replace('fusioncatcher.py','').strip(),),
                           " - New version:       %s"% (serverversion.replace('fusioncatcher.py','').strip(),),
                           "="*80,
                           ""
                          ]
                    file(info_file,'a').writelines([el.rstrip('\r\n')+'\n' for el in old])
                    print >>sys.stderr,'\n'.join(old)
        except:
            pass

    ################################################################################
    # Initialize pipeline
    ################################################################################
    import workflow

    #
    job = workflow.pipeline(
            log_filename       = log_file,
            checksums_filename = options.checksums_filename,
            hash_library       = options.hash,
            start_step         = options.start_step)

    ##############################################################################
    # SAVE EXTRA INFORMATION
    ##############################################################################

    os.system("set +e") # make sure that the shell scripts are still executed even if there are errors


    # save version of FusionCatcher
    job.add('printf',kind='program')
    job.add(('"\n================================================\n'+
             'Sotware version: %s\n'+
             '================================================\n\n\n"') % (
             parser.get_version(),),
             kind='parameter')
    job.add('>',info_file,kind='output')
    job.run()
    job.add('printf',kind='program')
    job.add('"\nSotware version: %s\n"' % (parser.get_version(),), kind='parameter')
    job.add('>>',log_file,kind='output')
    job.run()

    # create the temporary directory
    if job.run():
        if not os.path.isdir(tmp_dir) and not islink(tmp_dir):
            os.makedirs(tmp_dir)

    # check options supported by SORT command
    sort_parallel = False
    r = os.system("sort --help | grep 'parallel' > "+outdir('sort_help.txt'))
    if (not r) and (not empty(outdir('sort_help.txt'))) and len(file(outdir('sort_help.txt'),'r').readlines()) == 1:
        sort_parallel = True
    delete_file(outdir('sort_help.txt'))
    # check options supported by SORT command
    sort_buffer = None
    r = os.system("sort --help | grep 'buffer-size' > "+outdir('sort_help.txt'))
    if (not r) and (not empty(outdir('sort_help.txt'))) and len(file(outdir('sort_help.txt'),'r').readlines()) == 1:
        sort_buffer = "80%"
    delete_file(outdir('sort_help.txt'))
    # check options suppported by SORT command
    sort_lzop_compress = False
    r = os.system("sort --help | grep 'compress-program' > "+outdir('sort_help.txt')+ " ; lzop --help 2>/dev/null | grep -i 'compress' > " +outdir('lzop_help.txt'))
    if (not r) and ((not empty(outdir('sort_help.txt'))) and len(file(outdir('sort_help.txt'),'r').readlines()) == 1 and
        (not empty(outdir('lzop_help.txt'))) and len(file(outdir('lzop_help.txt'),'r').readlines()) >= 1):
        sort_lzop_compress = True
    delete_file(outdir('sort_help.txt'))
    delete_file(outdir('lzop_help.txt'))
    # check options suppported by SORT command
    sort_gzip_compress = False
    r = os.system("sort --help | grep 'compress-program' > "+outdir('sort_help.txt')+ " ; gzip --help 2>/dev/null | grep -i 'compress' > " +outdir('gzip_help.txt'))
    if (not r) and ((not empty(outdir('sort_help.txt'))) and len(file(outdir('sort_help.txt'),'r').readlines()) == 1 and
        (not empty(outdir('gzip_help.txt'))) and len(file(outdir('gzip_help.txt'),'r').readlines()) >= 1):
        sort_gzip_compress = True
    delete_file(outdir('sort_help.txt'))
    delete_file(outdir('gzip_help.txt'))


    # disable any compression done by SORT ===> FASTER
    sort_lzop_compress = False
    sort_gzip_compress = False

    # check if PIGZ is installed
    pigz = False
    r = os.system("pigz --version 2>/dev/null")
    if (not r):
        pigz = True

    # check if PXZ is installed
    pxz = False
    r = os.system("pxz --version 2>/dev/null")
    if (not r):
        pxz = True

    # save version of ENSEMBL used to analyse this data
    info(job,
         fromfile = datadir('version.txt'),
         tofile = info_file,
         top = ["===================================",
                "GENOME INFORMATION:",
                "==================================="],
         bottom = "\n")
#    job.add('printf',kind='program')
#    job.add('"\n============\nGENOME:\n============\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',datadir('version.txt'),kind='input')
#    job.add('>>',info_file,kind='output')
#    job.run()

    # save the genome fasta files used
    info(job,
         fromfile = datadir('genome_information.txt'),
         tofile = info_file,
         top = ["===================================",
                "Genome FASTA files:",
                "==================================="],
         bottom = "\n\n\n")
#    job.add('printf',kind='program')
#    job.add('"\nGenome FASTA files:\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',datadir('genome_information.txt'),kind='input')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    if not os.system("cat /etc/issue 2>&1 >/dev/null"):
        job.add('printf',kind='program')
        job.add('"\nLinux:\n------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('cat',kind='program')
        job.add('/etc/issue',kind='parameter')
        job.add('2>&1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

    # save version of Python used to analyze this data
    job.add('printf',kind='program')
    job.add('"\nPython:\n------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('python_version.py',kind='program')
    job.add('>>',info_file,kind='output')
    job.run(error_message=("Please, check if 'Python' (from "+
        "<http://python.org/>) is installed (or if 'configuration.cfg' file "+
        "is set up correctly)!"))

    # save version of BioPython used to analyze this data
    job.add('printf',kind='program')
    job.add('"\nBioPython:\n------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('biopython_version.py',kind='program')
    job.add('>>',info_file,kind='output')
    job.run(error_message=("Please, check if 'BioPython' (from "+
        "<http://biopython.org/>) is installed (or if 'configuration.cfg' file "+
        "is set up correctly)!"))

    # save version of BOWTIE used to analyze this data
    job.add('printf',kind='program')
    job.add('"\n===========\nTOOLS:\n===========\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('printf',kind='program')
    job.add('"\nSORT:\n------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('sort',kind='program')
    job.add('--version',kind='parameter')
    job.add('|',kind='parameter')
    job.add('head','-1',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('printf',kind='program')
    job.add('"\nBOWTIE:\n------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('bowtie',kind='program')
    job.add('--version',kind='parameter')
    job.add('|',kind='parameter')
    job.add('head','-2',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run(error_message=("Please, check if 'Bowtie' (from "+
        "<http://bowtie-bio.sourceforge.net/index.shtml>) is installed and it "+
        "is in the corresponding PATH!"))

    job.add('printf',kind='program')
    job.add('"\nPIGZ:\n------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    if pigz:
        job.add('pigz',kind='program')
        job.add('--version',kind='parameter')
        #job.add('2>&1',kind='parameter')
        job.add('2>>',info_file,kind='output')
        job.run()

    job.add('printf',kind='program')
    job.add('"\nPXZ:\n------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    if pxz:
        job.add('pxz',kind='program')
        job.add('--version',kind='parameter')
        job.add('2>&1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

    job.add('printf',kind='program')
    job.add('"\nfastq-dump (from SRA Toolkit):\n--------------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    if not os.system("which fastq-dump 2>&1 >/dev/null"):
        job.add('fastq-dump',kind='program')
        job.add('2>&1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('tail','-2',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

    job.add('printf',kind='program')
    job.add('"\nGNU Parallel:\n--------------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    if not os.system("which parallel 2>&1 >/dev/null"):
        job.add('parallel',kind='program')
        job.add('--version',kind='parameter')
        job.add('2>&1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

    job.add('printf',kind='program')
    job.add('"\nSAMTools:\n---------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    if not os.system("which samtools 2>&1 >/dev/null"):
        job.add('samtools',kind='program')
        job.add('2>&1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-3',kind='parameter')
        job.add('|',kind='parameter')
        job.add('tail','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

    job.add('printf',kind='program')
    job.add('"\nliftOver:\n---------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    if not os.system("which liftOver 2>&1 >/dev/null"):
        job.add('liftOver',kind='program')
        job.add('2>&1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

    job.add('printf',kind='program')
    job.add('"\nSeqTK:\n---------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('seqtk',kind='program')
    job.add('2>&1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('head','-3',kind='parameter')
    job.add('|',kind='parameter')
    job.add('tail','-1',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    job.add('printf',kind='program')
    job.add('"\nsed:\n---------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('sed',kind='program')
    job.add('--version',kind='parameter')
    job.add('2>&1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('head','-1',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()


    job.add('printf',kind='program')
    job.add('"\nawk:\n---------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('awk',kind='program')
    job.add('--version',kind='parameter')
    job.add('2>&1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('head','-1',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    if not options.skip_blat:
        # save version of BLAT used to analyze the data
        job.add('printf',kind='program')
        job.add('"\nBLAT:\n------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('blat',kind='program')
        job.add('|',kind='parameter')
        job.add('head','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if BLAT (from "+
            "<http://users.soe.ucsc.edu/~kent/src/> and "+
            "<http://hgdownload.cse.ucsc.edu/admin/exe/>) is installed correctly and it "+
            "is in the corresponding PATH (or if 'configuration.cfg' file is "+
            "set up correctly)!\n If there is no wish to use BLAT aligner then please "+
            "(re)run FusionCatcher using command line option '--skip-blat'.\n"+
            "Please, also read its commercial "+
            "license <http://www.kentinformatics.com/> if this applies in your case!"))

        # save version of faToTwoBit from BLAT used to analyze the data
        job.add('printf',kind='program')
        job.add('"\nfaToTwoBit (from BLAT toolbox):\n----------------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('faToTwoBit',kind='program')
        job.add('2>&1',kind='parameter')
        job.add('>','/dev/null',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if faToTwoBit (from BLAT toolbox, here: "+
            "<http://users.soe.ucsc.edu/~kent/src/> and "+
            "<http://hgdownload.cse.ucsc.edu/admin/exe/>) is installed correctly and it "+
            "is in the corresponding PATH (or if 'configuration.cfg' file is "+
            "set up correctly)!\n If there is no wish to use BLAT aligner then please "+
            "(re)run FusionCatcher using command line option '--skip-blat'.\n"+
            "Please, also read its commercial "+
            "license <http://www.kentinformatics.com/> if this applies in your case!"))

    if not options.skip_star:
        # save version of BLAT used to analyze the data
        job.add('printf',kind='program')
        job.add('"\nSTAR:\n------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('STAR',kind='program')
        job.add('--version',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if STAR (from "+
            "<https://code.google.com/p/rna-star/> and "+
            "<http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARreleases/Patches/>) is installed correctly and it "+
            "is in the corresponding PATH (or if 'configuration.cfg' file is "+
            "set up correctly)!\n If there is no wish to use STAR aligner then please "+
            "(re)run FusionCatcher using command line option '--skip-star'."))

    if not options.skip_bowtie2:
        # save version of BOWTIE2 used to analyze the data
        job.add('printf',kind='program')
        job.add('"\n\nBOWTIE2:\n------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('bowtie2',kind='program')
        job.add('--version',kind='parameter')
        job.add('2>&1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-2',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if BOWTIE2 (from "+
            "<http://bowtie-bio.sourceforge.net/bowtie2/index.shtml> "+
            "is installed correctly and it "+
            "is in the corresponding PATH (or if 'configuration.cfg' file is "+
            "set up correctly)!\n If there is no wish to use BOWTIE2 aligner then please "+
            "(re)run FusionCatcher using command line option '--skip-bowtie2'."))

    if not options.skip_bwa:
        # save version of BWA used to analyze the data
        job.add('printf',kind='program')
        job.add('"\n\nBWA:\n------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('bwa',kind='program')
        job.add('2>&1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-3',kind='parameter')
        job.add('|',kind='parameter')
        job.add('tail','-1',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if BWA (from "+
            "<http://bio-bwa.sourceforge.net/> "+
            "is installed correctly and it "+
            "is in the corresponding PATH (or if 'configuration.cfg' file is "+
            "set up correctly)!\n If there is no wish to use BWA aligner then please "+
            "(re)run FusionCatcher using command line option '--skip-bwa'."))


    if options.assembly:
        # save version of VELVET used to analyze the data
        job.add('printf',kind='program')
        job.add('"\nVELVET:\n------\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
        job.add('velvetg',kind='program')
        job.add('|',kind='parameter')
        job.add('head','-2',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if VELVETT (from "+
            "<http://www.ebi.ac.uk/~zerbino/velvet/>)"+
            "is installed and it "+
            "is in the corresponding PATH!"))
        job.add('velveth',kind='program')
        job.add('|',kind='parameter')
        job.add('head','-2',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run(error_message = ("Please, check if VELVETT (from "+
            "<http://www.ebi.ac.uk/~zerbino/velvet/>)"+
            "is installed and it "+
            "is in the corresponding PATH!"))

    # save the command line arguments again
    job.add('printf',kind='program')
    job.add('"\n\n\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    if job.run():
        clo = [ el + ' = ' + str(getattr(options,el))+'\n' for el in dir(options)
                if (not el.startswith('_')) and
                    type(getattr(options,el)).__name__ in ('str','bool','int','float')]
        clo.append("main_script = " + expand(sys.argv[0])+'\n')
        clo.append("main_script_version = " + version.replace('%prog','').strip()  +'\n')
        clo.insert(0,"Pipeline parameters:\n")
        clo.insert(1,"====================\n")
        clo.append("\n")
        clo.append("Current working directory:\n---------------------------\n%s\n" % (expand(os.getcwd()),))
        clo.append("\n")
        clo.append("Command line used for launching FusionCatcher:\n")
        clo.append("----------------------------------------------\n")
        clo.append("%s\n" % (' \\ \n'.join(sys.argv),))
        clo.append("----------------------------------------------\n")
        clo.append("\n\n\n")
        clo.append("Shebang for Python scripts:\n")
        clo.append("---------------------------\n")
        clo.append(file(expand(sys.argv[0]),'r').readline()+"\n")
        clo.append("\n\n\n")
        file(info_file,'a').writelines(clo)

        clo = []
        clo.append("#!/usr/bin/env bash\n")
        clo.append("cd '%s'\n" % (os.getcwd(),))
        clo.append("\n")
        clo.append("%s\n" % (' \\\n'.join(sys.argv),))
        clo.append("\n")
        file(outdir('restart.sh'),'w').writelines(clo)

        os.system("chmod u+x '%s'" % (outdir('restart.sh'),))

    info(job,
         fromfile = configfile,
         tofile = info_file,
         top = ["===================",
                "CONFIGURATION.CFG:",
                "==================="],
         bottom = "\n\n\n")

    if islink(data_dir):
        job.add('printf',kind='program')
        job.add('"\n============\nDATA DIRECTORY:\n============\n%s\nIt links to:\n%s\n\n\n"' % (data_dir,expand(os.readlink(data_dir[:-1] if data_dir.endswith(os.sep) else data_dir))),kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()
    else:
        job.add('printf',kind='program')
        job.add('"\n============\nDATA DIRECTORY:\n============\n%s\nIt is not a link!\n\n\n"' % (data_dir,),kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()


    if job.run():
        # test that the version of build data matches
        build_version = file(datadir('version.txt'),'r').readline().strip()
        build_version = build_version.lower().split('.py')
        old_build_version = 'fusioncatcher-build.py 0.99.3d beta--------'
        old_build_version = old_build_version.lower().split('.py')
        pipeline_version = parser.get_version()
        pipeline_version = pipeline_version.lower().strip().split('.py')
        if len(pipeline_version) > 1 and len(build_version) > 1 and (build_version[1].strip() == pipeline_version[1].strip() or old_build_version[1].strip() == build_version[1].strip()):
            print "Version of the data build matches the version of pipeline version!"
        else:
            print >>sys.stderr,"...................."
            print >>sys.stderr,"ERROR: The version of the data build does not match the version of this pipeline!"
            print >>sys.stderr,"Please, run again the 'fusioncatcher-build.py' in order to fix this!"
            print >>sys.stderr,"...................."
            job.close()
            sys.exit(1)

    # find available memory
    job.add('printf',kind='program')
    job.add('"\n============\nMEMORY:\n============\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('free',kind='program')
    job.add('-m',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('printf',kind='program')
    job.add('"\nTotal installed RAM memory = %d MB"' % (total_memory,),kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('printf',kind='program')
    job.add('"\n\n\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    ##############################################################################
    # STARTING with the input
    ##############################################################################

    if options.input_filename.endswith(','):
        print >>sys.stderr,"ERROR: The input list of files is ending with a comma!"
        print >>sys.stderr,"Please, remove the comma from the end of the list of files "
        print >>sys.stderr,"or check that there are no blanks, tabs, or spaces after the comma"
        print >>sys.stderr,"and re-run FusionCatcher from the very beggining!"
        job.close()
        sys.exit(1)

    list_input_files = sorted([expand(el.strip('"').strip("'")) if el.find('://') == -1 else el.strip('"').strip("'") for el in options.input_filename.split(',') if el])

    f = []
    urls = []
    new_input_output = outdir('input/')
    for element in list_input_files:
        if os.path.isdir(element):
            print "The input '%s' is a directory..." % (element,)
            f.extend([os.path.join(element,el) for el in os.listdir(element)])
        elif element.find('http://') != -1 or element.find('https://') != -1 or element.find('ftp://') != -1 or element.find('ssh://') != -1:
            if not os.path.isdir(new_input_output):
                os.makedirs(new_input_output)

            if element.find('ssh://') != -1:
                #scp -r user@your.server.example.com:/path/to/foo /home/user/Desktop/
                job.add('scp',kind='program')
                job.add('-r',kind='parameter')
                job.add('',element.replace('ssh://',''),kind='parameter')
                job.add('',new_input_output,kind='output')
                job.run()
            else:
                job.add('wget',kind='program')
                job.add('--no-check-certificate',kind='parameter')
                if element.endswith('/') or element.endswith(os.sep):
                    job.add('-r',kind='parameter')
                    job.add('-nd',kind='parameter')
                    job.add('-np',kind='parameter')
                    job.add('-A','fq,fastq,fq.gz,fastq.gz,fq.zip,fastq.zip,fq.gz,fastq.bz2,sra,fastq.xz,fq.xz',kind='parameter')
                job.add('',element,kind='parameter')
                job.add('-P',new_input_output,kind='output')
                job.run()
        else:
            f.append(element)
    if os.path.isdir(new_input_output):
        urls1 = [os.path.join(new_input_output,el) for el in os.listdir(new_input_output) if os.path.isfile(os.path.join(new_input_output,el)) and not el.startswith('.')]
        f.extend(urls1)
        urls.extend(urls1)
        urls2 = [os.path.join(new_input_output,el) for el in os.listdir(new_input_output) if os.path.isdir(os.path.join(new_input_output,el)) and not el.startswith('.')]
        for elx in urls2:
            urls3 = [os.path.join(elx,el) for el in os.listdir(elx) if os.path.isfile(os.path.join(elx,el)) and not el.startswith('.')]
            f.extend(urls3)
            urls.extend(urls3)
    list_input_files = [el for el in f if is_known_extension(el) and (not empty(el)) and (not el.startswith('.'))]
    list_input_files = list(set(list_input_files))
    list_input_files.sort()

    job.write(["Input files (which contain the short reads):"]+list_input_files)

    info(job,
         fromfile = None,
         tofile = info_file,
         top = ["",
                "Input files (which contain the short reads):",
                "--------------------------------------------"]+list_input_files,
         bottom = "\n\n\n",
         temp_path=temp_flag)

#    job.add('printf',kind='program')
#    job.add('"\n============\nInput files:\n============\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    for el in list_input_files:
#        job.add('printf',kind='program')
#        job.add('"%s\n"' % (el,),kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    if len(list_input_files) < 1:
        job.write("ERROR: No input valid files have been found (given input is: '%s' )!" % (options.input_filename,),stderr=True)
        job.close()
        sys.exit(1)

    # protect the input files by accidental deletion
    job.protect(list_input_files)

    # handle the SRA files
    new_list_input_files = []
    i = 0
    for input_file in list_input_files:
        if input_file.endswith(".sra"):

            i = i + 1

            outfile1 = outdir(os.path.basename(input_file)[:-4]+'_1.fastq')
            outfile2 = outdir(os.path.basename(input_file)[:-4]+'_2.fastq')
            job.add('fastq-dump',kind='program')
            job.add('--split-files',kind='parameter')
            job.add('',input_file,kind='input')
            job.add('-O',out_dir,kind='output',checksum='no')
            job.add('',outfile1,kind='output',command_line='no')
            job.add('',outfile2,kind='output',command_line='no')
            if job.run(error_message=("Please, check if 'fastq-dump' (from NCBI SRA "+
                "Toolgit <http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>) "+
                "is installed and it is in the corresponding PATH!")):
                # test if there are indeed two output files
                if empty(outfile2):
                    job.write("ERROR: The SRA input file '%s' does not contain paired-end reads!" % (input_file,),stderr=True)
                    sys.exit(1)

            outfile3 = outdir(os.path.basename(input_file)[:-4]+'_1.fq')
            outfile4 = outdir(os.path.basename(input_file)[:-4]+'_2.fq')
            job.add('sra2illumina.py',kind='program')
            #job.add('--tag_read_name','Z'+str(i)+'Z',kind='parameter')
            job.add('--input_1',outfile1,kind='input',temp_path=temp_flag)
            job.add('--input_2',outfile2,kind='input',temp_path=temp_flag)
            job.add('--output_1',outfile3,kind='output')
            job.add('--output_2',outfile4,kind='output')
            job.add('--link','change',kind='parameter')
            job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
            job.run()

            new_list_input_files.append(outfile3)
            new_list_input_files.append(outfile4)

        elif input_file.endswith(".bam") or input_file.endswith(".sam"):

            i = i + 1

            outfile1 = outdir(os.path.basename(input_file)[:-4]+'_1.fq')
            outfile2 = outdir(os.path.basename(input_file)[:-4]+'_2.fq')

#            job.add('samtools',kind='program')
#            job.add('view',input_file,kind='input')
#            job.add('|',kind='parameter')
# ##            job.add('grep',kind='parameter')
# ##            job.add('-v','^@',kind='parameter')
# ##            job.add('|',kind='parameter')
#            job.add('LC_ALL=C',kind='parameter')
#            job.add('sort',kind='parameter')
#            job.add('-t',"'\t'",kind='parameter')
#            job.add('-k','1,1',kind='parameter')
#            if sort_buffer:
#                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
#            if sort_parallel:
#                job.add('--parallel',options.processes,kind='parameter',checksum='no')
#            if sort_lzop_compress:
#                job.add('--compress-program','lzop',kind='parameter',checksum='no')
#            elif sort_gzip_compress:
#                job.add('--compress-program','gzip',kind='parameter',checksum='no')
#            job.add('-T',tmp_dir,kind='parameter',checksum='no')
#            job.add('|',kind='parameter')
#            job.add('awk',kind='parameter')
#            job.add("""-F"\\\\t" '{if ( $1==old1 ) { print "@"$1"\\n"$10"\\n+\\n"$11 > "%s"; print "@"old1"\\n"old10"\\n+\\n"old11 > "%s"; old1="" }; {old1=$1; old10=$10; old11=$11}; }'""" % (outfile1,outfile2),kind='parameter')
# ##            job.add("""-F"\\\\t" '{if ( $1==old ) { print old0"\\n"$0; old="" }; {old=$1; old0=$0}; }'""",kind='parameter')
# ##            job.add('|',kind='parameter')
# ##            job.add('awk',kind='parameter')
# ##            job.add("""'{if(NR%%2==0) {print "@"$1"\\n"$10"\\n+\\n"$11 > "%s.fq"} else {print "@"$1"\\n"$10"\\n+\\n"$11 > "%s.fq"}}'""" % (outfile1,outfile2),kind='parameter')
#            job.add('',outfile1,kind='output',command_line='no')
#            job.add('',outfile2,kind='output',command_line='no')
#            job.run()
#            #cat H716.sam | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > 1.fq
#            #cat H716.sam | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > 2.fq
#            # java -jar SamToFastq.jar INPUT=G28616.NCI-H2228.1.bam F=r1.fq F2=r2.fq NON_PF=True
#

            job.add('java',kind='program')
            job.add('-jar',os.path.join(confs.get('PICARD',''),'picard.jar'),kind='parameter')
            job.add('SamToFastq',kind='parameter')
            job.add('NON_PF=','True',kind='parameter',space='no')
            job.add('INPUT=',input_file,kind='input',space='no')
            job.add('F=',outfile1,kind='output',space='no')
            job.add('F2=',outfile2,kind='output',space='no')
            job.run()

            new_list_input_files.append(outfile1)
            new_list_input_files.append(outfile2)


        else:
            new_list_input_files.append(input_file)


    list_input_files = new_list_input_files[:]


    shuffled = False
    new_list_input_files = []
    for i,input_file in enumerate(list_input_files):
        output_file = None
        if input_file.lower().endswith('.gz') and (not input_file.lower().endswith('.tar.gz')):
            output_file = outdir('init-'+str(i)+'_'+os.path.basename(input_file)[:-3])
            # decompress
            if pigz:
                job.add('pigz',kind='program')
                job.add('-p',options.processes,kind='parameter',checksum='no')
            else:
                job.add('gzip',kind='program')
            job.add('-d',kind='parameter')
            job.add('-f',kind='parameter')
            job.add('-c',kind='parameter')
            job.add('',input_file,kind='input')
            job.add('>',output_file,kind='output')
            job.run()
        elif input_file.lower().endswith('.xz'):
            output_file = outdir('init-'+str(i)+'_'+os.path.basename(input_file)[:-3])
            # decompress
            # decompress
            if pxz:
                job.add('pxz',kind='program')
                job.add('-T',options.processes,kind='parameter',checksum='no')
            else:
                job.add('xz',kind='program')
            job.add('-dc',kind='parameter')
            job.add('',input_file,kind='input')
            job.add('>',output_file,kind='output')
            job.run()
        elif input_file.lower().endswith('.zip'):
            output_file = outdir('init-'+str(i)+'_'+os.path.basename(input_file)[:-4])
            # decompress
            job.add('unzip',kind='program')
            job.add('-p',kind='parameter')
            job.add('-o',kind='parameter')
            job.add('',input_file,kind='input')
            job.add('>',output_file,kind='output')
            job.run()
        elif input_file.lower().endswith('.bz2') and (not input_file.lower().endswith('.tar.bz2')):
            output_file = outdir('init-'+str(i)+'_'+os.path.basename(input_file)[:-4])
            # decompress
            job.add('bzip2',kind='program')
            job.add('-d',kind='parameter')
            job.add('-f',kind='parameter')
            job.add('-c',kind='parameter')
            job.add('',input_file,kind='input')
            job.add('>',output_file,kind='output')
            job.run()
        elif ( input_file.lower().endswith('.fq') or
               input_file.lower().endswith('.fastq')):
            output_file = outdir('init-'+str(i)+'_'+os.path.basename(input_file))
            # link
            job.link(input_file, output_file, temp_path=temp_flag)
        else:
            print >> sys.stderr,"ERROR: unknown extension of the input file! '%s'" % (input_file,)
            print >> sys.stderr,"Supported extension files are:"
            print >> sys.stderr," .fastq,"
            print >> sys.stderr," .fq,"
            print >> sys.stderr," .txt.zip,"
            print >> sys.stderr," .fastq.zip,"
            print >> sys.stderr," .fq.zip,"
            print >> sys.stderr," .txt.gz,"
            print >> sys.stderr," .fastq.gz,"
            print >> sys.stderr," .fq.gz,"
            print >> sys.stderr," .txt.bz2,"
            print >> sys.stderr," .fastq.bz2,"
            print >> sys.stderr," .fq.bz2,"
            print >> sys.stderr," .fastq.xz,"
            print >> sys.stderr," .fq.xz,"
            print >> sys.stderr," .bam,"
            print >> sys.stderr," .sam"
            sys.exit(1)

        new_list_input_files.append(output_file)
        
    if options.single_end and list_input_files:
        # single-end reads
        new_list_input_files = []
        for i,input_file in enumerate(list_input_files):
            output_file_1 = outdir('single-'+str(i)+'_'+os.path.basename(input_file))
            output_file_2 = outdir('single-'+str(i)+'_'+os.path.basename(input_file))

            # compute the read lengths for the input file
            job.add('lengths_reads.py',kind='program')
            job.add('--input',input_file,kind='input')
            job.add('--output',outdir('single','log_lengths_single_reads_%s.txt' % (str(i),)),kind='output')
            job.run()
            
            max_len_reads = 0
            if os.path.exists(outdir('single','log_lengths_single_reads_%s.txt' % (str(i),))):
                max_len_reads = int(file(outdir('single','log_lengths_single_reads_%s.txt' % (str(i),)),'r').readline().rstrip())
                if not options.bridges:
                    options.bridges = int(math.ceil(float(max_len_reads)/float(160)))
                if max_len_reads < options.sonication:
                    print >> sys.stderr,"ERROR: The input single-end reads are too short to be used by FusionCatcher!"
                    sys.exit(1)

            job.add('seqtk',kind='program')
            job.add('trimfq',kind='parameter')
            job.add('-q','0.25',kind='parameter')
            job.add('',input_file,kind='input',temp_path=temp_flag)
            job.add('>',outdir('single-read--%d.fq' %(i,)),kind='output')
            job.run()



            job.add('fragment_fastq.py',kind='program')
            job.add('-1',outdir('single-read--%d.fq' %(i,)),kind='input',temp_path=temp_flag)
            job.add('-2','-',kind='parameter')
            job.add('-f',output_file_1,kind='output')
            job.add('-r',output_file_2,kind='output')
            job.add('--window-size',options.trim_psl_3end_keep,kind='parameter')
            job.add('--step-size',options.trim_psl_3end_keep-length_anchor_minimum+1,kind='parameter')
            job.add('--threshold-read',options.trim_psl_3end_keep + 10,kind='parameter')
            job.add('--anchors',options.bridges,kind='parameter')
            job.add('--skip-short',options.trim_3end_keep,kind='parameter')
            job.add('--trim-n',kind='parameter')
            job.run()
            
            new_list_input_files.append(output_file_1)
            new_list_input_files.append(output_file_2)
            
        list_input_files = new_list_input_files[:]
        

    nif = len(list_input_files)
    if nif == 0:
        print >> sys.stderr,"ERROR: No input files found!"
        sys.exit(1)
    elif nif == 1 or options.skip_interleave_processing:
        pass
    else:
        list_input_files = new_list_input_files[:]
        new_list_input_files = []
        pairs = [(list_input_files[i-1],list_input_files[i]) for i in xrange(1,len(list_input_files),2)]

        i = -1
        final_i = i
        for (f,r) in pairs:
            i = i + 1
            # automatically remove adapters
            output_1_file = outdir(os.path.basename(f).replace('init-','init-noadapt-'))
            output_2_file = outdir(os.path.basename(r).replace('init-','init-noadapt-'))

            job.add('printf',kind='program')
            job.add('"\n\nFirst 8 lines of input FASTQ file: %s\n-------------------------\n"' % (r,),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()

            job.add('head',kind='program')
            job.add('-8',r,kind='input')
            job.add('>>',info_file,kind='output')
            job.run()


            job.add('printf',kind='program')
            job.add('"\n\nFirst 8 lines of input FASTQ file: %s\n-------------------------\n"' % (f,),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()

            job.add('head',kind='program')
            job.add('-8',f,kind='input')
            job.add('>>',info_file,kind='output')
            job.run()


            job.add('printf',kind='program')
            job.add('"\n\nLast 8 lines of input FASTQ file: %s\n-------------------------\n"' % (r,),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()

            job.add('tail',kind='program')
            job.add('-8',r,kind='input')
            job.add('>>',info_file,kind='output')
            job.run()


            job.add('printf',kind='program')
            job.add('"\n\nLast 8 lines of input FASTQ file: %s\n-------------------------\n"' % (f,),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()

            job.add('tail',kind='program')
            job.add('-8',f,kind='input')
            job.add('>>',info_file,kind='output')
            job.run()

            job.add('overlap.py',kind='program')
            job.add('--input_1',f,kind='input')
            job.add('--input_2',r,kind='input')
            job.add('--processes',options.processes,kind='parameter',checksum='no')
            job.add('--fail-gracefully',kind='parameter',checksum='no')
            job.add('--output',outdir('log_overlaps__%d.txt' % (i,)),kind='output')
            job.add('2>',outdir('log_overlaps_error__%d.txt' % (i,)),kind='parameter',checksum='no')
            job.run()

            sizes = False
            if os.path.exists(outdir('log_overlaps_error__%d.txt' % (i,))):
                # find out if it failed because of the different reads sizes
                sizes = True if [1 for line in file(outdir('log_overlaps_error__%d.txt' % (i,)),'r').readlines() if line.lower().find('different lengths') != -1] else False

            if job.iff(sizes, id = "#different_sizes_of_reads_overlap_%d#" % (i,)):
                # input FASTQ files contain reads of different sizes therefore they are paded so end up being the same size

                job.clean(outdir('log_overlaps_error__%d.txt' % (i,)),temp_path=temp_flag)
                job.clean(outdir('log_overlaps__%d.txt' % (i,)),temp_path=temp_flag)

                job.add('lengths_reads.py',kind='program')
                job.add('--input',f,kind='input')
                job.add('--output',outdir('log_lengths_original_reads_f_%d.txt' % (i,)),kind='output')
                job.run()

                job.add('lengths_reads.py',kind='program')
                job.add('--input',r,kind='input')
                job.add('--output',outdir('log_lengths_original_reads_r_%d.txt' % (i,)),kind='output')
                job.run()

                #job.add('cat',kind='program')
                #job.add('',outdir('log_lengths_original_reads_f_%d.txt' % (i,)),kind='input',temp_path=temp_flag)
                #job.add('',outdir('log_lengths_original_reads_r_%d.txt' % (i,)),kind='input',temp_path=temp_flag)
                #job.add('|',kind='parameter')
                #job.add('LC_ALL=C',kind='parameter')
                #job.add('sort',kind='parameter')
                #if sort_buffer:
                #    job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                #if sort_parallel:
                #    job.add('--parallel',options.processes,kind='parameter',checksum='no')
                #if sort_lzop_compress:
                #    job.add('--compress-program','lzop',kind='parameter',checksum='no')
                #elif sort_gzip_compress:
                #    job.add('--compress-program','gzip',kind='parameter',checksum='no')
                #job.add('-T',tmp_dir,kind='parameter',checksum='no')
                #job.add('-nr',kind='parameter')
                #job.add('|',kind='parameter')
                #job.add('uniq',kind='parameter')
                #job.add('>',outdir('log_lengths_original_reads_fr_%d.txt' % (i,)),kind='output')
                #job.run()


                x = outdir('log_lengths_original_reads_f_%d.txt' % (i,))
                y = outdir('log_lengths_original_reads_r_%d.txt' % (i,))

                if job.run():
                    n1 = 0
                    if os.path.isfile(x):
                        n1 = int(file(x,'r').readline().strip())
                    n2 = 0
                    if os.path.isfile(y):
                        n2 = int(file(y,'r').readline().strip())

                    file(outdir('log_lengths_original_reads_fr_%d.txt' % (i,)),'w').write(str(min(n1,n2)))

                job.clean(x,temp_path=temp_flag)
                job.clean(y,temp_path=temp_flag)

                ff = f[:]
                rr = r[:]
                f = outdir(os.path.basename(f).replace('init-','init-f-'))
                r = outdir(os.path.basename(r).replace('init-','init-r-'))

                # add A to the reads which are shorter in order to make all reads have the same length
                job.add('padding-fastq.py',kind='program')
                job.add('--input',ff,kind='input',temp_path=temp_flag)
                job.add('--output',f,kind='output')
                job.add('--size',outdir('log_lengths_original_reads_fr_%d.txt' % (i,)),kind='parameter',from_file='yes')
                job.run()
                job.add('padding-fastq.py',kind='program')
                job.add('--input',rr,kind='input',temp_path=temp_flag)
                job.add('--output',r,kind='output')
                job.add('--size',outdir('log_lengths_original_reads_fr_%d.txt' % (i,)),kind='parameter',from_file='yes')
                job.run()

                job.clean(outdir('log_lengths_original_reads_fr_%d.txt' % (i,)),temp_path=temp_flag)

                job.add('overlap.py',kind='program')
                job.add('--input_1',f,kind='input')
                job.add('--input_2',r,kind='input')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--output',outdir('log_overlaps_fr_%d.txt' % (i,)),kind='output')
                job.run()

                info(job,
                     fromfile = outdir('log_overlaps_fr_%d.txt' % (i,)),
                     tofile = info_file,
                     top = ["","","",
                            "Pair-reads overlappings (after padding):",
                            "----------------------------------------"],
                     bottom = "\n\n\n",
                     temp_path=temp_flag)

            else:
                info(job,
                     fromfile = outdir('log_overlaps__%d.txt' % (i,)),
                     tofile = info_file,
                     top = ["","","",
                            "Pair-reads overlappings:",
                            "------------------------"],
                     bottom = "\n\n\n",
                     temp_path=temp_flag)

                job.clean(outdir('log_overlaps_error__%d.txt' % (i,)),temp_path=temp_flag)


            if not options.skip_adapter_filtering:

    #            job.add('printf',kind='program')
    #            job.add('"\nPair-reads overlappings\n-------------------------\n"',kind='parameter')
    #            job.add('>>',info_file,kind='output')
    #            job.run()
    #            job.add('cat',kind='program')
    #            job.add('',outdir('log_overlaps.txt'),kind='input',temp_path=temp_flag)
    #            job.add('>>',info_file,kind='output')
    #            job.run()


                job.add('remove_adapter.py',kind='program')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--input_1',f,kind='input',temp_path=temp_flag)
                job.add('--input_2',r,kind='input',temp_path=temp_flag)
                job.add('--output_1',output_1_file,kind='output')
                job.add('--output_2',output_2_file,kind='output')
                job.add('--trim-n',options.mismatches+1,kind='parameter')
                job.add('--link','hard',kind='parameter',checksum='no')
                job.add('>',outdir('log_adapters_%d.txt' % (i,)),kind='parameter',checksum='no')
                job.add('2>&1',kind='parameter',checksum='no')
                job.run()

                info(job,
                     fromfile = outdir('log_adapters_%d.txt' % (i,)),
                     tofile = info_file,
                     top = ["",
                            "Adapters information:",
                            "--------------------"],
                     bottom = "\n\n\n",
                     temp_path = temp_flag)

#                job.add('printf',kind='program')
#                job.add('"\nAdapters\n---------\n"',kind='parameter')
#                job.add('>>',info_file,kind='output')
#                job.run()
#                job.add('cat',kind='program')
#                job.add('',outdir('log_adapters.txt'),kind='input',temp_path=temp_flag)
#                job.add('>>',info_file,kind='output')
#                job.run()
            else:
                job.link(f,output_1_file,temp_path=temp_flag)
                job.link(r,output_2_file,temp_path=temp_flag)

            in1 = output_1_file
            in2 = output_2_file
            ou3 = outdir(os.path.basename(in1).replace('init-','init-clear-'))
            ou4 = outdir(os.path.basename(in2).replace('init-','init-clear-'))
            # remove the reads marked as bad by Illumina
            job.add('remove-bad-illumina.py',kind='program')
            job.add('--input',in1,kind='input',temp_path=temp_flag)
            job.add('--output',ou3,kind='output')
            job.add('--link','hard',kind='parameter',checksum='no')
            job.add('2>',outdir('log_bad_illumina_1_%d.txt' % (i,)),kind='parameter',checksum='no')
            job.run()

            info(job,
                fromfile = outdir('log_bad_illumina_1_%d.txt' % (i,)),
                tofile = info_file,
                top = ["Reads (mate 1 from pair) removed because being marked as bad by Illumina:",
                       "-------------------------------------------------------------------------"],
                bottom = "\n\n\n",
                temp_path=temp_flag)

#            job.add('printf',kind='program')
#            job.add(('"\nReads (mate 1 reads) removed because being marked as bad by Illumina:\n'+
#                     '---------------------------------------------------------------------\n"'),kind='parameter')
#            job.add('>>',info_file,kind='output')
#            job.run()
#            job.add('cat',kind='program')
#            job.add('',outdir('log_bad_illumina_1.txt'),kind='input',temp_path=temp_flag)
#            job.add('>>',info_file,kind='output')
#            job.run()
#            job.add('printf',kind='program')
#            job.add('"\n\n\n"',kind='parameter')
#            job.add('>>',info_file,kind='output')
#            job.run()


            job.add('remove-bad-illumina.py',kind='program')
            job.add('--input',in2,kind='input',temp_path=temp_flag)
            job.add('--output',ou4,kind='output')
            job.add('--link','hard',kind='parameter',checksum='no')
            job.add('2>',outdir('log_bad_illumina_2_%d.txt' % (i,)),kind='parameter',checksum='no')
            job.run()

            info(job,
                 fromfile = outdir('log_bad_illumina_2_%d.txt' % (i,)),
                 tofile = info_file,
                 top = ["Reads (mate 2 from pair) removed because being marked as bad by Illumina:",
                        "-------------------------------------------------------------------------"],
                 bottom = "\n\n\n",
                 temp_path=temp_flag)

#            job.add('printf',kind='program')
#            job.add(('"Reads (read 2 from pair) removed because being marked as bad by Illumina:\n'+
#                     '----------------------------------------------------------------------------\n"'),kind='parameter')
#            job.add('>>',info_file,kind='output')
#            job.run()
#            job.add('cat',kind='program')
#            job.add('',outdir('log_bad_illumina_2.txt'),kind='input',temp_path=temp_flag)
#            job.add('>>',info_file,kind='output')
#            job.run()
#            job.add('printf',kind='program')
#            job.add('"\n\n\n"',kind='parameter')
#            job.add('>>',info_file,kind='output')
#            job.run()

            # add /1 and /2 in the end of the reads ids
            in1 = ou3
            in2 = ou4
            ou3 = outdir(os.path.basename(in1).replace('init-','init-sra-'))
            ou4 = outdir(os.path.basename(in2).replace('init-','init-sra-'))
            job.add('sra2illumina.py',kind='program')
            #job.add('--tag_read_name','Z'+str(i)+'Z',kind='parameter')
            job.add('--input_1',in1,kind='input',temp_path=temp_flag)
            job.add('--input_2',in2,kind='input',temp_path=temp_flag)
            job.add('--output_1',ou3,kind='output')
            job.add('--output_2',ou4,kind='output')
            job.add('--link','hard',kind='parameter')
            job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
            job.run()

            in1 = ou3
            in2 = ou4
            output_file = outdir(os.path.basename(in1).replace('init-','init-shuffle-').replace("_1.",".").replace("_2.","."))
#            job.add('shuffle.py',kind='program')
#            job.add('--input_1',in1,kind='input',temp_path=temp_flag)
#            job.add('--input_2',in2,kind='input',temp_path=temp_flag)
#            job.add('--output',output_file,kind='output')
#            job.run()
#            job.add('awk',kind='program')
#            job.add("""'{print; getline; print; getline; print; getline; print; getline < "%s"; print; getline < "%s"; print;  getline < "%s"; print;  getline < "%s"; print}' '%s'""" % (in2,in2,in2,in2,in1),kind='parameter')
#            job.add('',in1,kind='input',temp_path=temp_flag,command_line='no')
#            job.add('',in2,kind='input',temp_path=temp_flag,command_line='no')
#            job.add('>',output_file,kind='output')
#            job.run()
            #awk '{print; getline; print; getline; print; getline; print;  getline < "2.txt"; print;  getline < "2.txt"; print;  getline < "2.txt"; print;  getline < "2.txt"; print}' 1.txt
            job.add('seqtk',kind='program')
            job.add('mergepe',kind='parameter')
            job.add('',in1,kind='input',temp_path=temp_flag)
            job.add('',in2,kind='input',temp_path=temp_flag)
            job.add('>',output_file,kind='output')
            job.run()
            shuffled = True

            new_list_input_files.append(output_file)

        job.add('printf',kind='program')
        job.add('"\n\n\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

        # last file in case that the number of input files is odd
        if len(list_input_files) % 2 == 1:
            f = list_input_files[-1]
            output_file = outdir(os.path.basename(f).replace('init-','init-single-'))
            job.link(f,output_file,temp_path=temp_flag)
            new_list_input_files.append(output_file)


    job.add('printf',kind='program')
    job.add('"\nTotal Count of reads (from all FASTQ files given as input and before any read removal is done, i.e. quality filtering, pre-processing):\n--------------\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    list_input_files = new_list_input_files[:]
    new_list_input_files = []
    for i,input_file in enumerate(list_input_files):

        job.add('printf',kind='program')
        job.add('"%s = "' % (input_file,),kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

        job.add('cat',kind='program')
        job.add('',input_file,kind='input')
        job.add('|',kind='parameter')
        job.add("echo $((`wc -l`/4))",kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

        # convert the read names to Illumina Solexa version 1.5 format (i.e. end in /1 or /2)
        output_file = outdir(os.path.basename(input_file).replace('init-','init-head-'))
        job.add('solexa18to15.py',kind='program')
        job.add('--fail',kind='parameter')
        job.add('--input',input_file,kind='input',temp_path=temp_flag)
        job.add('--output',output_file,kind='output')
        job.add('--link','hard',kind='parameter')
        job.run()

        # convert the quality scores to Illumina Solexa version 1.5 format
        infile = output_file
        output_file = outdir(os.path.basename(infile).replace('init-','init-phred-'))
        job.add('phred.py',kind='program')
        job.add('--link','hard',kind='parameter')
        job.add('--input',infile,kind='input',temp_path=temp_flag)
        job.add('--output',output_file,kind='output')
        job.add('--input_type','auto-detect',kind='parameter')
        #job.add('--output_type','illumina-1.5',kind='parameter')
        job.add('--output_type','sanger',kind='parameter')
        job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
        job.run()

        new_list_input_files.append(output_file)

    job.add('printf',kind='program')
    job.add('"\n\n\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    output_file = outdir('origin.fq')
    # concatenate reads before trimming
    if len(list_input_files) > 1:
        #job.add('concatenate.py',kind='program')
        job.add('cat',kind='program')
        job.add_list('',new_list_input_files,kind='input',temp_path=temp_flag)
        if options.trimfq < 1:
            job.add('|',kind='parameter')
            job.add('seqtk',kind='program')
            job.add('trimfq',kind='parameter')
            job.add('-q',options.trimfq,kind='parameter')
            job.add('-',kind='parameter')
        job.add('>',output_file,kind='output')
        job.run()
    else:
        job.link(new_list_input_files[0], output_file, temp_path=temp_flag)

    # compute the read lengths for the input file
    job.add('lengths_reads.py',kind='program')
    job.add('--input',outdir('origin.fq'),kind='input')
    job.add('--output',outdir('log_lengths_original_reads.txt'),kind='output')
    job.add('--counts',outdir('log_counts_original_reads.txt'),kind='output')
    job.run()
    #cat snu16/reads_acgt.fq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq

    max_len_reads = 0
    if os.path.exists(outdir('log_lengths_original_reads.txt')):
        max_len_reads = int(file(outdir('log_lengths_original_reads.txt'),'r').readline().rstrip())
        if not options.bridges:
            options.bridges = int(math.ceil(float(max_len_reads)/float(160)))
        
    # save lengths reads
    info(job,
         fromfile = outdir('log_lengths_original_reads.txt'),
         tofile = info_file,
         top = ["Length of all original reads:",
                "-----------------------------"],
         bottom = "\n\n\n")

#    job.add('printf',kind='program')
#    job.add('"\nLength of Original Reads:\n-----------------\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('log_lengths_original_reads.txt'),kind='input',temp_path=temp_flag)
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
         
         
    if shuffled and (not options.skip_compress_ids):
        # lossy compression of the reads ids
        job.add('compress-reads-ids.py',kind='program')
        job.add('--input',outdir('origin.fq'),kind='input',temp_path=temp_flag)
        job.add('--output',outdir('original.fq'),kind='output')
        job.add('--count-reads',outdir('log_counts_original_reads.txt'),kind='output')
        job.add('--lowercase',kind='parameter')
        job.run()
    else:
        job.link(outdir('origin.fq'), outdir('original.fq'), temp_path=temp_flag)





    info(job,
         fromfile = outdir('log_counts_original_reads.txt'),
         tofile = info_file,
         top = ["Total counts of all original reads (after the reads marked by Illumina as bad have been removed):",
                "-------------------------------------------------------------------------------------------------"],
         bottom = "\n\n\n",
         temp_path = temp_flag)
#    job.add('printf',kind='program')
#    job.add(('"\nTotal Reads Counts (after the reads marked by Illumina as bad have been removed):\n'+
#                '---------------------------------------------------------------------------------\n"'),kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('log_counts_original_reads.txt'),kind='input',temp_path=temp_flag)
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    if options.sonication > 80 and options.sonication <= max_len_reads:
    
        job.add('printf',kind='program')
        job.add('"\n\nInput reads are broken up bioinformatically in smaller pieces due detection of very long reads!\n--------------------------------------------------------------------\n\n"',kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

        job.add('seqtk',kind='program')
        job.add('trimfq',kind='parameter')
        job.add('-q','0.25',kind='parameter')
        job.add('',outdir('original.fq'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('original-temp.fq'),kind='output')
        job.run()

        job.add('seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-1',outdir('original-temp.fq'),kind='input')
        job.add('>',outdir('or1.fq'),kind='output')
        job.run()
        
        job.add('seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-2',outdir('original-temp.fq'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('or2.fq'),kind='output')
        job.run()

        job.add('fragment_fastq.py',kind='program')
        job.add('-1',outdir('or1.fq'),kind='input',temp_path=temp_flag)
        job.add('-2',outdir('or2.fq'),kind='input',temp_path=temp_flag)
        job.add('-f',outdir('originala.fq'),kind='output')
        #job.add('-f',outdir('soni1.fq'),kind='output')
        #job.add('-r',outdir('soni2.fq'),kind='output')
        job.add('--window-size',options.trim_psl_3end_keep,kind='parameter')
        job.add('--step-size',options.trim_psl_3end_keep-length_anchor_minimum+1,kind='parameter')
        job.add('--threshold-read',options.trim_psl_3end_keep + 10,kind='parameter')
        job.add('--anchors',options.bridges,kind='parameter')
        job.add('--skip-short',options.trim_3end_keep,kind='parameter')
        job.add('--trim-n',kind='parameter')
        job.run()

#        job.add('seqtk',kind='program')
#        job.add('mergepe',kind='parameter')
#        job.add('',outdir('soni1.fq'),kind='input',temp_path=temp_flag)
#        job.add('',outdir('soni2.fq'),kind='input',temp_path=temp_flag)
#        job.add('>',outdir('originala.fq'),kind='output')
#        job.run()

        # compute the read lengths for the input file
        job.add('lengths_reads.py',kind='program')
        job.add('--input',outdir('originala.fq'),kind='input')
        job.add('--output',outdir('log_lengths_original_reads_final.txt'),kind='output')
        job.add('--counts',outdir('log_counts_original_reads_final.txt'),kind='output')
        job.run()
        
        max_len_reads = 0
        if os.path.exists(outdir('log_lengths_original_reads_final.txt')):
            max_len_reads = int(file(outdir('log_lengths_original_reads_final.txt'),'r').readline().rstrip())

        info(job,
             fromfile = outdir('log_counts_original_reads_final.txt'),
             tofile = info_file,
             top = ["Total counts of all input reads (after breaking up them bioinformatically):",
                    "-------------------------------------------------------------------------------------------------"],
             bottom = "\n\n\n",
             temp_path = temp_flag)

        # save lengths reads
        info(job,
             fromfile = outdir('log_lengths_original_reads_final.txt'),
             tofile = info_file,
             top = ["Length of all input reads (after breaking up them bioinformatically):",
                    "---------------------------------------------------------------------"],
             bottom = "\n\n\n")
    else:
        job.link(outdir('original.fq'), outdir('originala.fq'), temp_path=temp_flag)



    input_file = outdir('originala.fq')
    output_file = outdir('original-t5-t3.fq')
    if options.trim_5end > 0 or options.trim_3end > 0:
        # trim 5
#        job.add('trim_reads.py',kind='program')
#        job.add('--input',input_file,kind='input', temp_path='no')
#        job.add('--output',output_file,kind='output')
#        job.add('--trim_end','5',kind='parameter')
#        job.add('--trim_size',options.trim_5end,kind='parameter')
#        job.run()
        job.add('seqtk',kind='program')
        job.add('trimfq',kind='parameter')
        job.add('-l','1',kind='parameter')
        if options.trim_5end > 0:
            job.add('-b',options.trim_5end,kind='parameter')
        if options.trim_3end > 0:
            job.add('-e',options.trim_3end,kind='parameter')
        job.add('',input_file,kind='input', temp_path='no')
        job.add('>',output_file,kind='output')
        job.run()
    else:
        job.link(input_file, output_file, temp_path='no')

    #input_file = output_file
    #output_file = input_file[:input_file.rfind('.fq')]+'-t3.fq'
    #if options.trim_3end > 0:
    #    # trim 3
##        job.add('trim_reads.py',kind='program')
##        job.add('--input',input_file,kind='input', temp_path=temp_flag)
##        job.add('--output',output_file,kind='output')
##        job.add('--trim_end','3',kind='parameter')
##        job.add('--trim_size',options.trim_3end,kind='parameter')
##        job.run()
    #    job.add('seqtk',kind='program')
    #    job.add('trimfq',kind='parameter')
    #    job.add('-e',options.trim_5end,kind='parameter')
    #    job.add('',input_file,kind='input', temp_path=temp_flag)
    #    job.add('>',output_file,kind='output')
    #    job.run()
    #else:
    #    job.link(input_file, output_file, temp_path=temp_flag)


    input_file = output_file
    output_file = outdir('reads.fq')
    if options.trim_3end_keep > 0:
        # trim from 3-end to have the reads all the same length
        job.add('trim_reads.py',kind='program')
        job.add('--input',input_file,kind='input', temp_path=temp_flag)
        job.add('--output',output_file,kind='output')
        job.add('--trim_end','3',kind='parameter')
        job.add('--trim_n',kind='parameter')
        job.add('--final_size',options.trim_3end_keep,kind='parameter')
        job.run()
#        job.add('seqtk',kind='program')
#        job.add('trimfq',kind='parameter')
#        #job.add('-q','0',kind='parameter')
#        job.add('-l','1',kind='parameter')
#        job.add('-B',options.trim_3end_keep,kind='parameter')
#        job.add('',input_file,kind='input', temp_path=temp_flag)
#        job.add('>',output_file,kind='output')
#        job.run()
    else:
        job.link(input_file, output_file, temp_path=temp_flag)

    # compute the read lengths for the input file
    job.add('lengths_reads.py',kind='program')
    job.add('--input',outdir('reads.fq'),kind='input')
    job.add('--output',outdir('log_lengths_reads.txt'),kind='output')
    job.run()

    #job.add(kind='program')
    if os.path.exists(outdir('log_lengths_reads.txt')):
        len_reads = int(file(outdir('log_lengths_reads.txt'),'r').readline().rstrip('\r\n'))
        # reads shorter than this will be skipped from analysis, 34?
        minimum_length_short_read = len_reads
    if job.run():
        file(outdir('log_minimum_length_short_read.txt'),'w').write(str(minimum_length_short_read))


    min_len_reads = 0
    if os.path.exists(outdir('log_minimum_length_short_read.txt')):
        min_len_reads = int(file(outdir('log_minimum_length_short_read.txt'),'r').readline().rstrip())

    # save lengths reads
    info(job,
         fromfile = outdir('log_lengths_reads.txt'),
         tofile = info_file,
         top = ["Lengths of all reads after trimming:",
                "------------------------------------"],
         bottom = "\n\n\n")
#    job.add('printf',kind='program')
#    job.add('"\nLength Reads (after trimming):\n--------------------------\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('log_lengths_reads.txt'),kind='input')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    ##############################################################################
    # FILTERING - ambiguous + Bs + too short
    ##############################################################################


    if options.skip_b_filtering:
        #job.link(outdir('reads.fq'), outdir('reads_no-shorts.fq'), temp_path=temp_flag)
        job.add('seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-L',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
        job.add('',outdir('reads.fq'),kind='input',temp_path=temp_flag)
        #job.add('>',outdir('reads_no-shorts.fq'),kind='output')
        if (not options.all_reads_junction) and (not options.skip_interleave_processing):
            job.add('|',kind='parameter')
            job.add('seqtk',kind='parameter')
            job.add('dropse',kind='parameter')
            job.add('-',kind='parameter')
        job.add('>',outdir('reads_acgt.fq'),kind='output')
        job.run()
    else:
        job.add('seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-L',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
        job.add('',outdir('reads.fq'),kind='input',temp_path=temp_flag)
        #job.add('>',outdir('reads_no-shorts.fq'),kind='output')
        if (not options.all_reads_junction) and (not options.skip_interleave_processing):
            job.add('|',kind='parameter')
            job.add('seqtk',kind='parameter')
            job.add('dropse',kind='parameter')
            job.add('-',kind='parameter')
        job.add('|',kind='parameter')
        # fix Illumina "B"
        job.add('fastq_b2n.py',kind='parameter')
        #job.add('--input',input_file,kind='input',temp_path=temp_flag)
        job.add('--input','-',kind='parameter')
        job.add('--replacement','A',kind='parameter')
        job.add('--sanger',kind='parameter')
        job.add('--ambiguous',kind='parameter')
        job.add('--threshold',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file = 'yes')
        #job.add('--output',output_file,kind='output')
        #job.add('--output',outdir('reads_b2n2a.fq'),kind='output')
        job.add('--output','-',kind='output')
        #job.run()
        # filter out the reads with poly tail
        # trim the poly tails
        #job.add('trim_poly_tails.py',kind='program')
        job.add('|',kind='parameter')
        job.add('trim_poly_tails.py',kind='parameter')
        #job.add('--input',outdir('reads_b2n2a.fq'),kind='input',temp_path=temp_flag)
        job.add('--input','-',kind='parameter')
        job.add('--repeats','20',kind='parameter') # 12
        #job.add('--skip_reads',kind='parameter')
        job.add('--output','-',kind='output')
        job.add('2>>',outdir('info.txt'),kind='output')
        job.add('|',kind='parameter')
        job.add('seqtk',kind='parameter')
        job.add('seq',kind='parameter')
        job.add('-L',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
        #job.add('',outdir('reads_b2n2a.fq'),kind='input',temp_path=temp_flag)
        job.add('-',kind='parameter')
        #job.add('>',outdir('reads_no-shorts.fq'),kind='output')
        if (not options.all_reads_junction) and (not options.skip_interleave_processing):
            job.add('|',kind='parameter')
            job.add('seqtk',kind='parameter')
            job.add('dropse',kind='parameter')
            job.add('-',kind='parameter')
        job.add('>',outdir('reads_acgt.fq'),kind='output')
        job.run()


    #convert ambiguous nucleotides to As
    #job.add('fastq2acgt.py',kind='program')
    #job.add('--input',outdir('reads_b2n.fq'),kind='input',temp_path=temp_flag)
    #job.add('--output',outdir('reads_b2n2a.fq'),kind='output')
    #job.run()

    if pigz:
        job.add('pigz',kind='program')
        job.add('-p',options.processes,kind='parameter',checksum='no')
    else:
        job.add('gzip',kind='program')
    job.add('--fast',kind='parameter')
    job.add('-c',outdir('originala.fq'),kind='input',temp_path=temp_flag)
    job.add('>',outdir('originala.fq.gz'),kind='output')
    job.run()





    # remove reads shorter than a given threshold
#    job.add('remove_shorter_reads.py',kind='program')
#    job.add('--input',outdir('reads_b2n2a.fq'),kind='input',temp_path=temp_flag)
#    job.add('--threshold',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
#    job.add('--output',outdir('reads_no-shorts.fq'),kind='output')
#    job.run()
#    job.add('seqtk',kind='program')
#    job.add('seq',kind='parameter')
#    job.add('-L',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
#    job.add('',outdir('reads_b2n2a.fq'),kind='input',temp_path=temp_flag)
#    job.add('>',outdir('reads_no-shorts.fq'),kind='output')
#    job.run()

    # remove the reads which do not form a pair
##    if (not options.all_reads_junction) and (not options.skip_interleave_processing):
        # assumption all reads are interleaved
#        job.add('remove_single_reads.py',kind='program')
#        job.add('--input',outdir('reads_no-shorts.fq'),kind='input',temp_path=temp_flag)
#        job.add('--interleaved',kind='parameter')
#        job.add('--output',outdir('reads_acgt.fq'),kind='output')
#        job.add('--log',outdir('log_reads2_removed.txt'),kind='output')
#        job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
#        job.add('--processes',options.processes,kind='parameter',checksum='no')
#        job.add('2>',outdir('log_removed_single_reads1.txt'),kind='parameter',checksum='no')
#        job.run()

#        job.add('cat',kind='program')
#        job.add('',outdir('reads_no-shorts.fq'),kind='input',temp_path=temp_flag)
#        job.add('|',kind='parameter')
#        job.add('paste','- - - -',kind='parameter')
#        job.add('|',kind='parameter')
#        job.add('awk',kind='parameter')
#        job.add('-F"\\\\t"',kind='parameter')
#        job.add("""'{ n=length($1); if (olde=="/1" && substr($1,0,n-1)==old && substr($1,n-1,2)=="/2") {print old1"\\n"old2"\\n+\\n"old3"\\n"$1"\\n"$2"\\n+\\n"$4; old=""; count++} {old=substr($1,0,n-1); olde=substr($1,n-1,2); old1=$1; old2=$2; old3=$4}} END{print NR-2*count > "%s"}'""" % (outdir('log_removed_single_reads1.txt'),),kind='parameter')
#        job.add('',outdir('log_removed_single_reads1.txt'),kind='output',command_line='no')
#        job.add('>',outdir('reads_acgt.fq'),kind='output')
#        job.run()
##        job.add('seqtk',kind='program')
##        job.add('dropse',kind='parameter')
##        job.add('',outdir('reads_no-shorts.fq'),kind='input',temp_path=temp_flag)
##        job.add('>',outdir('reads_acgt.fq'),kind='output')
##        job.run()

#        job.add('printf',kind='program')
#        job.add(('"\n\nCount of short reads removed due to missing their mate read:\n'+
#                 '-----------------------------------------------------------------\n"'),kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
#        job.add('cat',kind='program')
#        job.add('',outdir('log_removed_single_reads1.txt'),kind='input',temp_path=temp_flag)
#        job.add('>>',info_file,kind='output')
#        job.run()
#        job.add('printf',kind='program')
#        job.add('"\n\n\n"',kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
##    else:
##        job.link(outdir('reads_no-shorts.fq'),outdir('reads_acgt.fq'),temp_path=temp_flag)


    job.add('cat',kind='program')
    job.add('',outdir('reads_acgt.fq'),kind='input')
    job.add('|',kind='parameter')
    job.add("echo $((`wc -l`/4))",kind='parameter')
    job.add('>>',outdir('log_removed_single_reads1.txt'),kind='output')
    job.run()

    info(job,
        fromfile = outdir('log_removed_single_reads1.txt'),
        tofile = info_file,
        top = ["\n\nCount of all short reads after removing reads due to missing their mate read:",
               "-----------------------------------------------------------------------------"],
        bottom = "\n\n\n")


    if job.iff(empty(outdir('reads_acgt.fq')),id = "#reads_acgt.fq#"):
        t = ["ERROR: Too many reads have been removed during the pre-filtering steps!",
             "Please, check that the input files are from a RNA-seq dataset with pair-reads ",
             "or that the input files are given correctly!"
             "Please, check that also the input reads have the same length!"
             ]
        t = '\n'.join(t)+'\n'
        print >>sys.stderr, t
        file(info_file,'a').write(t)
        file(log_file,'a').write(t)
        job.clean(outdir('original.fq'),temp_path=temp_flag)
        job.clean(outdir('original.fq.gz'),temp_path=temp_flag)
        job.clean(outdir('originala.fq'),temp_path=temp_flag)
        job.clean(outdir('originala.fq.gz'),temp_path=temp_flag)
        job.clean(outdir('log_lengths_original_reads.txt'),temp_path=temp_flag)
        job.clean(outdir('log_lengths_original_reads_plus.txt'),temp_path=temp_flag)
        job.clean(outdir('log_lengths_original_reads_final.txt'),temp_path=temp_flag)
        job.clean(outdir('log_lengths_reads.txt'),temp_path=temp_flag)
        job.clean(outdir('log_removed_single_reads1.txt'),temp_path=temp_flag)
        job.clean(outdir('log_minimum_length_short_read.txt'),temp_path=temp_flag)
        job.clean(outdir('reads_acgt.fq'),temp_path=temp_flag)
        job.close()
        sys.exit(1)

    no_reads = 0
    if os.path.isfile(outdir('log_removed_single_reads1.txt')):
        no_reads = int(file(outdir('log_removed_single_reads1.txt'),'r').readline().strip())

####
####
    if not options.skip_automatic_scaling:
        if max_len_reads < 60:
            options.skip_bowtie2 = True
            options.skip_bwa = True
            options.skip_star_bowtie = True
            if not is_optparse_provided(parser,'limit_star'):
                options.limit_star = 3 * (2**30)
            job.add('printf',kind='program')
            job.add(('"\nInput reads are too short (maxim found length is %d) in order to use BWA and therefore BOWTIE2 method is disabled automatically!\n"') % (max_len_reads,),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()
        if (not is_optparse_provided(parser,'mismatches_psl')) and (not is_optparse_provided(parser,'trim_psl_3end_keep')):
            if max_len_reads > 109:
                options.trim_psl_3end_keep = max_len_reads - 20
                options.mismatches_psl = math.ceil(float(options.trim_psl_3end_keep)/50)
                job.add('printf',kind='program')
                job.add('"\nAdjusted automatically mismatches_psl and trim_psl_3end_keep (%d,%d) because reads of maximum length of %d bp were found!.\n\n"' % (options.trim_psl_3end_keep,options.mismatches_psl,max_len_reads),kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()


        if (no_reads and
            (not is_optparse_provided(parser,'spanning_pairs')) and
            (not is_optparse_provided(parser,'spanning_reads')) and
            (not is_optparse_provided(parser,'length_anchor')) and
            (not options.sensitive) and
            (not options.mildly_sensitive) and
            (not options.highly_sensitive) and
            (not options.paranoid_sensitive)):
            if no_reads < 35000000 and no_reads >= 15000000 and max_len_reads < 60:
                spanning_pairs_bowtie = 3
                spanning_pairs_minimum = min([spanning_pairs_bowtie,spanning_pairs_blat,spanning_pairs_star,spanning_pairs_bowtie2,spanning_pairs_bwa])
                spanning_reads_bowtie = 2
                spanning_reads_minimum = min([spanning_reads_bowtie,spanning_reads_blat,spanning_reads_star,spanning_reads_bowtie2,spanning_reads_bwa])
                length_anchor_bowtie = 14
                length_anchor_minimum = min([length_anchor_bowtie,length_anchor_bowtie2,length_anchor_star,length_anchor_blat,length_anchor_bwa])
                job.add('printf',kind='program')
                job.add('"\nAdjusted automatically spanning_pairs, spanning_reads, and length_anchor (%d,%d,%d,%d) for count of reads [15000000,3500000) and reads shorter than 60bp.\n\n"' % (spanning_pairs_bowtie,spanning_reads_bowtie,length_anchor_bowtie,length_anchor_bwa),kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
            if no_reads < 15000000 and no_reads >= 1000000 and max_len_reads < 60:
                spanning_pairs_bowtie = 2
                spanning_pairs_minimum = min([spanning_pairs_bowtie,spanning_pairs_blat,spanning_pairs_star,spanning_pairs_bowtie2,spanning_pairs_bwa])
                spanning_reads_bowtie = 2
                spanning_reads_minimum = min([spanning_reads_bowtie,spanning_reads_blat,spanning_reads_star,spanning_reads_bowtie2,spanning_reads_bwa])
                length_anchor_bowtie = 13
                length_anchor_minimum = min([length_anchor_bowtie,length_anchor_bowtie2,length_anchor_star,length_anchor_blat,length_anchor_bwa])
                job.add('printf',kind='program')
                job.add('"\nAdjusted automatically spanning_pairs, spanning_reads, and length_anchor (%d,%d,%d,%d) for count of reads [1000000,15000000) and reads shorter than 60bp.\n\n"' % (spanning_pairs_bowtie,spanning_reads_bowtie,length_anchor_bowtie,length_anchor_bwa),kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
            if no_reads < 20000000 and no_reads >= 1000000 and max_len_reads > 74:
                spanning_pairs_bowtie = 2
                spanning_pairs_minimum = min([spanning_pairs_bowtie,spanning_pairs_blat,spanning_pairs_star,spanning_pairs_bowtie2,spanning_pairs_bwa])
                spanning_reads_bowtie = 2
                spanning_reads_minimum = min([spanning_reads_bowtie,spanning_reads_blat,spanning_reads_star,spanning_reads_bowtie2,spanning_reads_bwa])
                job.add('printf',kind='program')
                job.add('"\nAdjusted automatically spanning_pairs, spanning_reads, and length_anchor (%d,%d,%d,%d) for count of reads [1000000,20000000) and reads longer than 74 bp.\n\n"' % (spanning_pairs_bowtie,spanning_reads_bowtie,length_anchor_bowtie,length_anchor_bwa),kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
####
####

    if job.run():
        if (2 * length_anchor_minimum) > len_reads - 1:
            job.write(["ERROR: The length of the anchor (i.e. %d) is too long compared to the length of the reads (i.e. %d )" % (options.length_anchor,len_reads),
                       "Suggestion 1: Decrease the length of the anchor using '--anchor-fusion' option if possible!",
                       "Suggestion 2: Decrease the size of the trimming using '--5end' or '--5keep' option if possible!",
                      ], stdout = True, stderr = True, log = True)
            job.close()
            sys.exit(1)


#    if job.iff(not empty(outdir('log_reads2_removed.txt')),id = "#log_reads2_removed.txt#"):
#        r = float(file(outdir('log_reads2_removed.txt'),'r').readline())
#        if r > 0.7:
#            t = ["ERROR: Too many reads (that is %.3f%%) have been removed because they miss theirs read-mates!" % (r,),
#                 "Please, check that the input files are from a RNA-seq dataset with pair-reads or that the input files are given correctly!"
#                 ]
#            t = '\n'.join(t)+'\n'
#            print >>sys.stderr, t
#            file(info_file,'a').write(t)
#            file(log_file,'a').write(t)
#            job.close()
#            sys.exit(1)
#    job.clean(outdir('log_reads2_removed.txt'),temp_path=temp_flag)


    ##############################################################################
    # FILTERING - ribosomal DNA + mitochondrion
    ##############################################################################

    # find available memory
    job.add('printf',kind='program')
    job.add('"\n============\nMEMORY (before using BOWTIE):\n============\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('free',kind='program')
    job.add('-m',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()
    job.add('printf',kind='program')
    job.add('"\n\n\n"',kind='parameter')
    job.add('>>',info_file,kind='output')
    job.run()

    # map using the filter index (not aligned, unique alignment, multiple alignments)
    job.add('bowtie',kind='program')
    job.add('-t',kind='parameter')
    #job.add('-q',kind='parameter')
    #job.add('-a',kind='parameter')
    job.add('-v',options.filter_mismatches,kind='parameter') #options.mismatches
    job.add('-p',options.processes,kind='parameter',checksum='no')
    #job.add('-m','1',kind='parameter')
    job.add('-k','1',kind='parameter')
    #job.add('--solexa1.3-quals',kind='parameter')
    job.add('--phred33-quals',kind='parameter')
    job.add('--tryhard',kind='parameter')
    job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
    #job.add('--best',kind='parameter')
    #job.add('--strata',kind='parameter')
    if len_reads > 40 and options.trim_wiggle:
        job.add('--trim3',options.trim_wiggle,kind='parameter') # trim on the fly 5bp from 3' end
        job.add('--trim5',options.trim_wiggle,kind='parameter') # trim the 5
#    job.add('--suppress','1,2,3,4,5,6,7,8',kind='parameter')
    job.add('--un',outdir('reads-filtered_temp.fq'),kind='output') # here is the result
    job.add('--max',outdir('reads-filtered_temp_multiple.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
    if options.skip_mitochondrion_filtering:
        if os.path.isfile(datadir('rtrna_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('rtrna_index/'),kind='input')
    else:
        if os.path.isfile(datadir('rtrna_mt_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('rtrna_mt_index/'),kind='input')
    job.add('',outdir('reads_acgt.fq'),kind='input',temp_path=temp_flag)
    job.add('',outdir('reads-filtered.map'),kind='output',temp_path=temp_flag)
#    job.add('','/dev/null',kind='parameter')
    job.add('2>',outdir('log_bowtie_reads-filtered-out.stdout.txt'),kind='output',checksum='no')
    #job.add('2>&1',kind='parameter',checksum='no')
    job.run()


    info(job,
        fromfile = outdir('log_bowtie_reads-filtered-out.stdout.txt'),
        tofile = info_file,
        top = ["Mapping all input reads on rRNA and/or MT for filtering purposes:",
               "----------------------------------------------------------------"],
        bottom = "\n\n\n",
        temp_path = temp_flag)

#    job.add('printf',kind='program')
#    job.add(('"\n\nMapping all input reads on rRNA and/or MT for filtering purposes:\n'+
#             '-----------------------------------------------\n"'),kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('log_bowtie_reads-filtered-out.stdout.txt'),kind='input',temp_path=temp_flag)
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    # remove the reads which map on PHIX174 which is used for improving the quality of Illumina NGS
    if not options.skip_phix174_filtering:

        if job.iff(empty(outdir('reads-filtered_temp.fq')),id = "#reads-filtered_temp.fq#"):
            t = ["ERROR: Too many reads have been removed during the pre-filtering steps!",
                 "Please, check that the input files are from a RNA-seq dataset with pair-reads "
                 "or that the input files are given correctly!"
                 ]
            t = '\n'.join(t)+'\n'
            print >>sys.stderr, t
            file(info_file,'a').write(t)
            file(log_file,'a').write(t)
            job.close()
            sys.exit(1)

        # map using the filter index (not aligned, unique alignment, multiple alignments)
        job.add('bowtie',kind='program')
        job.add('-t',kind='parameter')
        #job.add('-q',kind='parameter')
        #job.add('-a',kind='parameter')
        job.add('-v','2',kind='parameter') #options.mismatches
        job.add('-p',options.processes,kind='parameter',checksum='no')
        #job.add('-m','1',kind='parameter')
        job.add('-k','1',kind='parameter')
        #job.add('--solexa1.3-quals',kind='parameter')
        job.add('--phred33-quals',kind='parameter')
        job.add('--tryhard',kind='parameter')
        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
        #job.add('--best',kind='parameter')
        #job.add('--strata',kind='parameter')
        if len_reads > 40 and options.trim_wiggle:
            job.add('--trim3',options.trim_wiggle,kind='parameter') # trim on the fly 5bp from 3' end
            job.add('--trim5',options.trim_wiggle,kind='parameter') # trim the 5
#        job.add('--suppress','1,2,3,4,5,6,7,8',kind='parameter')
        job.add('--un',outdir('reads-filtered_temp-phix.fq'),kind='output') # here is the result
        job.add('--max',outdir('reads-filtered_temp_multiple-phix.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
        if os.path.isfile(datadir('phix174_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('phix174_index/'),kind='input')
        job.add('',outdir('reads-filtered_temp.fq'),kind='input',temp_path=temp_flag)
        job.add('',outdir('reads-filtered-phix.map'),kind='output',temp_path=temp_flag)
        job.add('2>',outdir('log_bowtie_reads-filtered-out-phix.stdout.txt'),kind='parameter',checksum='no')
        #job.add('2>&1',kind='parameter',checksum='no')
        job.run()

        info(job,
            fromfile = outdir('log_bowtie_reads-filtered-out-phix.stdout.txt'),
            tofile = info_file,
            top = ["Mapping all input reads on Enterobacteria phage phiX174 genome (used for improving quality of Illumina NGS) for filtering purposes:",
                   "-----------------------------------------------------------------------------------------------------------------------------------"],
            bottom = "\n\n\n",
            temp_path = temp_flag)

#        job.add('printf',kind='program')
#        job.add(('"\n\nMapping all input reads on Enterobacteria phage phiX174 genome (used for improving quality of Illumina NGS) for filtering purposes:\n'+
#                 '----------------------------------------------------------------------------------------------------------------------------------------\n"'),kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
#        job.add('cat',kind='program')
#        job.add('',outdir('log_bowtie_reads-filtered-out-phix.stdout.txt'),kind='input',temp_path=temp_flag)
#        job.add('>>',info_file,kind='output')
#        job.run()
#        job.add('printf',kind='program')
#        job.add('"\n\n\n"',kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
    else:
        job.link(outdir('reads-filtered_temp.fq'),outdir('reads-filtered_temp-phix.fq'),temp_path=temp_flag)


    # remove the reads which map on HLA
    if not options.skip_hla_filtering:
    #if options.hla_filtering:

        if job.iff(empty(outdir('reads-filtered_temp-phix.fq')),id = "#reads-filtered_temp-phix.fq#"):
            t = ["ERROR: Too many reads have been removed during the pre-filtering steps!",
                 "Please, check that the input files are from a RNA-seq dataset with pair-reads ",
                 "or that the input files are given correctly!"
                 ]
            t = '\n'.join(t)+'\n'
            print >>sys.stderr, t
            file(info_file,'a').write(t)
            file(log_file,'a').write(t)
            job.close()
            sys.exit(1)

        job.add('printf',kind='program')
        job.add(('"Mapping all input reads on HLA human histocompatibility complexes:\n'+
                 '-------------------------------------------------------------------\n"'),kind='parameter')
        job.add('>>',info_file,kind='output')
        job.run()

        # map using the filter index (not aligned, unique alignment, multiple alignments)
        job.add('bowtie',kind='program')
        job.add('-t',kind='parameter')
        #job.add('-q',kind='parameter')
        #job.add('-a',kind='parameter')
        job.add('-v',options.filter_mismatches,kind='parameter') #options.mismatches
        job.add('-p',options.processes,kind='parameter',checksum='no')
        #job.add('-m','1',kind='parameter')
        job.add('-k','1',kind='parameter')
        #job.add('--solexa1.3-quals',kind='parameter')
        job.add('--phred33-quals',kind='parameter')
        job.add('--tryhard',kind='parameter')
        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
        #job.add('--best',kind='parameter')
        #job.add('--strata',kind='parameter')
        if len_reads > 40 and options.trim_wiggle:
            job.add('--trim3',options.trim_wiggle,kind='parameter') # trim on the fly 5bp from 3' end
            job.add('--trim5',options.trim_wiggle,kind='parameter') # trim the 5
        job.add('--suppress','1,2,4,5,6,7,8',kind='parameter')
        job.add('--un',outdir('reads-filtered_temp-hla.fq'),kind='output') # here is the result
        job.add('--max',outdir('reads-filtered_temp_multiple-hla.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
        if os.path.isfile(datadir('hla_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('hla_index/'),kind='input')
        job.add('',outdir('reads-filtered_temp-phix.fq'),kind='input',temp_path=temp_flag)
#        job.add('',outdir('reads-filtered-hla.map'),kind='output',temp_path=temp_flag)
#        job.add('2>',outdir('log_bowtie_reads-filtered-out-hla.stdout.txt'),kind='parameter',checksum='no')
        job.add('2>>',info_file,kind='parameter',checksum='no')
        #job.add('2>&1',kind='parameter',checksum='no')
#        job.run()
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('-c',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        job.add('-rn',kind='parameter')
        job.add('-k','1,1',kind='parameter')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('cat',kind='parameter')
        job.add('|',kind='parameter')
        job.add('head','-20',kind='parameter')
        job.add('>',outdir('hla_statistics.txt'),kind='output')
        job.run()


        info(job,
             fromfile = outdir('hla_statistics.txt'),
             tofile = info_file,
             top = ["","","",
                    "HLA types (histocompatibility complexes) found to have reads mapped on their sequences:",
                    "---------------------------------------------------------------------------------------",
                    "Reads_count\tHLA type"],
             bottom = "\n\n\n",
             temp_path=temp_flag)

    else:
        job.link(outdir('reads-filtered_temp-phix.fq'),outdir('reads-filtered_temp-hla.fq'),temp_path=temp_flag)

    # remove the reads which do not form a pair
    if (not options.all_reads_junction) and (not options.skip_interleave_processing):
        # assumption all reads are not interleaved
        job.add('cat',kind='program')
        job.add('',outdir('reads-filtered_temp-hla.fq'),kind='input',temp_path=temp_flag)
        job.add('|',kind='parameter')
        job.add('paste','- - - -',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        job.add('-k','1,1',kind='parameter')
        job.add('-t',"'\t'",kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('tr',kind='parameter')
        job.add('"\\t"',kind='parameter')
        job.add('"\\n"',kind='parameter')
        job.add('|',kind='parameter')
        job.add('seqtk',kind='parameter')
        job.add('dropse',kind='parameter')
        #job.add('-',kind='parameter')
        job.add('>',outdir('reads-filtered.fq'),kind='output')
        job.run()

        job.add('cat',kind='program')
        job.add('',outdir('reads-filtered.fq'),kind='input')
        job.add('|',kind='parameter')
        job.add("echo $((`wc -l`/4))",kind='parameter')
        job.add('>>',outdir('log_removed_single_reads2.txt'),kind='output')
        job.run()

#        job.add('awk',kind='parameter')
#        job.add('-F"\\\\t"',kind='parameter')
#        job.add("""'{n=length($1); if (olde=="/1" && substr($1,0,n-1)==old && substr($1,n-1,2)=="/2") {print old1"\\n"old2"\\n+\\n"old3"\\n"$1"\\n"$2"\\n+\\n"$4; old=""; count++} {old=substr($1,0,n-1); olde=substr($1,n-1,2); old1=$1; old2=$2; old3=$4}} END {print NR-2*count > "%s"}'""" % (outdir('log_removed_single_reads2.txt'),),kind='parameter')
#        job.add('',outdir('log_removed_single_reads2.txt'),kind='output',command_line='no')
#        job.add('>',outdir('reads-filtered.fq'),kind='output')
#        job.run()

#        job.add('tr',kind='parameter')
#        job.add('"\\t"',kind='parameter')
#        job.add('"\\n"',kind='parameter')
#        job.add('|',kind='parameter')
#        job.add('remove_single_reads.py',kind='parameter')
#        ##job.add('--input',outdir('reads-filtered_temp-hla.fq'),kind='input',temp_path=temp_flag)
#        job.add('--interleaved',kind='parameter')
#        job.add('--input','-',kind='input')
#        job.add('--output',outdir('reads-filtered.fq'),kind='output')
#        job.add('--log',outdir('log_reads1_removed.txt'),kind='output')
#        job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
#        job.add('--processes',options.processes,kind='parameter',checksum='no')
#        job.add('2>',outdir('log_removed_single_reads2.txt'),kind='parameter',checksum='no')
#        job.run()


        info(job,
             fromfile = outdir('log_removed_single_reads2.txt'),
             tofile = info_file,
             top = ["\n\nCount of all short reads after removing reads due to missing their mate read:",
                    "-----------------------------------------------------------------------------"],
             bottom = "\n\n\n",
             temp_path = temp_flag)

#        job.add('printf',kind='program')
#        job.add(('"\n\nCount of short reads removed due to missing their mate read:\n'+
#                 '-----------------------------------------------------------------\n"'),kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
#        job.add('cat',kind='program')
#        job.add('',outdir('log_removed_single_reads2.txt'),kind='input',temp_path=temp_flag)
#        job.add('>>',info_file,kind='output')
#        job.run()
#        job.add('printf',kind='program')
#        job.add('"\n\n\n"',kind='parameter')
#        job.add('>>',info_file,kind='output')
#        job.run()
    else:
        job.link(outdir('reads-filtered_temp-hla.fq'),outdir('reads-filtered.fq'),temp_path=temp_flag)

    if job.run():
        if not empty(outdir('log_reads1_removed.txt')):
            r = float(file(outdir('log_reads1_removed.txt'),'r').readline())
            if r > 0.7:
                t = ["ERROR: Too many reads (that is %.3f%%) have been removed because they miss theirs read-mates!" % (r,),
                     "Please, check that the input files are from a RNA-seq dataset with pair-reads or that the input files are given correctly!"
                     ]
                t = '\n'.join(t)+'\n'
                print >>sys.stderr, t
                file(info_file,'a').write(t)
                file(log_file,'a').write(t)
                job.close()
                sys.exit(1)
    job.clean(outdir('log_reads1_removed.txt'),temp_path=temp_flag)

    job.add('cat',kind='program')
    job.add('',outdir('reads-filtered.fq'),kind='input')
    job.add('|',kind='parameter')
    job.add("echo $((`wc -l`/4))",kind='parameter')
    job.add('>',outdir('count_reads_left_after_filtering.txt'),kind='output')
    job.run()


    info(job,
         fromfile = outdir('count_reads_left_after_filtering.txt'),
         tofile = info_file,
         top = ["Total Reads Counts (after the all filtering steps):",
                "---------------------------------------------------"],
         bottom = "\n\n\n")




#    job.add('printf',kind='program')
#    job.add(('"\nTotal Reads Counts (after the all filtering steps):\n'+
#                '--------------------------------------------------\n"'),kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('count_reads_left_after_filtering.txt'),kind='input')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    if job.iff(empty(outdir('reads-filtered.fq')),id = "#reads-filtered.fq#"):
        t = ["ERROR: Too many reads have been removed during the pre-filtering steps!",
             "Please, check that the input files are from a RNA-seq dataset with pair-reads "
             "or that the input files are given correctly!"
             ]
        t = '\n'.join(t)+'\n'
        print >>sys.stderr, t
        file(info_file,'a').write(t)
        file(log_file,'a').write(t)
        job.close()
        sys.exit(1)


    ##############################################################################
    # MAPPING short reads against the genome
    ##############################################################################
    # map using the genome index (not aligned, unique alignment, multiple alignments); results in MAP BOWTIE format
    job.add('bowtie',kind='program')
    job.add('-t',kind='parameter')
    #job.add('-q',kind='parameter')
    job.add('-a',kind='parameter')
    #job.add('-k','2',kind='parameter')
    job.add('-v','1',kind='parameter') # options.filter_mismatches # stjude
    job.add('-p',options.processes,kind='parameter',checksum='no')
    job.add('-m','1',kind='parameter')
    #job.add('-k','1',kind='parameter')
    job.add('--suppress','5,6,7',kind='parameter')
    #job.add('--solexa1.3-quals',kind='parameter')
    job.add('--phred33-quals',kind='parameter')
    job.add('--best',kind='parameter')
    job.add('--strata',kind='parameter')
    job.add('--tryhard',kind='parameter') # ??? really necessary? stjude
    job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
    job.add('--un',outdir('reads_filtered_not-mapped-genome.fq'),kind='output')
    job.add('--max',outdir('reads-filtered_multiple-mappings-genome.fq'),kind='output') # if this is missing then these reads are going to '--un'
    if os.path.isfile(datadir('genome_index','.1.ebwtl')):
        job.add('--large-index',kind='parameter')
    job.add('',datadir('genome_index/'),kind='input')
    job.add('',outdir('reads-filtered.fq'),kind='input')
    job.add('',outdir('reads_filtered_genome.map'),kind='output') # <== best mappings on genome #######
    job.add('2>',outdir('log_bowtie_reads_mapped-genome.stdout.txt'),kind='parameter',checksum='no')
    #job.add('2>&1',kind='parameter',checksum='no')
    job.run()

    info(job,
         fromfile = outdir('log_bowtie_reads_mapped-genome.stdout.txt'),
         tofile = info_file,
         top = ["Mapping the filtered reads on genome:",
                "-------------------------------------"],
         bottom = "\n\n\n")



#    job.add('printf',kind='program')
#    job.add(('"\n\nMapping the filtered reads on genome:\n'+
#             '-----------------------------------------------\n"'),kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('log_bowtie_reads_mapped-genome.stdout.txt'),kind='input')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    # extract the names of the short reads which mapped on the genome
    job.add('cut',kind='program')
    job.add('-f','1',kind='parameter')
    job.add('',outdir('reads_filtered_genome.map'),kind='input')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    if sort_buffer:
        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
    if sort_parallel:
        job.add('--parallel',options.processes,kind='parameter',checksum='no')
    if sort_lzop_compress:
        job.add('--compress-program','lzop',kind='parameter',checksum='no')
    elif sort_gzip_compress:
        job.add('--compress-program','gzip',kind='parameter',checksum='no')
    job.add('-T',tmp_dir,kind='parameter',checksum='no')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('list-names-reads-filtered_genome.txt'),kind='output')
    job.run()



    if not options.split_seqtk_subseq:
        #extract the short reads which mapped on genome
        job.add('extract_short_reads.py',kind='program')
        job.add('--input',outdir('reads-filtered.fq'),kind='input')
        job.add('--list',outdir('list-names-reads-filtered_genome.txt'),kind='input')
        job.add('--output',outdir('reads_filtered_unique-mapped-genome.fq'),kind='output')
        job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
        job.run(error_message = ("If this fails (again?) due to a memory error (e.g. not enough free memory) then lowering the "+
                                 "buffer size for specifically this script might help. This can be done by using the FusionCatcher's "+
                                 "command line option '--extra-buffer-size "+str(int(options.extract_buffer_size)/2)+"' ."))
    elif options.split_seqtk_subseq == 1:
        #extract the short reads which mapped on genome
        job.add('seqtk',kind='program')
        job.add('subseq',kind='parameter')
        job.add('',outdir('reads-filtered.fq'),kind='input')
        job.add('',outdir('list-names-reads-filtered_genome.txt'),kind='input')
        job.add('>',outdir('reads_filtered_unique-mapped-genome.fq'),kind='output')
        job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
            "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))
    elif options.split_seqtk_subseq > 1:
        #extract the short reads which mapped on genome
        job.add('seqtk-subseq.sh',kind='program')
        job.add('',options.split_seqtk_subseq,kind='parameter')
        job.add('',outdir('reads-filtered.fq'),kind='input')
        job.add('',outdir('list-names-reads-filtered_genome.txt'),kind='input')
        job.add('',outdir('reads_filtered_unique-mapped-genome.fq'),kind='output')
        job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
            "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))

    job.clean(outdir('reads-filtered.fq'),temp_path=temp_flag)
    job.clean(outdir('list-names-reads-filtered_genome.txt'),temp_path=temp_flag)

        #IDEA
        # cat list-names-reads-filtered_genome.txt | gnu_parallel --part -k -j1 --block 10G seqtk subseq reads-filtered.fq - > reads_filtered_unique-mapped-genome.fq


    # sed 's/^/^@/g' ids.txt > aha.txt ; cat /apps/reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq | paste - - - - | /apps/tools/parallel-20140822/src/parallel --pipe --no-notice grep -f aha.txt - | tr "\t" "\n" > result.fq

    ##############################################################################
    # MAPPING short reads which do not map on genome against the transcriptome
    ##############################################################################

    if job.iff(empty(outdir('reads_filtered_not-mapped-genome.fq')),id = "#reads_filtered_not-mapped-genome.fq#"):
        job.add('echo',kind='program')
        job.add('-n','""',kind='parameter')
        job.add('>',outdir('reads_filtered_not-mapped-genome_transcriptome.map'),kind='output')
        job.run()

        job.add('printf',kind='program')
        job.add('"\nMapping on transcriptome the filtered reads which did not map on genome:\n------------------------------------------------------------------------\nNo reads and no alignments!\n\n\n"', kind='parameter')
        job.add('>>',log_file,kind='output')
        job.run()

    else:
        # map using the transcript index (not mapped, unique alignment, multiple alignments)
        job.add('bowtie',kind='program')
        job.add('-t',kind='parameter')
        #job.add('-q',kind='parameter')
        job.add('-a',kind='parameter')
        job.add('-v',options.mismatches,kind='parameter')
        job.add('-p',options.processes,kind='parameter',checksum='no')
        #job.add('--solexa1.3-quals',kind='parameter')
        job.add('--phred33-quals',kind='parameter')
        job.add('--suppress','5,6,7',kind='parameter')
        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
        job.add('--tryhard',kind='parameter')
        job.add('--best',kind='parameter')
        job.add('--strata',kind='parameter')
        job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.fq'),kind='output') # <== reads which do not map on transcriptome and genome! #######
        job.add('--max',outdir('reads_filtered_genome-transcriptome_multiple.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
        if os.path.isfile(datadir('transcripts_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('transcripts_index/'),kind='input')
        job.add('',outdir('reads_filtered_not-mapped-genome.fq'),kind='input')
        job.add('',outdir('reads_filtered_not-mapped-genome_transcriptome.map'),kind='output')
        job.add('2>',outdir('log_bowtie_reads_not-mapped-genome_but_mapped-transcriptome.stdout.txt'),kind='parameter',checksum='no')
        #job.add('2>&1',kind='parameter',checksum='no')
        job.run()

        info(job,
             fromfile = outdir('log_bowtie_reads_not-mapped-genome_but_mapped-transcriptome.stdout.txt'),
             tofile = info_file,
             top = ["Mapping on transcriptome the filtered reads which did not map on genome:",
                    "------------------------------------------------------------------------"],
             bottom = "\n\n\n")

#    job.add('printf',kind='program')
#    job.add(('"\n\nMapping on transcriptome the filtered reads which did not map on genome:\n'+
#                  '------------------------------------------------------------------------\n"'),kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('log_bowtie_reads_not-mapped-genome_but_mapped-transcriptome.stdout.txt'),kind='input')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('printf',kind='program')
#    job.add('"\n\n\n"',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    # extract ids of short reads which mapped on the transcriptome
    job.add('cut',kind='program')
    job.add('-f','1',kind='parameter')
    job.add('',outdir('reads_filtered_not-mapped-genome_transcriptome.map'),kind='input')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    if sort_buffer:
        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
    if sort_parallel:
        job.add('--parallel',options.processes,kind='parameter',checksum='no')
    if sort_lzop_compress:
        job.add('--compress-program','lzop',kind='parameter',checksum='no')
    elif sort_gzip_compress:
        job.add('--compress-program','gzip',kind='parameter',checksum='no')
    job.add('-T',tmp_dir,kind='parameter',checksum='no')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('list-names-reads-filtered_not-mapped-genome_mapped-transcriptome.txt'),kind='output')
    job.run()


    if not options.split_seqtk_subseq:
        #extract the short reads which mapped on the transcriptome and do not map on genome
        job.add('extract_short_reads.py',kind='program')
        job.add('--input',outdir('reads_filtered_not-mapped-genome.fq'),kind='input')
        job.add('--list',outdir('list-names-reads-filtered_not-mapped-genome_mapped-transcriptome.txt'),kind='input')
        job.add('--output',outdir('reads_filtered_not-mapped-genome_mapped-transcriptome.fq'),kind='output')
        job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
        job.run(error_message = ("If this fails (again?) due to a memory error (e.g. not enough free memory) then lowering the "+
                                 "buffer size for specifically this script might help. This can be done by using the FusionCatcher's "+
                                 "command line option '--extra-buffer-size "+str(int(options.extract_buffer_size)/2)+"' ."))
    elif options.split_seqtk_subseq == 1:
        #extract the short reads which mapped on the transcriptome and do not map on genome
        job.add('seqtk',kind='program')
        job.add('subseq',kind='parameter')
        job.add('',outdir('reads_filtered_not-mapped-genome.fq'),kind='input')
        job.add('',outdir('list-names-reads-filtered_not-mapped-genome_mapped-transcriptome.txt'),kind='input')
        job.add('>',outdir('reads_filtered_not-mapped-genome_mapped-transcriptome.fq'),kind='output')
        job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
            "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))
    elif options.split_seqtk_subseq > 1:
        #extract the short reads which mapped on the transcriptome and do not map on genome
        job.add('seqtk-subseq.sh',kind='program')
        job.add('',options.split_seqtk_subseq,kind='parameter')
        job.add('',outdir('reads_filtered_not-mapped-genome.fq'),kind='input')
        job.add('',outdir('list-names-reads-filtered_not-mapped-genome_mapped-transcriptome.txt'),kind='input')
        job.add('',outdir('reads_filtered_not-mapped-genome_mapped-transcriptome.fq'),kind='output')
        job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
            "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))

    job.clean(outdir('reads_filtered_not-mapped-genome.fq'),temp_path=temp_flag)
    job.clean(outdir('list-names-reads-filtered_not-mapped-genome_mapped-transcriptome.txt'),temp_path=temp_flag)

    ##############################################################################
    # MAPPING short reads (which map uniquely on genome) against the transcriptome
    ##############################################################################

    # map using the transcript index (not mapped, unique alignment, multiple alignments)
    job.add('bowtie',kind='program')
    job.add('-t',kind='parameter')
    #job.add('-q',kind='parameter')
    job.add('-a',kind='parameter')
    job.add('-v',options.mismatches,kind='parameter')
    job.add('-p',options.processes,kind='parameter',checksum='no')
    job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
    #job.add('--solexa1.3-quals',kind='parameter')
    job.add('--phred33-quals',kind='parameter')
    job.add('--best',kind='parameter')
    job.add('--tryhard',kind='parameter')
    job.add('--strata',kind='parameter')
    job.add('--suppress','5,6,7',kind='parameter')
    if os.path.isfile(datadir('transcripts_index','.1.ebwtl')):
        job.add('--large-index',kind='parameter')
    job.add('',datadir('transcripts_index/'),kind='input')
    job.add('',outdir('reads_filtered_unique-mapped-genome.fq'),kind='input')
    job.add('',outdir('reads_filtered_unique-mapped-genome_transcriptome_temp.map'),kind='output')
    job.add('2>',outdir('log_bowtie_reads_unique-mapped-genome_mapped-transcriptome.stdout.txt'),kind='parameter',checksum='no')
    #job.add('2>&1',kind='parameter',checksum='no')
    job.run()

    info(job,
         fromfile = outdir('log_bowtie_reads_unique-mapped-genome_mapped-transcriptome.stdout.txt'),
         tofile = info_file,
         top = ["Mapping on transcriptome the filtered reads which map uniquely on genome:",
                "------------------------------------------------------------------------"],
         bottom = "\n\n\n")

    # filter the mapped reads on transcriptome wich mapped also on genome using mismatches
    job.add('remove_reads_genome_transcriptome.py',kind='program')
    job.add('--input_map_1',outdir('reads_filtered_genome.map'),kind='input',temp_path=temp_flag)
    job.add('--input_map_2',outdir('reads_filtered_unique-mapped-genome_transcriptome_temp.map'),kind='input',temp_path=temp_flag)
    job.add('--mismatches_column','5',kind='parameter')
    job.add('--output',outdir('reads_filtered_unique-mapped-genome_transcriptome.map'),kind='output')
    job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
    job.run()

    # extract the names of the short reads which mapped on the transcriptome
    job.add('cut',kind='program')
    job.add('-f','1',kind='parameter')
    job.add('',outdir('reads_filtered_unique-mapped-genome_transcriptome.map'),kind='input')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    if sort_buffer:
        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
    if sort_parallel:
        job.add('--parallel',options.processes,kind='parameter',checksum='no')
    if sort_lzop_compress:
        job.add('--compress-program','lzop',kind='parameter',checksum='no')
    elif sort_gzip_compress:
        job.add('--compress-program','gzip',kind='parameter',checksum='no')
    job.add('-T',tmp_dir,kind='parameter',checksum='no')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('list-names-reads-filtered_unique-mapped-genome_mapped-transcriptome.txt'),kind='output')
    job.run()

#    job.add('printf',kind='program')
#    job.add('"Count reads left mapping transcriptome and genome (before filtering out those with better mappings on genome):\n---------------------------------------------------\n "',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('reads_filtered_unique-mapped-genome.fq'),kind='input')
#    job.add('|',kind='parameter')
#    job.add("echo $((`wc -l`/4))",kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    if not options.split_seqtk_subseq:
        #extract the short reads which mapped on the transcriptome and do not map on genome
        job.add('extract_short_reads.py',kind='program')
        job.add('--input',outdir('reads_filtered_unique-mapped-genome.fq'),kind='input')
        job.add('--list',outdir('list-names-reads-filtered_unique-mapped-genome_mapped-transcriptome.txt'),kind='input')
        job.add('--output',outdir('reads_filtered_unique-mapped-genome_mapped-transcriptome.fq'),kind='output')
        job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
        job.run(error_message = ("If this fails (again?) due to a memory error (e.g. not enough free memory) then lowering the "+
                                 "buffer size for specifically this script might help. This can be done by using the FusionCatcher's "+
                                 "command line option '--extra-buffer-size "+str(int(options.extract_buffer_size)/2)+"' ."))
    elif options.split_seqtk_subseq == 1:
        #extract the short reads which mapped on the transcriptome and do not map on genome
        job.add('seqtk',kind='program')
        job.add('subseq',kind='parameter')
        job.add('',outdir('reads_filtered_unique-mapped-genome.fq'),kind='input')
        job.add('',outdir('list-names-reads-filtered_unique-mapped-genome_mapped-transcriptome.txt'),kind='input')
        job.add('>',outdir('reads_filtered_unique-mapped-genome_mapped-transcriptome.fq'),kind='output')
        job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
            "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))
    elif options.split_seqtk_subseq > 1:
        #extract the short reads which mapped on the transcriptome and do not map on genome
        job.add('seqtk-subseq.sh',kind='program')
        job.add('',options.split_seqtk_subseq,kind='parameter')
        job.add('',outdir('reads_filtered_unique-mapped-genome.fq'),kind='input')
        job.add('',outdir('list-names-reads-filtered_unique-mapped-genome_mapped-transcriptome.txt'),kind='input')
        job.add('',outdir('reads_filtered_unique-mapped-genome_mapped-transcriptome.fq'),kind='output')
        job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
            "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))
        
    job.clean(outdir('reads_filtered_unique-mapped-genome.fq'),temp_path=temp_flag)
    job.clean(outdir('list-names-reads-filtered_unique-mapped-genome_mapped-transcriptome.txt'),temp_path=temp_flag)
        
#    job.add('printf',kind='program')
#    job.add('"Count reads left mapping on transcriptome and genome after filtering out those with better mappings on genome:\n---------------------------------------------------\n "',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('reads_filtered_unique-mapped-genome_mapped-transcriptome.fq'),kind='input')
#    job.add('|',kind='parameter')
#    job.add("echo $((`wc -l`/4))",kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    # group reads which map on transcriptome in one FASTQ file
    #job.add('concatenate.py',kind='program')
    job.add('cat',kind='program')
    job.add('',outdir('reads_filtered_not-mapped-genome_mapped-transcriptome.fq'),kind='input',temp_path=temp_flag)
    job.add('',outdir('reads_filtered_unique-mapped-genome_mapped-transcriptome.fq'),kind='input',temp_path=temp_flag)
    job.add('>',outdir('reads_filtered_mapped-transcriptome.fq'),kind='output')
    job.run()


#    job.add('printf',kind='program')
#    job.add('"Total count reads mapping on transcriptome (mapping on genome and not mapping on genome):\n---------------------------------------------------\n "',kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()
#    job.add('cat',kind='program')
#    job.add('',outdir('reads_filtered_mapped-transcriptome.fq'),kind='input')
#    job.add('|',kind='parameter')
#    job.add("echo $((`wc -l`/4))",kind='parameter')
#    job.add('>>',info_file,kind='output')
#    job.run()

    # group reads' mappings on transcriptome in one MAP file
    #job.add('concatenate.py',kind='program')
    job.add('cat',kind='program')
    job.add('',outdir('reads_filtered_not-mapped-genome_transcriptome.map'),kind='input',temp_path=temp_flag)
    job.add('',outdir('reads_filtered_unique-mapped-genome_transcriptome.map'),kind='input',temp_path=temp_flag)
#    job.add('-',kind='parameter')  # <== best mappings on transcriptome #######
    #job.add('',outdir('reads_filtered_transcriptome.map'),kind='output')  # <== best mappings on transcriptome #######
    #job.run()
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    if sort_buffer:
        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
    if sort_parallel:
        job.add('--parallel',options.processes,kind='parameter',checksum='no')
    if sort_lzop_compress:
        job.add('--compress-program','lzop',kind='parameter',checksum='no')
    elif sort_gzip_compress:
        job.add('--compress-program','gzip',kind='parameter',checksum='no')
    job.add('-T',tmp_dir,kind='parameter',checksum='no')
    #job.add('-s',kind='parameter') # stable sort
    job.add('-t',"'\t'",kind='parameter')
    job.add('-k','1,1',kind='parameter')
    #job.add('',outdir('reads_filtered_transcriptome.map'),kind='input',temp_path = temp_flag)
    job.add('>',outdir('reads_filtered_transcriptome_sorted-read.map'),kind='output')
#    job.add('|',kind='parameter')
#    if pigz:
#        job.add('pigz',kind='parameter')
#        job.add('-p',options.processes,kind='parameter',checksum='no')
#    else:
#        job.add('gzip',kind='parameter')
#    job.add('>',outdir('reads_filtered_transcriptome_sorted-read.map.gz'),kind='output')
    job.run()

    if ( (options.homolog > 0 or (not options.skip_ambiguous_filtering)) and
         job.iff(not empty(outdir('reads_filtered_mapped-transcriptome.fq')),id = "#reads_filtered_mapped-transcriptome.fq#")
         ):
        ##############################################################################
        # ALL POSSIBLE MAPPINGS of all short reads on transcriptome
        ##############################################################################
        # map against the transcriptome the short reads which do not map on the genome
        # map using the transcript index (not mapped, unique alignment, multiple alignments)
        classic_1 = True
        job.add('bowtie',kind='program')
        job.add('-t',kind='parameter')
        #job.add('-q',kind='parameter')
        #job.add('-a',kind='parameter')
        job.add('-k','200',kind='parameter')
        #job.add('-v',options.mismatches,kind='parameter')
        job.add('-v', options.ambiguous_mismatches,kind='parameter') # stjude
        ##job.add('--strata',kind='parameter') # stjude ??????
        job.add('-p',options.processes,kind='parameter',checksum='no')
        #job.add('--solexa1.3-quals',kind='parameter')
        job.add('--phred33-quals',kind='parameter')
        #job.add('--best',kind='parameter')
        if classic_1:
            job.add('--suppress','2,4,5,6,7,8',kind='parameter') # original
        else:
            job.add('--suppress','2,4,5,6,7',kind='parameter') # stjude
        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
        if os.path.isfile(datadir('transcripts_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('transcripts_index/'),kind='input')
        job.add('',outdir('reads_filtered_mapped-transcriptome.fq'),kind='input',temp_path = temp_flag)
        #job.add('',outdir('reads_filtered_all-possible-mappings-transcriptome.map'),kind='output') # <== best mappings on transcriptome ####### XXX
        job.add('2>',outdir('log_bowtie_reads_filtered_all-possible-mappings-transcriptome_map.stdout.txt'),kind='output',checksum='no')
        #job.add('2>&1',kind='parameter',checksum='no')
        #job.run() #XXX
        #job.add('|',kind='parameter') # XXX
        # sort the reads' all possible mappings on transcriptome by reads name
        #job.add('LC_ALL=C',kind='parameter')
        job.add('|',kind='parameter') # XXX
        if classic_1:
            job.add('sed',kind='parameter')
            job.add("""'s/\\tEN.*\;EN/\\tEN/'""",kind='parameter')
        else:
            job.add('awk',kind='parameter')
            job.add("""'{print $1"\\t"substr($2,index($2,";")+1)"\\t"gsub(">","",$3)}'""",kind='parameter')
        job.add('|',kind='parameter') # XXX
        job.add('uniq',kind='parameter') # XXX
        #job.add('|',kind='parameter') # XXX
        job.add('>',outdir('reads_filtered_all-possible-mappings-transcriptome.map'),kind='output') # <== best mappings on transcriptome ####### XXX
        job.run()
        #
        job.add('LC_ALL=C',kind='program') # XXX
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        #job.add('-s',kind='parameter') # stable sort
        job.add('-t',"'\t'",kind='parameter')
        job.add('-k','1,1',kind='parameter')
        job.add('',outdir('reads_filtered_all-possible-mappings-transcriptome.map'),kind='input',temp_path=temp_flag) # XXX
        #job.add('>',outdir('reads_filtered_all-possible-mappings-transcriptome_sorted.map'),kind='output')
        #job.run()
#        job.add('|',kind='parameter') # XXX
#        job.add('grep',kind='parameter') # remove the ENSG09 genes
#        job.add('-v',kind='parameter')
#        job.add('-F',kind='parameter')
#        job.add('-f',datadir('custom_genes_mark.txt'),kind='input')
        job.add('|',kind='parameter') # XXX
        job.add('uniq',kind='parameter') # XXX
        #job.add('|',kind='parameter') # XXX
        job.add('>',outdir('reads_filtered_all-possible-mappings-transcriptome_sorted.map'),kind='output')
        job.run() #XXX
        # find the homolog genes using the reads
        job.add('find_homolog_genes.py',kind='program')
        job.add('--input',outdir('reads_filtered_all-possible-mappings-transcriptome_sorted.map'),kind='input',temp_path=temp_flag) # XXX
        #job.add('--input','-',kind='parameter')
        #job.add('--reads',outdir('log_number_of_reads_processed.txt'),kind='parameter',from_file='yes')
        job.add('--input_exons',datadir('exons.txt'),kind='input')
        job.add('--filter',datadir('custom_genes_mark.txt'),kind='input')
        job.add('--processes',options.processes,kind='parameter')
        job.add('--reads','1',kind='parameter')
        if not classic_1:
            job.add('--d1',kind='parameter') # stjude -- only 1 mismatch away
        #job.add('--output_offending_reads',outdir('list_offending_reads.txt'),kind='output')
        job.add('--output_offending_pair_reads',outdir('list_offending_reads_.txt'),kind='output')
        if options.homolog > 0:
            job.add('--output',outdir('list_candidates_ambiguous_homologous_genes_1.txt'),kind='output') # <== list of genes that might be homologous #######
        else:
            job.add('--output',outdir('list_candidates_ambiguous_homologous_genes_1.txt'),kind='output',temp_path=temp_flag) # <== list of genes that might be homologous #######
        job.run()

        job.add('LC_ALL=C',kind='program')
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('',outdir('list_offending_reads_.txt'),kind='input',temp_path=temp_flag)
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('>',outdir('list_offending_reads.txt'),kind='output')
        job.run()

        info(job,
             fromfile = outdir('log_bowtie_reads_filtered_all-possible-mappings-transcriptome_map.stdout.txt'),
             tofile = info_file,
             top = ["All mappings of reads (not mapped on genome + mapped uniquely on genome + mapped on transcriptome) on transcriptome",
                    "-------------------------------------------------------------------------------------------------------------------"],
             bottom = "\n\n\n",
             temp_path = temp_flag
             )


        # get the number of reads mapping on transcriptome and genome
        crgt_ = 0
        if job.run():
            tag = '# reads processed:'
            crgt_ = [line.rstrip('\r\n') for line in file(outdir('log_bowtie_reads_unique-mapped-genome_mapped-transcriptome.stdout.txt'),'r') if line.lower().find(tag)!=-1]
            crgt_ = float(crgt_.pop(0).split(tag)[1].strip())
            crgt = int(crgt_ * options.homolog)
            file(outdir('log_number_of_reads_processed.txt'),'w').write(str(crgt))

        info(job,
            fromfile = outdir('log_number_of_reads_processed.txt'),
            tofile = info_file,
            top = ["Threshold for pairs of genes to be marked as 'similar_reads' (reference number: %d)" % (int(crgt_),),
                   "-----------------------------------------------------------------------------------------"],
            bottom = "\n\n\n")

        ##############################################################################
        # ALL MAPPINGS of short reads (which map multiple times on genome) against the transcriptome
        ##############################################################################
        if job.iff(not empty(outdir('reads-filtered_multiple-mappings-genome.fq')),id="#reads-filtered_multiple-mappings-genome.fq#"):
            classic_2 = True
            job.add('bowtie',kind='program')
            job.add('-p',options.processes,kind='parameter',checksum='no')
            job.add('-t',kind='parameter')
            #job.add('-q',kind='parameter')
            #job.add('-a',kind='parameter')
            job.add('-k','200',kind='parameter')
            job.add('-v',options.ambiguous_mismatches,kind='parameter')
            #job.add('--strata',kind='parameter') # stjude
            #job.add('--solexa1.3-quals',kind='parameter')
            job.add('--phred33-quals',kind='parameter')
            #job.add('--best',kind='parameter')
            if classic_2:
                job.add('--suppress','2,4,5,6,7,8',kind='parameter')
            else:
                job.add('--suppress','2,4,5,6,7',kind='parameter')
            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
            if os.path.isfile(datadir('transcripts_index','.1.ebwtl')):
                job.add('--large-index',kind='parameter')
            job.add('',datadir('transcripts_index/'),kind='input')
            job.add('',outdir('reads-filtered_multiple-mappings-genome.fq'),kind='input',temp_path=temp_flag)
            #job.add('',outdir('reads_filtered_all-possible-mappings-transcriptome_multiple.map'),kind='output') # <== best mappings on transcriptome #######
            job.add('2>',outdir('log_bowtie_reads_filtered_all-possible-mappings-transcriptome_multiple_map.stdout.txt'),kind='output',checksum='no')
            # job.add('2>&1',kind='parameter',checksum='no') #XXX
            #job.run()
            job.add('|',kind='parameter') # XXX
            # sort the reads' all possible mappings on transcriptome by reads name
            #job.add('awk',kind='program')
            #job.add("""awk '{print $1"\\t"substr($2,index($2,";")+1)"\\t"gsub(">","",$3)}'""",outdir('reads_filtered_all-possible-mappings-transcriptome.map'),kind='input',temp_path=temp_flag)
            if classic_2:
                job.add('sed',kind='parameter')
                job.add("""'s/\\tEN.*\;EN/\\tEN/'""",kind='parameter')
            else:
                job.add('awk',kind='parameter')
                job.add("""'{print $1"\\t"substr($2,index($2,";")+1)"\\t"gsub(">","",$3)}'""",kind='parameter')
            job.add('|',kind='parameter') # XXX
            job.add('uniq',kind='parameter') # XXX
            #job.add('|',kind='parameter') # XXX
            job.add('>',outdir('reads_filtered_all-possible-mappings-transcriptome_multiple.map'),kind='output') # <== best mappings on transcriptome #######
            job.run()

            job.add('LC_ALL=C',kind='program') # XXX
            #job.add('LC_ALL=C',kind='parameter')
            job.add('sort',kind='parameter')
            if sort_buffer:
                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
            if sort_parallel:
                job.add('--parallel',options.processes,kind='parameter',checksum='no')
            if sort_lzop_compress:
                job.add('--compress-program','lzop',kind='parameter',checksum='no')
            elif sort_gzip_compress:
                job.add('--compress-program','gzip',kind='parameter',checksum='no')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            job.add('-t',"'\t'",kind='parameter')
            job.add('-k','1,1',kind='parameter')
            job.add('',outdir('reads_filtered_all-possible-mappings-transcriptome_multiple.map'),kind='input',temp_path=temp_flag) # XXX
            #job.add('>',outdir('reads_filtered_all-possible-mappings-transcriptome_multiple_sorted.map'),kind='output') # XXX
            #job.run() # XXX
            job.add('|',kind='parameter') # XXX
            job.add('uniq',kind='parameter') # XXX
            job.add('>',outdir('reads_filtered_all-possible-mappings-transcriptome_multiple_sorted.map'),kind='output') # XXX
            job.run() #XXX
#            job.add('grep',kind='parameter') # remove the ENSG09 genes
#            job.add('-v',kind='parameter')
#            job.add('-F',kind='parameter')
#            job.add('-f',datadir('custom_genes_mark.txt'),kind='input')
#            job.add('|',kind='parameter') # XXX
            # find the homolog genes using the reads
            job.add('find_homolog_genes.py',kind='program')
            job.add('--input',outdir('reads_filtered_all-possible-mappings-transcriptome_multiple_sorted.map'),kind='input',temp_path=temp_flag)
            #job.add('--input','-',kind='parameter')
            job.add('--reads','1',kind='parameter')
            if not classic_2:
                job.add('--d1',kind='parameter') # stjude -- only 1 mismatch away # using this requires bowtie '--suppress','4,5,6,7' instead of '--suppress','4,5,6,7,8'
            job.add('--input_exons',datadir('exons.txt'),kind='input')
            job.add('--filter',datadir('custom_genes_mark.txt'),kind='input')
            job.add('--processes',options.processes,kind='parameter')
            job.add('--output',outdir('list_candidates_ambiguous_homologous_genes_2.txt'),kind='output') # <== list of genes that might be homologous #######
            job.run()

            info(job,
                fromfile = outdir('log_bowtie_reads_filtered_all-possible-mappings-transcriptome_multiple_map.stdout.txt'),
                tofile = info_file,
                top = ["Mapping all short reads (which already are mapping multiple times on genome) on transcriptome:",
                       "----------------------------------------------------------------------------------------------"],
                bottom = "\n\n\n",
                temp_path = temp_flag)

        else:
            job.add('echo',kind='program')
            job.add('-n','""',kind='parameter')
            job.add('>',outdir('list_candidates_ambiguous_homologous_genes_2.txt'),kind='output')
            job.run()


        # join the found homolog genes using the reads
        job.add('join_homolog_genes.py',kind='program')
        job.add('--input_1',outdir('list_candidates_ambiguous_homologous_genes_1.txt'),kind='input',temp_path=temp_flag)
        job.add('--input_2',outdir('list_candidates_ambiguous_homologous_genes_2.txt'),kind='input',temp_path=temp_flag)
        job.add('--reads',outdir('log_number_of_reads_processed.txt'),kind='parameter',from_file='yes')
        job.add('--all',outdir('all_ambiguous_genes.txt'),kind='output')
        if options.homolog > 0:
            job.add('--output',outdir('list_candidates_ambiguous_homologous_genes.txt'),kind='output') # <== list of genes that might be homologous #######
        else:
            job.add('--output',outdir('list_candidates_ambiguous_homologous_genes.txt'),kind='output',temp_path=temp_flag) # <== list of genes that might be homologous #######
        job.run()


        if not options.skip_ambiguous_filtering:
            # remove the offending reads from the transcriptome mapping
#            job.add('reads_from_map.py',kind='program')
#            job.add('--input_reads',outdir('list_offending_reads.txt'),kind='input',temp_path=temp_flag)
#            job.add('--input_map',outdir('reads_filtered_transcriptome_sorted-read.map'),kind='input')
#            job.add('--output_map',outdir('reads_filtered_transcriptome_sorted-read_no-offending-reads.map'),kind='output')
#            job.add('--operation','remove',kind='parameter')
#            job.run()
            job.add('LC_ALL=C',kind='program')
            job.add('join',kind='parameter')
            job.add('-1','1',kind='parameter')
            job.add('-2','1',kind='parameter')
            job.add('-v','2',kind='parameter')
            job.add('-t',"'\t'",kind='parameter')
            job.add('',outdir('list_offending_reads.txt'),kind='input',temp_path=temp_flag)
            job.add('',outdir('reads_filtered_transcriptome_sorted-read.map'),kind='input')
            job.add('>',outdir('reads_filtered_transcriptome_sorted-read_no-offending-reads.map'),kind='output')
            job.run()
        else:
            job.link(outdir('reads_filtered_transcriptome_sorted-read.map'), outdir('reads_filtered_transcriptome_sorted-read_no-offending-reads.map'),temp_path='no')
            job.clean(outdir('list_offending_reads.txt'),temp_path=temp_flag)
    else:
        job.link(outdir('reads_filtered_transcriptome_sorted-read.map'), outdir('reads_filtered_transcriptome_sorted-read_no-offending-reads.map'),temp_path='no')


    ##############################################################################
    # FIND FUSION GENES
    ##############################################################################

    # find ALL fusion genes and transcripts where the offending reads have been removed
    # (offending reads = reads which map at least on two different genes)
    job.add('find_fusion_genes_map.py',kind='program')
    job.add('--input',outdir('reads_filtered_transcriptome_sorted-read_no-offending-reads.map'),kind='input',temp_path=temp_flag)
    job.add('--input_hugo',datadir('genes_symbols.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_no-offending-reads.txt'),kind='output')
    job.add('--output_fusion_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='output')
    if options.reads_preliminary_fusions:
        job.add('--output_fusion_reads_split',outdir('pre-fusion'),kind='output')
    #job.add('--output_fusion_reads_simple',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads_only-ids.txt'),kind='output')
    job.add('--output_missing_mate_reads',outdir('candidate_fusion-genes_missing_mates.txt'),kind='output')
    job.run()

    if options.reads_preliminary_fusions:
        job.add('concatenate.py',kind='program')
        job.add('-f',outdir('pre-fusion'),kind='input')
        job.add('-',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('>',outdir('pre-fusion_ids.txt'),kind='output')
        job.run()

        job.add('seqtk',kind='program')
        job.add('subseq',kind='parameter')
        job.add('',outdir('originala.fq.gz'),kind='input')
        job.add('',outdir('pre-fusion_ids.txt'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('originala-pre-fusion.fq'),kind='output')
        job.run()

        parts = [el.strip() for el in file(outdir('pre-fusion'),'r').readlines()]
        for par in parts:
            job.add('seqtk',kind='program')
            job.add('subseq',kind='parameter')
            job.add('',outdir('originala-pre-fusion.fq'),kind='input')
            job.add('',par,kind='input',temp_path=temp_flag)
            job.add('>',par+'.fq',kind='output')
            job.run()
            
        job.clean(outdir('originala-pre-fusion.fq'),temp_path=temp_flag)
        
    # label fusion genes -- no protein product
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_no-offending-reads.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','no_protein',kind='parameter')
    job.add('--similar_gene_symbols',kind='parameter')
    job.add('--filter_genes',datadir('genes_with_no_proteins.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_01.txt'),kind='output')
    job.run()
    # label fusion genes -- paralogs
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_01.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','paralogs',kind='parameter')
    job.add('--filter_gene_pairs',datadir('paralogs.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_02.txt'),kind='output')
    job.run()
    # label fusion genes -- potential readthrough
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_02.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','adjacent',kind='parameter')
    job.add('--filter_gene_pairs',datadir('adjacent_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_03.txt'),kind='output')
    job.run()
    # label fusion genes -- fully overlapping in Ensembl
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_03.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ensembl_fully_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ensembl_fully_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_04.txt'),kind='output')
    job.run()
    # label fusion genes -- partially overlapping in Ensembl
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_04.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ensembl_partially_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ensembl_partially_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_05.txt'),kind='output')
    job.run()
    # label fusion genes -- overlapping and on same strand in Ensembl
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_05.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ensembl_same_strand_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ensembl_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_06.txt'),kind='output')
    job.run()
    # label fusion genes -- similar region
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_06.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','similar_reads',kind='parameter')
    job.add('--filter_gene_pairs',outdir('list_candidates_ambiguous_homologous_genes.txt'),kind='input',temp_path=temp_flag)
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_07.txt'),kind='output')
    job.run()
    # label fusion genes -- minimum distance between genes on the same strand
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_07.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','short_distance',kind='parameter')
    job.add('--min_dist_gene_gene',options.min_dist,kind='parameter')
    job.add('--min_dist_gene_gene_database',datadir('exons.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_08.txt'),kind='output')
    job.run()
    # label fusion genes -- minimum distance between genes on the same strand
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_08.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','distance1000bp',kind='parameter')
    job.add('--min_dist_gene_gene','1000',kind='parameter')
    job.add('--min_dist_gene_gene_database',datadir('exons.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_09.txt'),kind='output')
    job.run()
    # label fusion genes -- minimum distance between genes on the same strand
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_09.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','distance10kbp',kind='parameter')
    job.add('--min_dist_gene_gene','10000',kind='parameter')
    job.add('--min_dist_gene_gene_database',datadir('exons.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_10.txt'),kind='output')
    job.run()
    # label fusion genes -- minimum distance between genes on the same strand
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_10.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','distance100kbp',kind='parameter')
    job.add('--min_dist_gene_gene','100000',kind='parameter')
    job.add('--min_dist_gene_gene_database',datadir('exons.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_11.txt'),kind='output')
    job.run()
    # label fusion genes -- pseudogenes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_11.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','pseudogene',kind='parameter')
    job.add('--filter_genes',datadir('pseudogenes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_12.txt'),kind='output')
    job.run()
    # label fusion genes -- rRNA (again)
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_12.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','rrna',kind='parameter')
    job.add('--filter_genes',datadir('rrnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_13.txt'),kind='output')
    job.run()
    # label fusion genes -- tRNA
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_13.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','trna',kind='parameter')
    job.add('--filter_genes',datadir('trnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_14.txt'),kind='output')
    job.run()
    # label fusion genes -- miRNA
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_14.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','mirna',kind='parameter')
    job.add('--filter_genes',datadir('mirnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_15.txt'),kind='output')
    job.run()
    # label fusion genes -- lincRNA
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_15.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','lincrna',kind='parameter')
    job.add('--filter_genes',datadir('lincrnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_16.txt'),kind='output')
    job.run()
    # label fusion genes -- MT
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_16.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','mt',kind='parameter')
    job.add('--filter_genes',datadir('mt.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_17.txt'),kind='output')
    job.run()
    # label fusion genes -- snoRNA
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_17.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','snorna',kind='parameter')
    job.add('--filter_genes',datadir('snornas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_18.txt'),kind='output')
    job.run()
    # label fusion genes -- snRNA
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_18.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','snrna',kind='parameter')
    job.add('--filter_genes',datadir('snrnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_19.txt'),kind='output')
    job.run()
    # label fusion genes -- Y RNAs
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_19.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','yrna',kind='parameter')
    job.add('--filter_genes',datadir('yrnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_20.txt'),kind='output')
    job.run()
    # label fusion genes -- 7SK RNAs
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_20.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','7skrna',kind='parameter')
    job.add('--filter_genes',datadir('7skrnas.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_21.txt'),kind='output')
    job.run()
    # label fusion genes -- antisense
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_21.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','antisense',kind='parameter')
    job.add('--filter_genes',datadir('antisenses.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_22.txt'),kind='output')
    job.run()
    # label fusion genes -- paralogs
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_22.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','pair_pseudo_genes',kind='parameter')
    job.add('--filter_gene_pairs',datadir('pairs_pseudogenes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_23.txt'),kind='output')
    job.run()
    # label fusion genes -- ribosomal proteins
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_23.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ribosomal_protein',kind='parameter')
    job.add('--filter_genes',datadir('ribosomal_proteins.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_24.txt'),kind='output')
    job.run()
    # label fusion genes -- oncogenes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_24.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','oncogene',kind='parameter')
    job.add('--filter_genes',datadir('oncogenes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_25.txt'),kind='output')
    job.run()
    # label fusion genes -- known fusions
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_25.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','known_fusion',kind='parameter')
    job.add('--filter_gene_pairs',datadir('known.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_26.txt'),kind='output')
    job.run()
    # label fusion genes -- cosmic
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_26.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','cosmic',kind='parameter')
    job.add('--filter_gene_pairs',datadir('cosmic.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_27.txt'),kind='output')
    job.run()
    # label fusion genes -- ChimerDB 2.0
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_27.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','chimerdb2',kind='parameter')
    job.add('--filter_gene_pairs',datadir('chimerdb2.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_28.txt'),kind='output')
    job.run()
    # label fusion genes -- CGP
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_28.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','cgp',kind='parameter')
    job.add('--filter_gene_pairs',datadir('cgp.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_29.txt'),kind='output')
    job.run()
    # label fusion genes -- ConjoinG
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_29.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','conjoing',kind='parameter')
    job.add('--filter_gene_pairs',datadir('conjoing.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_30.txt'),kind='output')
    job.run()
    # label fusion genes -- TICdb
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_30.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ticdb',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ticdb.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_31.txt'),kind='output')
    job.run()
    # label fusion genes -- RP11-... genes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_31.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','rp11_gene',kind='parameter')
    job.add('--filter_genes',datadir('rp11.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_32.txt'),kind='output')
    job.run()
    # label fusion genes -- CTA-... genes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_32.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','cta_gene',kind='parameter')
    job.add('--filter_genes',datadir('cta.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_33.txt'),kind='output')
    job.run()
    # label fusion genes -- CTB-... genes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_33.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ctb_gene',kind='parameter')
    job.add('--filter_genes',datadir('ctb.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_34.txt'),kind='output')
    job.run()
    # label fusion genes -- CTD-... genes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_34.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ctd_gene',kind='parameter')
    job.add('--filter_genes',datadir('ctd.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_35.txt'),kind='output')
    job.run()
    # label fusion genes -- CTC-... genes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_35.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ctc_gene',kind='parameter')
    job.add('--filter_genes',datadir('ctc.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_36.txt'),kind='output')
    job.run()
    # label fusion genes -- RP??-... genes
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_36.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','rp_gene',kind='parameter')
    job.add('--filter_genes',datadir('rp.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_37.txt'),kind='output')
    job.run()
    # label fusion genes -- found in healthy samples
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_37.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','healthy',kind='parameter')
    job.add('--filter_gene_pairs',datadir('healthy.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_38.txt'),kind='output')
    job.run()
    # label fusion genes -- CACG
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_38.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','cacg',kind='parameter')
    job.add('--filter_gene_pairs',datadir('cacg.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_39.txt'),kind='output')
    job.run()
    # label fusion genes -- fully overlapping in UCSC
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_39.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ucsc_fully_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ucsc_fully_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_40.txt'),kind='output')
    job.run()
    # label fusion genes -- partially overlapping in UCSC
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_40.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ucsc_partially_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ucsc_partially_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_41.txt'),kind='output')
    job.run()
    # label fusion genes -- overlapping and on same strand in UCSC
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_41.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ucsc_same_strand_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('ucsc_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_42.txt'),kind='output')
    job.run()
    # label fusion genes -- fully overlapping in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_42.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','refseq_fully_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('refseq_fully_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_43.txt'),kind='output')
    job.run()
    # label fusion genes -- partially overlapping in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_43.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','refseq_partially_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('refseq_partially_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_44.txt'),kind='output')
    job.run()
    # label fusion genes -- overlapping and on same strand in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_44.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','refseq_same_strand_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('refseq_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_45.txt'),kind='output')
    job.run()
    # label fusion genes -- duplicated genes from DGD database
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_45.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','duplicates',kind='parameter')
    job.add('--filter_gene_pairs',datadir('dgd.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_46.txt'),kind='output')
    job.run()
    # label fusion genes -- TCGA
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_46.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','tcga',kind='parameter')
    job.add('--filter_gene_pairs',datadir('tcga.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_47.txt'),kind='output')
    job.run()
    # label fusion genes -- BodyMap2
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_47.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','bodymap2',kind='parameter')
    job.add('--filter_gene_pairs',datadir('bodymap2.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_48.txt'),kind='output')
    job.run()
    # label fusion genes -- Metazoa
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_48.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','metazoa',kind='parameter')
    job.add('--filter_genes',datadir('metazoa.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_49.txt'),kind='output')
    job.run()
    # label fusion genes -- cell lines
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_49.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','cell_lines',kind='parameter')
    job.add('--filter_gene_pairs',datadir('celllines.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_50.txt'),kind='output')
    job.run()
    # label fusion genes -- ambiguous (only if the abguous counts > supporting pairs)
    job.add('label_ambiguous_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_50.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','ambiguous',kind='parameter')
    job.add('--factor','20',kind='parameter') # 15
    job.add('--input_ambiguous',outdir('all_ambiguous_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_51.txt'),kind='output')
    job.run()
    # label fusion genes -- fully overlapping in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_51.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','gencode_fully_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('gencode_fully_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_52.txt'),kind='output')
    job.run()
    # label fusion genes -- partially overlapping in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_52.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','gencode_partially_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('gencode_partially_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_53.txt'),kind='output')
    job.run()
    # label fusion genes -- overlapping and on same strand in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_53.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','gencode_same_strand_overlapping',kind='parameter')
    job.add('--filter_gene_pairs',datadir('gencode_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_54.txt'),kind='output')
    job.run()
    # label fusion genes -- overlapping and on same strand in RefSeq
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_54.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','prostates',kind='parameter')
    job.add('--filter_gene_pairs',datadir('prostates.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_55.txt'),kind='output')
    job.run()
    # label fusion genes -- non-tumor cell lines
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_55.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','non-tumor_cell-lines',kind='parameter')
    job.add('--filter_gene_pairs',datadir('non-tumor_cells.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_56.txt'),kind='output')
    job.run()
    # label with focus the fusions which are given by the user
    if options.focus_fusions:
        job.add('label_fusion_genes.py',kind='program')
        job.add('--input',outdir('candidate_fusion-genes_56.txt'),kind='input',temp_path=temp_flag)
        job.add('--label','focus',kind='parameter')
        job.add('--filter_gene_pairs',options.focus_fusions,kind='input')
        job.add('--output_fusion_genes',outdir('candidate_fusion-genes_57.txt'),kind='output')
        job.run()
    else:
        job.link(outdir('candidate_fusion-genes_56.txt'),outdir('candidate_fusion-genes_57.txt'),temp_path=temp_flag)
    # label fusion genes -- minimum distance between genes on the same strand
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_57.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','distance200kbp',kind='parameter')
    job.add('--min_dist_gene_gene','200000',kind='parameter')
    job.add('--min_dist_gene_gene_database',datadir('exons.txt'),kind='input')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_100.txt'),kind='output')
    job.run()
    #
    last_candidate_file = outdir('candidate_fusion-genes_100.txt')
    if options.label_file:
        title = options.label_title.strip().split(',')
        files = options.label_file.strip().split(',')
        thres = None
        if options.label_threshold:
            thres = options.label_threshold.strip().split(',')
        ainput = outdir('candidate_fusion-genes_custom___0.txt')
        job.link(last_candidate_file,
                 ainput,
                 temp_path = temp_flag)
        aout = ainput[:]
        for i in xrange(len(title)):
            ain = ainput.replace("___0.txt","___%s.txt" % (i,))
            aout = ainput.replace("___0.txt","___%s.txt" % (i+1,))
            job.add('label_fusion_genes.py',kind='program')
            job.add('--input',ain,kind='input',temp_path=temp_flag)
            job.add('--label',title[i],kind='parameter')
            if thres and thres != '0':
                job.add('--filter_gene_pairs_threshold',thres[i],kind='parameter')
            job.add('--filter_gene_pairs',files[i],kind='output')
            job.add('--output_fusion_genes',aout,kind='output')
            job.run()
        job.link(aout,
                 outdir('candidate_fusion-genes_custom___last.txt'),
                 temp_path = temp_flag)
    else:
        job.link(last_candidate_file,
                 outdir('candidate_fusion-genes_custom___last.txt'),
                 temp_path = temp_flag)
    # label fusion genes -- banned
    job.add('label_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_custom___last.txt'),kind='input',temp_path=temp_flag)
    job.add('--label','banned',kind='parameter')
    job.add('--filter_gene_pairs',datadir('banned.txt'),kind='output')
    job.add('--output_fusion_genes',outdir('candidate_fusion-genes_last.txt'),kind='output')
    job.run()

    ##############################################################################
    # FIND EXON-EXON JUNCTION IN FUSION GENES
    ##############################################################################
    # extract the relevant fusion genes from the found list for further analysis
    job.add('extract_fusion_genes.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_last.txt'),kind='input',temp_path=temp_flag)
    job.add('--input_fusion_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='output')
    job.add('--threshold_pairs',spanning_pairs_minimum,kind='parameter') # considers only the fusion genes candidates with more than 3 paired-end reads
    job.add('--threshold_count',options.spanning_pairs_count,kind='parameter')
    if options.biotypes_more:
        job.add('--skip_labels',options.biotypes+','+options.biotypes_more,kind='parameter') # skips the fusion genes candidates which are labeled
    else:
        job.add('--skip_labels',options.biotypes,kind='parameter') # skips the fusion genes candidates which are labeled
    if not options.skip_known_fusions:
        job.add('--allowed_labels','known_fusion,cosmic,ticdb,cgp',kind='parameter') # it allows the known fusions to be considered for further analysis
    if options.focus_fusions:
        job.add('--further_labels','focus',kind='parameter') # it allows the focus fusion genes to pass even when they are under the threshold
    job.add('--output',outdir('candidate_fusion-genes_exon-exon.txt'),kind='output')
    job.add('--output_fusion',outdir('candidate_fusion-genes_further.txt'),kind='output')
    job.add('--output_fusion_reads',outdir('candidate_fusion-genes_further_paired-reads.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('-F',kind='parameter')
    job.add('-f',datadir('custom_genes_mark.txt'),kind='input')
    job.add('',outdir('candidate_fusion-genes_further.txt'),kind='input')
    job.add('>',outdir('candidate_fusion-genes_further_mark.txt'),kind='output')
    job.run(successful_exit_status=(0,1))


    # save preliminary list of candidate fusion genes
    job.add('add_ambiguous_counts.py',kind='program')
    job.add('--input',outdir('candidate_fusion-genes_further.txt'),kind='input')
    job.add('--input_ambiguous',outdir('all_ambiguous_genes.txt'),kind='input')
    job.add('--output',outdir('preliminary-list_candidate-fusion-genes.txt'),kind='output')
    job.run()

    candidates = True
    if job.iff((empty(outdir('candidate_fusion-genes_exon-exon.txt')) or
                empty(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.fq'))),
                id = "#no-candidate-fusion-genes-found-1#"):
        candidates = False
        t = ["="*80,
             "WARNING: No candidate fusion genes have been found (due to no paired-reads",
             "         being found which support any possible fusion gene)!",
             "="*80
            ]
        job.write(t, stderr=True)
        if job.run():
            file(info_file,'a').writelines([el.rstrip('\r\n')+'\n' for el in [""]+t+[""]])


#        if (not candidates) and (not options.keep_viruses):
#            t = ["="*80,
#                 "WARNING: Viruses statistics and filtering has been skipped due to no fusions beeing found!",
#                 "         If one wants to have even in cases when no fusion genes are found then always ",
#                 "         (re)run FusionCatcher using '--keep-viruses-alignments' command line option!",
#                 "="*80
#                ]
#            job.write(t, stderr=True)
#            if job.run():
#                file(info_file,'a').writelines([el.rstrip('\r\n')+'\n' for el in [""]+t+[""]])


        # summary the exon-exon mappings
        # just get the header
        job.add('build_report_fusions_map.py',kind='program')
        job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BOWTIE.txt'), kind='output')
        job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BOWTIE.zip'), kind='output')
        job.run()


#        job.clean(outdir('reads_filtered_transcriptome_sorted-read.map'),temp_path=temp_flag)
        job.clean(outdir('original.fq.gz'),temp_path=temp_flag)
        job.clean(outdir('originala.fq.gz'),temp_path=temp_flag)
        job.clean(outdir('candidate_fusion-genes_missing_mates.txt'),temp_path=temp_flag)
        job.clean(outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),temp_path=temp_flag)
        #job.clean(outdir('candidate_fusion-genes_exon-exon.txt'))

    #if candidates or options.keep_viruses: # this has been removed in order to get always the viruses statistics
    if options.skip_viruses_filtering and (not options.keep_viruses) and (not candidates):
        t = ["="*80,
             "WARNING: Viruses statistics and filtering has been skipped due to no fusions!",
             "         being found and use of command line option '--skip-vir'!",
             "="*80
            ]
        job.write(t, stderr=True)
        if job.run():
            file(info_file,'a').writelines([el.rstrip('\r\n')+'\n' for el in [""]+t+[""]])

    if job.iff( not empty(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.fq')),
                id = "#no-candidate-fusion-genes-found-1A#"):
        #
        # FIND the FUSION POINT
        #
        if options.mismatches < 2: # stjude # before was 3
            ##############################################################################
            # FILTER the unmapped reads
            ##############################################################################

            # map on transcriptome again
            job.add('bowtie',kind='program')
            job.add('-t',kind='parameter')
            #job.add('-q',kind='parameter')
            job.add('-a',kind='parameter')
            job.add('-v','2',kind='parameter') #options.mismatches # stjude # before was 3
            job.add('-p',options.processes,kind='parameter',checksum='no')
            job.add('-m','1',kind='parameter')
            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
            job.add('--tryhard',kind='parameter')
            job.add('--best',kind='parameter')
            job.add('--strata',kind='parameter')
            job.add('--suppress','5,6,7',kind='parameter')
            job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_more.fq'),kind='output') # here is the result
            job.add('--max',outdir('reads_filtered_not-mapped_multiple.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
            if os.path.isfile(datadir('transcripts_index','.1.ebwtl')):
                job.add('--large-index',kind='parameter')
            job.add('',datadir('transcripts_index/'),kind='input')
            job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.fq'),kind='input',temp_path=temp_flag)
            #job.add('',outdir('reads-unmapped-filtered-trans.map'),kind='output')
            job.add('2>',outdir('log_bowtie_reads-unmapped-filtered-out-transcriptome.stdout.txt'),kind='output',checksum='no')
            #job.add('2>&1',kind='parameter',checksum='no') # XXX
            #job.run()
            job.add('|',kind='parameter')
            job.add('LC_ALL=C',kind='parameter')
            job.add('sort',kind='parameter')
            if sort_buffer:
                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
            if sort_parallel:
                job.add('--parallel',options.processes,kind='parameter',checksum='no')
            if sort_lzop_compress:
                job.add('--compress-program','lzop',kind='parameter',checksum='no')
            elif sort_gzip_compress:
                job.add('--compress-program','gzip',kind='parameter',checksum='no')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            #job.add('-s',kind='parameter') # stable sort
            job.add('-t',"'\t'",kind='parameter')
            job.add('-k','1,1',kind='parameter')
            #job.add('',outdir('reads-unmapped-filtered-trans.map'),kind='input',temp_path = temp_flag)
            job.add('>',outdir('reads-unmapped-filtered-trans.sorted.map'),kind='output')
            job.run()

            info(job,
                fromfile = outdir('log_bowtie_reads-unmapped-filtered-out-transcriptome.stdout.txt'),
                tofile = info_file,
                top = ["Mapping all reads (which do not map on genome and do not map on transcriptome) on transcriptome for filtering purposes:",
                       "-----------------------------------------------------------------------------------------------------------------------"],
                bottom = "\n\n\n",
                temp_path = temp_flag)

#            if pigz:
#                job.add('pigz',kind='program')
#                job.add('-p',options.processes,kind='parameter',checksum='no')
#            else:
#                job.add('gzip',kind='program')
#            job.add('-d',kind='parameter')
#            job.add('-c',outdir('reads_filtered_transcriptome_sorted-read.map.gz'),kind='input', temp_path = temp_flag)
#            job.add('|',kind='parameter')
            job.add('LC_ALL=C',kind='program')
            job.add('sort',kind='parameter')
            if sort_buffer:
                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
            if sort_parallel:
                job.add('--parallel',options.processes,kind='parameter',checksum='no')
            if sort_lzop_compress:
                job.add('--compress-program','lzop',kind='parameter',checksum='no')
            elif sort_gzip_compress:
                job.add('--compress-program','gzip',kind='parameter',checksum='no')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            #job.add('-s',kind='parameter') # stable sort
            job.add('-t',"'\t'",kind='parameter')
            job.add('-k','1,1',kind='parameter')
            job.add('-m',kind='parameter')
            job.add('',outdir('reads-unmapped-filtered-trans.sorted.map'),kind='input',temp_path = temp_flag)
            job.add('',outdir('reads_filtered_transcriptome_sorted-read.map'),kind='input',temp_path = temp_flag)
#            job.add('-',kind='parameter')
            job.add('>',outdir('reads_filtered_transcriptome_sorted-read_end.map'),kind='output')
            job.run()
#            job.add('|',kind='parameter')
#            if pigz:
#                job.add('pigz',kind='parameter')
#                job.add('-p',options.processes,kind='parameter',checksum='no')
#            else:
#                job.add('gzip',kind='parameter')
#            job.add('>',outdir('reads_filtered_transcriptome_sorted-read_end.map.gz'),kind='output')
#            job.run()

            # map on genome again
            job.add('bowtie',kind='program')
            job.add('-t',kind='parameter')
            #job.add('-q',kind='parameter')
            job.add('-a',kind='parameter')
            job.add('-v','2',kind='parameter') #options.mismatches # stjude # it was 3 before
            job.add('-p',options.processes,kind='parameter',checksum='no')
            job.add('-m','1',kind='parameter')
            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
            job.add('--tryhard',kind='parameter')
            job.add('--best',kind='parameter')
            job.add('--strata',kind='parameter')
#            job.add('--suppress','1,2,3,4,5,6,7,8',kind='parameter')
            job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end.fq'),kind='output') # here is the result
            job.add('--max',outdir('reads_filtered_not-mapped_multiple_end.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
            if os.path.isfile(datadir('genome_index','.1.ebwtl')):
                job.add('--large-index',kind='parameter')
            job.add('',datadir('genome_index/'),kind='input')
            job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_more.fq'),kind='input',temp_path=temp_flag)
            job.add('',outdir('reads-unmapped-filtered-geno.map'),kind='output',temp_path=temp_flag)
#            job.add('','/dev/null',kind='parameter')
            job.add('2>',outdir('log_bowtie_reads-unmapped-filtered-out-genome.stdout.txt'),kind='output',checksum='no')
            #job.add('2>&1',kind='parameter',checksum='no')
            job.run()

            job.clean(outdir('log_bowtie_reads-unmapped-filtered-out-genome.stdout.txt'),temp_path=temp_flag)

        else:
            job.link(outdir('reads_filtered_transcriptome_sorted-read.map'),
                     outdir('reads_filtered_transcriptome_sorted-read_end.map'),
                     temp_path = temp_flag)
            job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.fq'),
                     outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end.fq'),
                     temp_path = temp_flag,
                     kind = 'soft')


        # filter out the reads with poly tail
        # trim the poly tails
        job.add('trim_poly_tails.py',kind='program')
        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end.fq'),kind='input',temp_path=temp_flag)
        job.add('--repeats','12',kind='parameter') # 12
        #job.add('--skip_reads',kind='parameter')
        job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f2b.fq'),kind='output')
        job.run()

        # clip the low quality ends
        job.add('clip_quality.py',kind='program')
        job.add('--processes',options.processes,kind='parameter',checksum='no')
        job.add('-t',options.trim_quality,kind='parameter') # below Q5 trimming starts
        job.add('--score-type','sanger',kind='parameter')
        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f2b.fq'),kind='input',temp_path=temp_flag)
        job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f2.fq'),kind='output')
        job.run()

        # remove reads shorter than a given threshold
#        job.add('remove_shorter_reads.py',kind='program')
#        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f2.fq'),kind='input',temp_path=temp_flag)
#        job.add('--threshold',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
#        job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f3.fq'),kind='output')
#        job.run()
        job.add('seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-L',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f2.fq'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f3.fq'),kind='output')
        job.run()

        if not options.skip_str_filtering:
            # remove STR reads
            job.add('remove_str.py',kind='program')
            job.add('--processes',options.processes,kind='parameter',checksum='no')
            job.add('--threshold','2.1',kind='parameter',checksum='no')
            job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f3.fq'),kind='input',temp_path = temp_flag)
            job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4.fq'),kind='output')
    #        job.add('--str',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4-str.fq'),kind='output',temp_path = temp_flag)
            job.add('--log',outdir('log_reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4-str.txt'),kind='output')
            job.run()
            job.add('cat',kind='program')
            job.add('',outdir('log_reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4-str.txt'),kind='input',temp_path = temp_flag)
            job.add('>>',info_file,kind='output')
            job.run()
            job.add('printf',kind='program')
            job.add('"\n\n\n"',kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()
        else:
            job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f3.fq'),
                     outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4.fq'),
                     temp_path=temp_flag)


        # remove the reads which map on viruses genomes downloaded from NCBI
        if not options.skip_viruses_filtering:

            if options.keep_viruses:
                # convert the quality scores to Illumina Solexa version 1.5 format
                job.add('phred.py',kind='program')
                job.add('--link','soft',kind='parameter')
                job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4.fq'),kind='input')
                job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4-sam.fq'),kind='output')
                job.add('--input_type','auto-detect',kind='parameter')
                job.add('--output_type','sanger',kind='parameter')
                job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
                job.run()

                job.add('bowtie',kind='program')
                job.add('-t',kind='parameter')
                #job.add('-q',kind='parameter')
                job.add('-v',options.filter_mismatches,kind='parameter') #options.mismatches
                job.add('-p',options.processes,kind='parameter',checksum='no')
                #job.add('-m','1',kind='parameter')
                job.add('-a',kind='parameter')
                job.add('--best',kind='parameter')
                job.add('--strata',kind='parameter')
                job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                job.add('--sam',kind='parameter')
                #job.add('--tryhard',kind='parameter')
                #job.add('--best',kind='parameter')
                #job.add('--strata',kind='parameter')
                if len_reads > 40 and options.trim_wiggle:
                    job.add('--trim3',options.trim_wiggle,kind='parameter') # trim on the fly 5bp from 3' end
                    job.add('--trim5',options.trim_wiggle,kind='parameter') # trim the 5
                if os.path.isfile(datadir('viruses_index','.1.ebwtl')):
                    job.add('--large-index',kind='parameter')
                job.add('',datadir('viruses_index/'),kind='input')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4-sam.fq'),kind='input',temp_path=temp_flag)
                #job.add('',outdir('reads-mapped-on-viruses.sam'),kind='output')
                job.add('2>',outdir('log_sam-viruses.stdout.txt'),kind='output',checksum='no')
                job.add('|',kind='parameter')
                job.add('awk',"""'$3 == "*" { next } { print }'""",kind='parameter')
                job.add('|',kind='parameter')
                if pigz:
                    job.add('pigz',kind='parameter')
                    job.add('-p',options.processes,kind='parameter',checksum='no')
                else:
                    job.add('gzip',kind='parameter')
                job.add('--fast',kind='parameter')
                job.add('>',outdir('reads-mapped-on-viruses.sam.gz'),kind='output')
                job.run()

                info(job,
                    fromfile = outdir('log_sam-viruses.stdout.txt'),
                    tofile = info_file,
                    top = ["Mapping all reads on viruses/bacteria genomes in order to generate SAM file:",
                           "----------------------------------------------------------------------------"],
                    bottom = "\n\n\n",
                    temp_path = temp_flag)

                # remove unmapped reads
                #samtools view -hS -F 4 mapped_unmapped.sam > mapped_only.sam

            job.add('printf',kind='program')
            job.add(('"\n\nMapping all input reads on viruses genomes database for filtering purposes:\n'+
                     '--------------------------------------------------------------------------------\n"'),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()

            job.add('bowtie',kind='program')
            job.add('-t',kind='parameter')
            #job.add('-q',kind='parameter')
            #job.add('-a',kind='parameter')
            job.add('-v',options.filter_mismatches,kind='parameter') #options.mismatches
            job.add('-p',options.processes,kind='parameter',checksum='no')
            #job.add('-m','1',kind='parameter')
            job.add('-k','1',kind='parameter')
            #job.add('--solexa1.3-quals',kind='parameter')
            job.add('--phred33-quals',kind='parameter')
            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
#            job.add('--tryhard',kind='parameter')
            #job.add('--best',kind='parameter')
            #job.add('--strata',kind='parameter')
            if len_reads > 40 and options.trim_wiggle:
                job.add('--trim3',options.trim_wiggle,kind='parameter') # trim on the fly 5bp from 3' end
                job.add('--trim5',options.trim_wiggle,kind='parameter') # trim the 5
            job.add('--suppress','1,2,4,5,6,7,8',kind='parameter') # originally was: '2,3,4,5,6,7,8'
            job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f5.fq'),kind='output') # here is the result
            job.add('--max',outdir('reads-filtered_temp_multiple-viruses.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
            if os.path.isfile(datadir('viruses_index','.1.ebwtl')):
                job.add('--large-index',kind='parameter')
            job.add('',datadir('viruses_index/'),kind='input')
            job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4.fq'),kind='input',temp_path=temp_flag)
            #job.add('',outdir('reads-filtered-viruses.map'),kind='output') # XXX
            job.add('2>>',info_file,kind='parameter',checksum='no')
            #job.add('2>&1',kind='parameter',checksum='no')
            #job.run()
            job.add('|',kind='parameter')
            job.add('LC_ALL=C',kind='parameter')
            job.add('sort',kind='parameter')
            if sort_buffer:
                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
            if sort_parallel:
                job.add('--parallel',options.processes,kind='parameter',checksum='no')
            if sort_lzop_compress:
                job.add('--compress-program','lzop',kind='parameter',checksum='no')
            elif sort_gzip_compress:
                job.add('--compress-program','gzip',kind='parameter',checksum='no')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            #job.add('',outdir('reads-filtered-viruses.map'),kind='input',temp_path=temp_flag) # XXX
            job.add('|',kind='parameter')
            job.add('uniq',kind='parameter')
            job.add('-c',kind='parameter')
            job.add('|',kind='parameter')
            job.add('LC_ALL=C',kind='parameter')
            job.add('sort',kind='parameter')
            job.add('-rn',kind='parameter')
            job.add('-k','1,1',kind='parameter')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            job.add('>',outdir('viruses_bacteria_temp.txt'),kind='output')
            job.run()


            info(job,
                 fromfile = outdir('viruses_bacteria_temp.txt'),
                 tofile = info_file,
                 top = ["","","",
                        "Viruses found to have reads mapped on their genomes:",
                        "----------------------------------------------------",
                        "Reads_count\tOrganism"],
                 bottom = "\n\n\n")

            job.add('printf',kind='program')
            job.add('"Counts_of_mapping_reads\tVirus/Bacteria/Phage\n"', kind='parameter')
            job.add('>',outdir('viruses_bacteria_header.txt'),kind='output')
            job.run()

            job.add('cat',kind='program')
            job.add('',outdir('viruses_bacteria_header.txt'),kind='input',temp_path=temp_flag)
            job.add('',outdir('viruses_bacteria_temp.txt'),kind='input',temp_path=temp_flag)
            job.add('>',outdir('viruses_bacteria_phages.txt'),kind='output')
            job.run()

#            job.add('printf',kind='program')
#            job.add(('"\n\n\nViruses found to have reads mapped on their genomes:\n'+
#                     '---------------------------------------------------------\n'+
#                     'Reads_count\tOrganism\n"'),kind='parameter')
#            job.add('>>',info_file,kind='output')
#            job.run()
#            job.add('cat',kind='program')
#            job.add('',outdir('viruses_statistics.txt'),kind='input',temp_path=temp_flag)
#            job.add('>>',info_file,kind='output')
#            job.run()
#            job.add('printf',kind='program')
#            job.add('"\n\n\n"',kind='parameter')
#            job.add('>>',info_file,kind='output')
#            job.run()

        else:
            job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f4.fq'),
                     outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f5.fq'),
                     temp_path=temp_flag)


        # filter -- map on genome again the trimmed reads with 11 bp from 3' end
        job.add('bowtie',kind='program')
        job.add('-t',kind='parameter')
        #job.add('-q',kind='parameter')
#        job.add('-a',kind='parameter')
        job.add('-v','1',kind='parameter') #options.mismatches
        job.add('-p',options.processes,kind='parameter',checksum='no')
#        job.add('-m','1',kind='parameter')
        job.add('--tryhard',kind='parameter')
        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
#        job.add('--best',kind='parameter')
#        job.add('--strata',kind='parameter')
        if len_reads > 40:
            job.add('--trim5', length_anchor_minimum - 4, kind='parameter') # trim 11 bp on the fly
        else:
            job.add('--trim5', length_anchor_minimum - 7, kind='parameter') # trim 11 bp on the fly
#        job.add('--suppress','1,2,3,4,5,6,7,8',kind='parameter')
        job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end2.fq'),kind='output') # here is the result
        job.add('--max',outdir('reads_filtered_not-mapped_multiple_end2.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
        if os.path.isfile(datadir('genome_index','.1.ebwtl')):
            job.add('--large-index',kind='parameter')
        job.add('',datadir('genome_index/'),kind='input')
        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end-f5.fq'),kind='input',temp_path=temp_flag)
        job.add('',outdir('reads-unmapped-filtered-geno_last.map'),kind='output',temp_path=temp_flag)
#        job.add('','/dev/null',kind='parameter')
        job.add('2>',outdir('log_bowtie_reads-unmapped-filtered-out-genome_last.stdout.txt'),kind='output',checksum='no')
        #job.add('2>&1',kind='parameter',checksum='no')
        job.run()

        info(job,
            fromfile = outdir('log_bowtie_reads-unmapped-filtered-out-genome_last.stdout.txt'),
            tofile = info_file,
            top = ["Mapping all trimmed unmapped reads again on genome for filtering purposes:",
                   "-------------------------------------------------------------------------"],
            bottom = "\n\n\n",
            temp_path = temp_flag)


        job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end2.fq'),
                 outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),
                 temp_path=temp_flag)
#        # filter -- map on genome again the trimmed reads with 11 bp from 5' end
#        job.add('bowtie',kind='program')
#        job.add('-t',kind='parameter')
#        job.add('-q',kind='parameter')
##        job.add('-a',kind='parameter')
#        job.add('-v','1',kind='parameter') #options.mismatches
#        job.add('-p',options.processes,kind='parameter',checksum='no')
##        job.add('-m','1',kind='parameter')
#        job.add('--tryhard',kind='parameter')
#         job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
##        job.add('--best',kind='parameter')
##        job.add('--strata',kind='parameter')
#        job.add('--trim3','13',kind='parameter') # trim 11 bp on the fly
#        job.add('--suppress','2,3,4,5,6,7,8',kind='parameter')
#        job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='output') # here is the result
#        job.add('--max',outdir('reads_filtered_not-mapped_multiple_final.fq'),kind='output',temp_path=temp_flag) # if this is missing then these reads are going to '--un'
#        job.add('',datadir('genome_index/'),kind='input')
#        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_end2.fq'),kind='input',temp_path=temp_flag)
#        job.add('',outdir('reads-unmapped-filtered-geno_last.map'),kind='output',temp_path=temp_flag)
#        job.add('>',outdir('log_bowtie_reads-unmapped-filtered-out-genome_last.stdout.txt'),kind='parameter',checksum='no')
#        job.add('2>&1',kind='parameter',checksum='no')
#        job.run()

    if job.iff(empty(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq')),
               id="#no-candidate-fusion-genes-found-2#"):

        if candidates:
            t = ["="*80,
                 "WARNING: No candidate fusion genes have been found (no unmapped reads found)!",
                 "="*80
                ]
            job.write(t, stderr=True)
            if job.run():
                file(info_file,'a').writelines([el.rstrip('\r\n')+'\n' for el in [""]+t+[""]])

            # summary the exon-exon mappings
            # just get the header
            job.add('build_report_fusions_map.py',kind='program')
            job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BOWTIE.txt'), kind='output')
            job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BOWTIE.zip'), kind='output')
            job.run()

            if options.keep_unmapped_reads:
                job.add('echo',kind='program')
                job.add('-n','""',kind='parameter')
                job.add('',outdir('unmapped-reads.fq.gz'), kind='output')
                job.run()
        candidates = False

        job.clean(outdir('reads_filtered_transcriptome_sorted-read.map'),temp_path=temp_flag)
        #
        #job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),
        #         outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq')
        #         )

    elif (not candidates) and options.keep_unmapped_reads:


        # convert FASTQ illumina to sanger
#        job.add('seqtk',kind='program')
#        job.add('seq',kind='parameter')
#        job.add('-Q64',kind='parameter')
#        job.add('-V',kind='parameter')
#        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='input',temp_path=temp_flag)
        #job.add('>',outdir('unmapped-reads.fq'),kind='output')
        #job.run()
        #job.add('|',kind='parameter')
        if pigz:
            job.add('pigz',kind='program')
            job.add('-p',options.processes,kind='parameter',checksum='no')
        else:
            job.add('gzip',kind='parameter')
        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='input',temp_path=temp_flag)
        job.add('--fast',kind='program')
        #job.add('-c',outdir('unmapped-reads.fq'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('unmapped-reads.fq.gz'),kind='output')
        job.run()


#        job.add('phred.py',kind='program')
#        job.add('--link','hard',kind='parameter')
#        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='input',temp_path=temp_flag)
#        job.add('--output',outdir('unmapped-reads.fq'),kind='output')
#        job.add('--input_type','illumina',kind='parameter')
#        job.add('--output_type','sanger',kind='parameter')
#        job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
#        job.run()
#
#        if pigz:
#            job.add('pigz',kind='program')
#            job.add('-p',options.processes,kind='parameter',checksum='no')
#        else:
#            job.add('gzip',kind='program')
#        job.add('-c',outdir('unmapped-reads.fq'),kind='input',temp_path=temp_flag)
#        job.add('>',outdir('unmapped-reads.fq.gz'),kind='output')
#        job.run()




    if candidates:

        # extract reads ids
        job.add('awk',kind='program')
        job.add("'NR%4==1 {print substr($0,2)}'",outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='input')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),kind='output')
        job.run()

        # add also the mates to the unmapped reads (the mates may be mapping just fine)
        job.add('awk',kind='program')
        job.add("""'{n=length($0); r=substr($0,1,n-1); print r"1"; print r"2"}'""",kind='parameter')
        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),kind='input',temp_path=temp_flag if not options.keep_unmapped_reads else 'no')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final2.txt'),kind='input')
        job.run()

        job.add('cat',kind='program')
        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final2.txt'),kind='input')
        job.add('',outdir('candidate_fusion-genes_further_paired-reads.txt'),kind='input')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        if sort_buffer:
            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
        if sort_parallel:
            job.add('--parallel',options.processes,kind='parameter',checksum='no')
        if sort_lzop_compress:
            job.add('--compress-program','lzop',kind='parameter',checksum='no')
        elif sort_gzip_compress:
            job.add('--compress-program','gzip',kind='parameter',checksum='no')
        job.add('-T',tmp_dir,kind='parameter',checksum='no')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('>',outdir('original_important.txt'),kind='output')
        job.run()

        if not options.split_seqtk_subseq:
            job.add('extract_short_reads.py',kind='program')
            job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
            job.add('--input',outdir('originala.fq.gz'),kind='input')
            job.add('--list',outdir('original_important.txt'),kind='input')
            job.add('--output',outdir('original_important.fq.gz'),kind='output')
            job.run(error_message = ("If this fails (again?) due to a memory error (e.g. not enough free memory) then lowering the "+
                                     "buffer size for specifically this script might help. This can be done by using the FusionCatcher's "+
                                     "command line option '--extra-buffer-size "+str(int(options.extract_buffer_size)/2)+"' ."))
        else:
            job.add('seqtk',kind='program')
            job.add('subseq',kind='parameter')
            job.add('',outdir('originala.fq.gz'),kind='input')
            job.add('',outdir('original_important.txt'),kind='input')
            job.add('|',kind='parameter')
            if pigz:
                job.add('pigz',kind='parameter')
                job.add('-p',options.processes,kind='parameter',checksum='no')
            else:
                job.add('gzip',kind='parameter')
            job.add('--fast',kind='parameter')
            #job.add('-c',outdir('unmapped-reads.fq'),kind='input',temp_path=temp_flag)
            job.add('>',outdir('original_important.fq.gz'),kind='output')
            job.run()

        job.clean(outdir('originala.fq.gz'),temp_path=temp_flag)

        # extract the line with reads which are important further from the transcriptome mapping
#        job.add('reads_from_map.py',kind='program')
#        job.add('--input_reads',outdir('original_important.txt'),kind='input',temp_path=temp_flag)
#        job.add('--input_map',outdir('reads_filtered_transcriptome_sorted-read_end.map'),kind='input',temp_path=temp_flag)
#        job.add('--output_map',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind='output')
#        job.add('--operation','extract',kind='parameter')
#        job.run()
        job.add('LC_ALL=C',kind='program')
        job.add('join',kind='parameter')
        job.add('-1','1',kind='parameter')
        job.add('-2','1',kind='parameter')
        job.add('-t',"'\t'",kind='parameter')
        job.add('',outdir('original_important.txt'),kind='input',temp_path=temp_flag)
        job.add('',outdir('reads_filtered_transcriptome_sorted-read_end.map'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind='output')
        job.run()



        if options.keep_unmapped_reads:

            #
#            job.add('extract_reads_ids.py',kind='program')
#            job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='input',temp_path=temp_flag)
#            job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),kind='output')
#            job.run()

#            extract the short reads which mapped on the transcriptome and do not map on genome
#            job.add('extract_short_reads.py',kind='program')
#            job.add('--input',outdir('original.fq.gz'),kind='input')
#            job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),kind='input',temp_path=temp_flag)
#            job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_long.fq'),kind='output')
#            job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
#            job.run(error_message = ("If this fails due to a memory error then lowering the "+
#                         "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
#                         "of FusionCatcher and running it again might help!"))

#            job.add('phred.py',kind='program')
#            job.add('--link','hard',kind='parameter')
#            job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_long.fq'),kind='input',temp_path=temp_flag)
#            job.add('--output',outdir('unmapped-reads.fq'),kind='output')
#            job.add('--input_type','illumina',kind='parameter')
#            job.add('--output_type','sanger',kind='parameter')
#            job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
#            job.run()

            if not options.split_seqtk_subseq:
                job.add('extract_short_reads.py',kind='program')
                job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                job.add('--input',outdir('original_important.fq.gz'),kind='input')
                job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),kind='input')
                job.add('--output','-',kind='parameter',checksum='no')
            else:
                job.add('seqtk',kind='program')
                job.add('subseq',kind='parameter')
                job.add('',outdir('original_important.fq.gz'),kind='input')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),kind='input')
            #job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_long.fq'),kind='output')
            #job.run()
            # convert FASTQ illumina to sanger
            #job.add('seqtk',kind='program')
#            job.add('|',kind='parameter')
#            job.add('seqtk',kind='parameter')
#            job.add('seq',kind='parameter')
#            job.add('-Q64',kind='parameter')
#            job.add('-V',kind='parameter')
#            job.add('-',kind='parameter')
            #job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_long.fq'),kind='input',temp_path=temp_flag)
            #job.add('>',outdir('unmapped-reads.fq'),kind='output')
            #job.run()
            job.add('|',kind='parameter')
            if pigz:
                job.add('pigz',kind='parameter')
                job.add('-p',options.processes,kind='parameter',checksum='no')
            else:
                job.add('gzip',kind='parameter')
            job.add('--fast',kind='parameter')
            #job.add('-c',outdir('unmapped-reads.fq'),kind='input',temp_path=temp_flag)
            job.add('>',outdir('unmapped-reads.fq.gz'),kind='output')
            job.run()

            job.clean(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.txt'),temp_path=temp_flag)

        if not options.all_reads_junction:
            job.add('remove_reads_exon_exon_fastq.py',kind='program')
            job.add('--input_fastq',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),kind='input',temp_path=temp_flag)
            job.add('--input_fusions',outdir('candidate_fusion-genes_exon-exon.txt'),kind='input')
            job.add('--input_transcriptome',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind = 'input')
            job.add('--output_fastq',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_plus.fq'),kind='output')
            job.add('--log',info_file,kind='output',checksum='no')
            job.run()
        else:
            job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),
                     outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_plus.fq'),
                     temp_path=temp_flag)

        # generate the exon-exon junctions
        job.add('generate_exon-exon_junctions.py',kind='program')
        job.add('--input_fusion_genes',outdir('candidate_fusion-genes_exon-exon.txt'),kind='input')
        job.add('--input_fasta_transcripts',datadir('transcripts.fa'),kind='input')
        job.add('--input_database_transcripts',datadir('transcripts.txt'),kind='input')
        job.add('--overlap_read',length_anchor_bowtie,kind='parameter') # :-)
        job.add('--unique_cut_sequences_same_pair',kind='parameter') #added
        job.add('--length_reads_filename',outdir('log_lengths_reads.txt'),kind='input')
        job.add('--output_cut_junction',outdir('exon-exon_junction_cut.fa'),kind='output')
        job.add('--output_count_seq',outdir('exon-exon_junction_cut__seq.txt'),kind='output')
        job.add('--output_count_nuc',outdir('exon-exon_junction_cut__nuc.txt'),kind='output')
        job.run()

        nucleotides_ee = int(file(outdir('exon-exon_junction_cut__nuc.txt'),"r").readline().strip())

        if nucleotides_ee > options.limit_bowtie:

            job.add('split-fasta.py',kind='program')
            job.add('--size',outdir('exon-exon_junction_cut__nuc.txt'),kind='input')
            job.add('--seqs',outdir('exon-exon_junction_cut__seq.txt'),kind='input')
            job.add('--threshold',options.limit_bowtie,kind='parameter')
            job.add('-i',outdir('exon-exon_junction_cut.fa'),kind='input')
            job.add('-o',outdir('exon-exon_junction_cut_split.fa'),kind='output')
            job.run()

            parts = [el.strip() for el in file(outdir('exon-exon_junction_cut_split.fa'),'r').readlines()]
            for i,part in enumerate(parts):
                # map the reads which do not align anywhere on the exon-exon junctions from fusion-genes
                # build index
                job.add('bowtie-build',kind='program')
                job.add('-f',kind='parameter')
                job.add('--quiet',kind='parameter')
#                job.add('--ntoa',kind='parameter')
                job.add('--offrate','1',kind='parameter')
                job.add('--ftabchars','5',kind='parameter')
                job.add('',part,kind='input')
                job.add('',part+'_dir/',kind='output')
                job.run()
                # map using the exon-exon fusion genes index (all possible mappings)
                job.add('bowtie',kind='program')
                job.add('-t',kind='parameter')
                #job.add('-q',kind='parameter')
                job.add('-a',kind='parameter')
                job.add('-v',options.mismatches,kind='parameter')
                job.add('-p',options.processes,kind='parameter',checksum='no')
                if os.path.isfile(os.path.join(part+'_dir','.1.ebwtl')):
                    job.add('--large-index',kind='parameter')
                job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                #job.add('--solexa1.3-quals',kind='parameter')
                job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq'),kind='output') # here is the result
                job.add('--tryhard',kind='parameter')
                #job.add('--best',kind='parameter')
                job.add('',part+'_dir/',kind='input',temp_path=temp_flag)
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_plus.fq'),kind='input',temp_path=temp_flag)
                #job.add('',outdir('reads_mapped-exon-exon-fusion-genes.map'),kind='output') # <== mappings on exon-exon junctions for fusion genes ####### # XXX
                job.add('2>',outdir('log_bowtie_reads_mapped-exon-exon-fusion-genes_map.stdout.txt.'+str(i)),kind='output',checksum='no')
                #job.add('2>&1',kind='parameter',checksum='no')
                #job.run()
                job.add('|',kind='parameter')
                # sort the reads' mappings on exon-exon by reference sequence, i.e.
                # gene-gene,transcript-transcript,exon-exon
                job.add('LC_ALL=C',kind='parameter')
                job.add('sort',kind='parameter')
                if sort_buffer:
                    job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                if sort_parallel:
                    job.add('--parallel',options.processes,kind='parameter',checksum='no')
                if sort_lzop_compress:
                    job.add('--compress-program','lzop',kind='parameter',checksum='no')
                elif sort_gzip_compress:
                    job.add('--compress-program','gzip',kind='parameter',checksum='no')
                job.add('-T',tmp_dir,kind='parameter',checksum='no')
                #job.add('-s',kind='parameter') # stable sort
                job.add('-t',"'\t'",kind='parameter')
                job.add('-k','3,3',kind='parameter')
                #job.add('',outdir('reads_mapped-exon-exon-fusion-genes.map'),kind='input',temp_path = temp_flag) # XXX
                job.add('>',outdir('reads_mapped-exon-exon-fusion-genes_sorted-ref.map.'+str(i)),kind='output',dest_list='exonexon')
                job.run()

                job.clean(outdir('log_bowtie_reads_mapped-exon-exon-fusion-genes_map.stdout.txt.'+str(i)),temp_path=temp_flag)

            job.sink(job.exonexon, outdir('exonexon.txt'))

            job.add('concatenate.py',kind='program')
            job.add('-f',outdir('exonexon.txt'),kind='input',temp_path=temp_flag)
            job.add('',outdir('reads_mapped-exon-exon-fusion-genes_sorted-ref.map'),kind='output')
            job.run()

#            for ft in job.exonexon:
#                job.clean(ft,temp_path=temp_flag)
            job.clean(job.exonexon,temp_path=temp_flag)
            
        else:
            # map the reads which do not align anywhere on the exon-exon junctions from fusion-genes
            # build index
            job.add('bowtie-build',kind='program')
            job.add('-f',kind='parameter')
            job.add('--quiet',kind='parameter')
#            job.add('--ntoa',kind='parameter')
            job.add('--offrate','1',kind='parameter')
            job.add('--ftabchars','5',kind='parameter')
            job.add('',outdir('exon-exon_junction_cut.fa'),kind='input')
            job.add('',outdir('exon-exon_fusion-genes/'),kind='output')
            job.run()
            # map using the exon-exon fusion genes index (all possible mappings)
            job.add('bowtie',kind='program')
            job.add('-t',kind='parameter')
            #job.add('-q',kind='parameter')
            job.add('-a',kind='parameter')
            job.add('-v',options.mismatches,kind='parameter')
            job.add('-p',options.processes,kind='parameter',checksum='no')
            if os.path.isfile(outdir('exon-exon_fusion-genes','.1.ebwtl')):
                job.add('--large-index',kind='parameter')
            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
            #job.add('--solexa1.3-quals',kind='parameter')
            job.add('--un',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq'),kind='output') # here is the result
            job.add('--tryhard',kind='parameter')
            #job.add('--best',kind='parameter')
            job.add('',outdir('exon-exon_fusion-genes/'),kind='input',temp_path=temp_flag)
            job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final_plus.fq'),kind='input',temp_path=temp_flag)
            #job.add('',outdir('reads_mapped-exon-exon-fusion-genes.map'),kind='output') # <== mappings on exon-exon junctions for fusion genes ####### # XXX
            job.add('2>',outdir('log_bowtie_reads_mapped-exon-exon-fusion-genes_map.stdout.txt'),kind='output',checksum='no')
            #job.add('2>&1',kind='parameter',checksum='no')
            #job.run()
            job.add('|',kind='parameter')
            # sort the reads' mappings on exon-exon by reference sequence, i.e.
            # gene-gene,transcript-transcript,exon-exon
            job.add('LC_ALL=C',kind='parameter')
            job.add('sort',kind='parameter')
            if sort_buffer:
                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
            if sort_parallel:
                job.add('--parallel',options.processes,kind='parameter',checksum='no')
            if sort_lzop_compress:
                job.add('--compress-program','lzop',kind='parameter',checksum='no')
            elif sort_gzip_compress:
                job.add('--compress-program','gzip',kind='parameter',checksum='no')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            #job.add('-s',kind='parameter') # stable sort
            job.add('-t',"'\t'",kind='parameter')
            job.add('-k','3,3',kind='parameter')
            #job.add('',outdir('reads_mapped-exon-exon-fusion-genes.map'),kind='input',temp_path = temp_flag) # XXX
            job.add('>',outdir('reads_mapped-exon-exon-fusion-genes_sorted-ref.map'),kind='output')
            job.run()

            job.clean(outdir('log_bowtie_reads_mapped-exon-exon-fusion-genes_map.stdout.txt'),temp_path=temp_flag)

        # more filtering -- remove the reads from the exon-exon junctions which
        # have the pair read mapping on a totally different gene than those
        # involved in the exon-exon junction
        job.add('remove_reads_exon_exon_map.py',kind='program')
        if not options.all_reads_junction:
            job.add('--only_pairs',kind='parameter')
        job.add('--input_exon_exon',outdir('reads_mapped-exon-exon-fusion-genes_sorted-ref.map'),kind='input',temp_path=temp_flag)
        job.add('--input_transcriptome',
                outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),
                kind = 'input',
                temp_path = 'no' if (((not options.skip_blat) or (not options.skip_star) or (not options.skip_bowtie2) or (not options.skip_bwa)) and (not options.all_reads_junction)) else temp_flag)
        job.add('--output',outdir('reads_mapped-exon-exon-fusion-genes_sorted-ref_filtered.map'),kind='output')
        job.run()

        # analyze the exon-exon mappings
        job.add('analyze_exon-exon_mappings.py',kind='program')
        job.add('--input',outdir('reads_mapped-exon-exon-fusion-genes_sorted-ref_filtered.map'),kind='input',temp_path=temp_flag)
        job.add('--input_hugo',datadir('genes_symbols.txt'),kind='input')
        job.add('--output',outdir('candidate_fusion-genes_exon-exon-junctions_summary.txt'),kind='output')
        job.add('--output_henrik',outdir('candidate_fusion-genes_exon-exon-junctions_reads-positions.txt'),kind='output')
        job.run()

        # summary the exon-exon mappings
        job.add('build_report_fusions_map.py',kind='program')
        job.add('--suporting_unique_reads',spanning_reads_bowtie,kind='parameter')
        job.add('--anchor2',length_anchor2,kind='parameter')
        job.add('--input_exons',datadir('exons.txt'),kind='input')
        job.add('--input_candidate_fusion_genes',outdir('candidate_fusion-genes_further.txt'),kind='input')
        job.add('--input_fusion_summary',outdir('candidate_fusion-genes_exon-exon-junctions_summary.txt'),kind='input',temp_path=temp_flag)
        job.add('--input_fusion_summary_more',outdir('candidate_fusion-genes_exon-exon-junctions_reads-positions.txt'),kind='input',temp_path=temp_flag)
        job.add('--input_candidate_fusions_missing_mates',outdir('candidate_fusion-genes_missing_mates.txt'),kind='input',temp_path=temp_flag)
        job.add('--input_fasta_juncs',outdir('exon-exon_junction_cut.fa'),kind='input',temp_path=temp_flag)
#        if options.reads_preliminary_fusions:
#            job.add('--output_all_candidate_fusion_genes_reads',outdir('pre-fusion'),kind='output')
        job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
        if (not options.skip_blat) or (not options.skip_star) or (not options.skip_bowtie2) or (not options.skip_bwa):
            job.add('--input_fastq',outdir('original_important.fq.gz'),kind='input') # needed for BLAT later
            job.add('--input_candidate_fusion_genes_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='input')
        else:
            job.add('--input_fastq',outdir('original_important.fq.gz'),kind='input',temp_path=temp_flag)
            job.add('--input_candidate_fusion_genes_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='input',temp_path=temp_flag)
        if options.psl_visualization and not empty(datadir('genome.2bit')):
            job.add('--input_genome_2bit',datadir('genome.2bit'),kind='input')
            job.add('--psl_alignment_type','web',kind='parameter')
        if options.sam_visualization:
            job.add('--input_genome_bowtie2',datadir('genome_index2/index'),kind='input')
            job.add('--sam_alignment','20',kind='parameter')
            job.add('--threads',options.processes,kind='parameter')
        if options.assembly:
            job.add('--velvet',kind='parameter')
        job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BOWTIE.txt'),kind='output')
        job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BOWTIE.zip'),kind='output')
        job.run()


##################################################################################
##################################################################################
##################################################################################
    # Find fusion genes using BLAT/STAR/BOWTIE2
##################################################################################
##################################################################################
##################################################################################
##################################################################################

    if ((not options.skip_blat) or (not options.skip_star) or (not options.skip_bowtie2) or (not options.skip_bwa)) and candidates:
        # generate the gene-gene junctions in FASTA format
        job.add('generate_gene-gene_junctions.py',kind='program')
        job.add('--input',outdir('candidate_fusion-genes_exon-exon.txt'),kind='input',temp_path=temp_flag)
        job.add('--input_database',datadir('genes.fa'),kind='input')
        job.add('--input_exons',datadir('exons.txt'),kind='input')
        job.add('--reverse',kind='parameter')
        job.add('--longest',outdir('gene-gene_longest.txt'),kind='output')
        job.add('--output',outdir('gene-gene.fa'),kind='output')
        job.add('--output_genes',outdir('gene-gene_unique.fa'),kind='output')
        job.add('--output_genes_count_seq',outdir('gene-gene_unique__seq.txt'),kind='output')
        job.add('--output_genes_count_nuc',outdir('gene-gene_unique__nuc.txt'),kind='output')
        job.run()

        nucleotides_ggu = int(file(outdir('gene-gene_unique__nuc.txt'),'r').readline().strip())

        job.add('du',kind='program')
        job.add('-b',outdir('gene-gene.fa'),kind='input')
        job.add('|',kind='parameter')
        job.add('cut',kind='parameter')
        job.add('-f','1',kind='parameter')
        job.add('>',outdir('gene-gene__nuc.txt'),kind='output')
        job.run()

        nucleotides_gg = int(file(outdir('gene-gene__nuc.txt'),'r').readline().strip())


        job.add('grep',kind='program')
        job.add('-c',kind='parameter')
        job.add("'^>'",outdir('gene-gene.fa'),kind='input')
        job.add('>',outdir('gene-gene__seq.txt'),kind='output')
        job.run(successful_exit_status=(0,1))

        sequences_gg = int(file(outdir('gene-gene__seq.txt'),'r').readline().strip())

        if options.trim_psl:
            # link
            job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq'),
                     outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl.fq'),
                     temp_path = temp_flag)
        else:
#            job.add('extract_reads_ids.py',kind='program')
#            job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq'),kind='input',temp_path=temp_flag)
#            job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.txt'),kind='output')
#            job.run()
            # extract reads ids
            job.add('awk',kind='program')
            job.add("'NR%4==1 {print substr($0,2)}'",outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq'),kind='input',temp_path=temp_flag)
            job.add('|',kind='parameter')
            job.add('uniq',kind='parameter')
            job.add('|',kind='parameter')
            job.add('LC_ALL=C',kind='parameter')
            job.add('sort',kind='parameter')
            if sort_buffer:
                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
            if sort_parallel:
                job.add('--parallel',options.processes,kind='parameter',checksum='no')
            if sort_lzop_compress:
                job.add('--compress-program','lzop',kind='parameter',checksum='no')
            elif sort_gzip_compress:
                job.add('--compress-program','gzip',kind='parameter',checksum='no')
            job.add('-T',tmp_dir,kind='parameter',checksum='no')
            job.add('|',kind='parameter')
            job.add('uniq',kind='parameter')
            job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.txt'),kind='output')
            job.run()

#            extract the short reads which mapped on the transcriptome and do not map on genome
#            job.add('extract_short_reads.py',kind='program')
#            job.add('--input',outdir('original.fq.gz'),kind='input')
#            job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.txt'),kind='input',temp_path=temp_flag)
#            job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_original.fq'),kind='output')
#            job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
#            job.run(error_message = ("If this fails due to a memory error then lowering the "+
#                         "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
#                         "of FusionCatcher and running it again might help!"))

            if min(max_len_reads,options.trim_psl_3end_keep) > min_len_reads + 9:
                # add also the reads of the paired reads which support the candidate fusion genes (some of them might overlap the fusion junction)
                job.add('cat',kind='program')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.txt'),kind='input',temp_path=temp_flag)
                job.add('',outdir('candidate_fusion-genes_further_paired-reads.txt'),kind='input')
                job.add('|',kind='parameter')
                job.add('uniq',kind='parameter')
                job.add('|',kind='parameter')
                job.add('LC_ALL=C',kind='parameter')
                job.add('sort',kind='parameter')
                if sort_buffer:
                    job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                if sort_parallel:
                    job.add('--parallel',options.processes,kind='parameter',checksum='no')
                if sort_lzop_compress:
                    job.add('--compress-program','lzop',kind='parameter',checksum='no')
                elif sort_gzip_compress:
                    job.add('--compress-program','gzip',kind='parameter',checksum='no')
                job.add('-T',tmp_dir,kind='parameter',checksum='no')
                job.add('|',kind='parameter')
                job.add('uniq',kind='parameter')
                job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome2.txt'),kind='output')
                job.run()
            else:
                job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.txt'),
                         outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome2.txt'),
                         temp_path = temp_flag)

            if not options.split_seqtk_subseq:
                job.add('extract_short_reads.py',kind='program')
                job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                job.add('--input',outdir('original_important.fq.gz'),kind='input')
                job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome2.txt'),kind='input')
                job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_temp.fq'),kind='output')
                job.run(error_message = ("If this fails (again?) due to a memory error (e.g. not enough free memory) then lowering the "+
                                         "buffer size for specifically this script might help. This can be done by using the FusionCatcher's "+
                                         "command line option '--extra-buffer-size "+str(int(options.extract_buffer_size)/2)+"' ."))
            elif options.split_seqtk_subseq == 1:
                job.add('seqtk',kind='program')
                job.add('subseq',kind='parameter')
                job.add('',outdir('original_important.fq.gz'),kind='input')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome2.txt'),kind='input')
                job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_temp.fq'),kind='output')
                job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
                           "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))
            elif options.split_seqtk_subseq > 1:
                job.add('seqtk-subseq.sh',kind='program')
                job.add('',options.split_seqtk_subseq,kind='parameter')
                job.add('',outdir('original_important.fq.gz'),kind='input')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome2.txt'),kind='input')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_temp.fq'),kind='output')
                job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
                           "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))


            job.clean(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome2.txt'),temp_path=temp_flag)

            # link
            #job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_original.fq'),
            #         outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl.fq'),
            #         temp_path=temp_flag)

            job.add('fastq_b2n.py',kind='program')
            job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_temp.fq'),kind='input',temp_path=temp_flag)
            job.add('--replacement','A',kind='parameter')
            job.add('--ambiguous',kind='parameter')
            job.add('--sanger',kind='parameter')
            job.add('--threshold',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file = 'yes')
            job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl.fq'),kind='output')
            job.run()


        job.add('trim_poly_tails.py',kind='program')
        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl.fq'),kind='input',temp_path=temp_flag)
        job.add('--repeats','12',kind='parameter')
        #job.add('--skip_reads',kind='parameter')
        #job.add('--replace',kind='parameter') # test eml4
        job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus_b.fq'),kind='output')
        job.run()

        # clip the low quality ends
        job.add('clip_quality.py',kind='program')
        job.add('--processes',options.processes,kind='parameter',checksum='no')
        job.add('-t',options.trim_quality,kind='parameter') # below Q5 trimming starts
        job.add('--score-type','sanger',kind='parameter')
        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus_b.fq'),kind='input',temp_path=temp_flag)
        job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus_temp.fq'),kind='output')
        job.run()

        if options.trim_psl_5end and options.trim_5end > 0:
            # trim 5
            job.add('seqtk',kind='program')
            job.add('trimfq',kind='parameter')
            job.add('-l','1',kind='parameter')
            job.add('-b',options.trim_5end,kind='parameter')
            job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus_temp.fq'),kind='input', temp_path='no')
            job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus.fq'),kind='output')
            job.run()
        else:
            job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus_temp.fq'),
                     outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus.fq'),
                     temp_path = temp_flag)

        # remove reads shorter than a given threshold
#        job.add('remove_shorter_reads.py',kind='program')
#        job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus.fq'),kind='input',temp_path=temp_flag)
#        job.add('--threshold',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
#        job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq'),kind='output')
#        job.run()
        job.add('seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-L',outdir('log_minimum_length_short_read.txt'),kind='parameter',from_file='yes')
        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-plus.fq'),kind='input',temp_path=temp_flag)
        job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq'),kind='output')
        job.run()


        if job.iff(empty(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq')),
                   id = "#reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq#"):
            candidates = False

            t = ["="*80,
                 "WARNING: No candidate fusion genes have been found using BLAT/STAR/BOWTIE2/BWA!",
                 "="*80
                ]
            job.write(t, stderr=True)
            if job.run():
                file(info_file,'a').writelines([el.rstrip('\r\n')+'\n' for el in [""]+t+[""]])

            if options.keep_unmapped_reads:
                job.add('echo',kind='program')
                job.add('-n','""',kind='parameter')
                job.add('>',outdir('unmapped-reads.fq.gz'), kind='output')
                job.run()

            # summary the exon-exon mappings
            # just get the header
            if not options.skip_blat:
                job.add('build_report_fusions_map.py',kind='program')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BLAT.txt'), kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BLAT.zip'), kind='output')
                job.run()

            if not options.skip_star:
                job.add('build_report_fusions_map.py',kind='program')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_STAR.txt'), kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_STAR.zip'), kind='output')
                job.run()

            if not options.skip_bowtie2:
                job.add('build_report_fusions_map.py',kind='program')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BOWTIE2.txt'), kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BOWTIE2.zip'), kind='output')
                job.run()

            if not options.skip_bwa:
                job.add('build_report_fusions_map.py',kind='program')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BWA.txt'), kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BWA.zip'), kind='output')
                job.run()

        else:

            job.add('printf',kind='program')
            job.add(('"\nCounts of reads before BLAT/STAR/BOWTIE2/BWA alignment:\n'+
                     '----------------------------------------------------------\n"'),kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()
            job.add('cat',kind='program')
            job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq'),kind='input')
            job.add('|',kind='parameter')
            job.add("echo $((`wc -l`/4))",kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()
            job.add('printf',kind='program')
            job.add('"\n\n\n"',kind='parameter')
            job.add('>>',info_file,kind='output')
            job.run()

            # trim the reads given as input to BLAT
            input_file = outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq')
            output_file = outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq')
            if options.trim_psl_3end_keep > 25 and ((options.trim_psl_3end_keep < max_len_reads and (not options.trim_psl)) or (options.trim_psl and min_len_reads > options.trim_psl_3end_keep)):
                # trim from 3-end to have the reads all the same length
#                job.add('trim_reads.py',kind='program')
#                job.add('--input',input_file,kind='input',temp_path=temp_flag)
#                job.add('--output',output_file,kind='output')
#                job.add('--trim_end','3',kind='parameter')
#                job.add('--final_size',options.trim_psl_3end_keep,kind='parameter')
#                job.run()
#
                # original
#                job.add('seqtk',kind='program')
#                job.add('trimfq',kind='parameter')
#                job.add('-l','1',kind='parameter')
#                #job.add('-q','0',kind='parameter')
#                job.add('-B',options.trim_psl_3end_keep,kind='parameter')
#                job.add('',input_file,kind='input', temp_path=temp_flag)
#                job.add('>',output_file,kind='output')
#                job.run()

                # try to trim only the unmapped read (do not trim the paired reads supporting the fusions)
                if not options.split_seqtk_subseq:
                    job.add('extract_short_reads.py',kind='program')
                    job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                    job.add('--input',input_file,kind='input')
                    job.add('--list',outdir('candidate_fusion-genes_further_paired-reads.txt'),kind='input')
                    job.add('--output',outdir('reads_not-for-trimming.fq'),kind='output')
                    job.run(error_message = ("If this fails (again?) due to a memory error (e.g. not enough free memory) then lowering the "+
                                             "buffer size for specifically this script might help. This can be done by using the FusionCatcher's "+
                                             "command line option '--extra-buffer-size "+str(int(options.extract_buffer_size)/2)+"' ."))
                else:
                    job.add('seqtk',kind='program')
                    job.add('subseq',kind='parameter')
                    job.add('',input_file,kind='input')
                    job.add('',outdir('candidate_fusion-genes_further_paired-reads.txt'),kind='input')
                    job.add('>',outdir('reads_not-for-trimming.fq'),kind='output')
                    job.run(error_message=("ERROR: Most likely this fails because there is not enough free RAM memory for running SEQTK SUBSEQ tool <https://github.com/lh3/seqtk> on this computer. "+
                                    "Please, try to (i) run it on a server/computer with larger amount of memory, or (ii) using command line option '--no-seqtk-subseq' !"))

                if not options.split_seqtk_subseq:
                    job.add('extract_short_reads.py',kind='program')
                    job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                    job.add('--input',input_file,kind='input',temp_path=temp_flag)
                    job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final2.txt'),kind='input')
                    job.add('--output','-',kind='parameter',checksum='no')
                else:
                    job.add('seqtk',kind='program')
                    job.add('subseq',kind='parameter')
                    job.add('',input_file,kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final2.txt'),kind='input')
                job.add('|',kind='parameter')
                job.add('seqtk',kind='parameter')
                job.add('trimfq',kind='parameter')
                job.add('-l','1',kind='parameter')
                #job.add('-q','0',kind='parameter')
                job.add('-B',options.trim_psl_3end_keep,kind='parameter')
                job.add('-',kind='parameter')
                job.add('>',outdir('reads_for-trimming.fq'),kind='output')
                job.run()


                job.add('cat',kind='program')
                job.add('',outdir('reads_not-for-trimming.fq'),kind='input',temp_path=temp_flag)
                job.add('',outdir('reads_for-trimming.fq'),kind='input',temp_path=temp_flag)
                job.add('>',output_file,kind='output')
                job.run()
            else:
                job.link(input_file, output_file, temp_path=temp_flag)

            if min_len_reads >= 47 and not options.skip_prefiltering_psl:
                # pre-filter the reads for BLAT because BLAT is slow
                # test to see if there are parts of the reads which can be mapped on gene-gene.fa using Bowtie which is faster
                # if yes then pass them to BLAT

                piece = 23
#                job.add('trim_reads.py',kind='program')
#                job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
#                job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_2.fq'),kind='output')
#                job.add('--trim_end','3',kind='parameter')
#                job.add('--final_size',piece,kind='parameter')
#                job.run()
                job.add('seqtk',kind='program')
                job.add('trimfq',kind='parameter')
                job.add('-l','1',kind='parameter')
                #job.add('-q','0',kind='parameter')
                job.add('-B',piece,kind='parameter')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
                job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_2.fq'),kind='output')
                job.run()


#                job.add('trim_reads.py',kind='program')
#                job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
#                job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_3.fq'),kind='output')
#                job.add('--trim_end','5',kind='parameter')
#                job.add('--final_size',piece,kind='parameter')
#                job.run()
                job.add('seqtk',kind='program')
                job.add('trimfq',kind='parameter')
                job.add('-l','1',kind='parameter')
                #job.add('-q','0',kind='parameter')
                job.add('-E',piece,kind='parameter')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
                job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_3.fq'),kind='output')
                job.run()


                # group reads which mapped partially
                #job.add('concatenate.py',kind='program')
                job.add('cat',kind='program')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_2.fq'),kind='input',temp_path=temp_flag)
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_3.fq'),kind='input',temp_path=temp_flag)
                job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.fq'),kind='output')
                job.run()

                if min_len_reads >= 59:
#                    job.add('trim_reads.py',kind='program')
#                    job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
#                    job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_4.fq'),kind='output')
#                    job.add('--trim_end','5',kind='parameter')
#                    job.add('--trim_size','7',kind='parameter')
#                    job.run()
                    job.add('seqtk',kind='program')
                    job.add('trimfq',kind='parameter')
                    job.add('-l','1',kind='parameter')
                    job.add('-b','7',kind='parameter')
                    job.add('-B',piece,kind='parameter')
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
                    #job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_4.fq'),kind='output')
                    #job.run()


#                    job.add('trim_reads.py',kind='program')
#                    job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_4.fq'),kind='input',temp_path=temp_flag)
#                    job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_5.fq'),kind='output')
#                    job.add('--trim_end','3',kind='parameter')
#                    job.add('--final_size',piece,kind='parameter')
#                    job.run()
                    #job.add('seqtk',kind='program')
                    #job.add('|',kind='parameter')
                    #job.add('seqtk',kind='parameter')
                    #job.add('trimfq',kind='parameter')
                    #job.add('-q','0',kind='parameter')
                    #job.add('-B',piece,kind='parameter')
                    #job.add('-',piece,kind='parameter')
                    #job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_4.fq'),kind='input')
                    job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_5.fq'),kind='output')
                    job.run()


#
                    job.add('seqtk',kind='program')
                    job.add('trimfq',kind='parameter')
                    job.add('-l','1',kind='parameter')
                    job.add('-e','7',kind='parameter')
                    job.add('-E',piece,kind='parameter')
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
                    #job.add('|',kind='parameter')
                    #job.add('seqtk',kind='parameter')
                    #job.add('trimfq',kind='parameter')
                    #job.add('-E',piece,kind='parameter')
                    #job.add('-',piece,kind='parameter')
                    job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_6.fq'),kind='output')
                    job.run()


                    # group reads which mapped partially
                    #job.add('concatenate.py',kind='program')
                    job.add('cat',kind='program')
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_5.fq'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_6.fq'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.fq'),kind='input',temp_path=temp_flag)
                    job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_x.fq'),kind='output')
                    job.run()
                else:
                    job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.fq'),
                             outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_x.fq'),
                             temp_path=temp_flag)

                # remove reads shorter than a given threshold
#                job.add('remove_shorter_reads.py',kind='program')
#                job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_x.fq'),kind='input',temp_path=temp_flag)
#                job.add('--threshold',piece,kind='parameter')
#                job.add('--output',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_final.fq'),kind='output')
#                job.run()
                job.add('seqtk',kind='program')
                job.add('seq',kind='parameter')
                job.add('-L',piece,kind='parameter')
                job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_x.fq'),kind='input',temp_path=temp_flag)
                job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_final.fq'),kind='output')
                job.run()

                if nucleotides_ggu > options.limit_bowtie:

                    job.add('split-fasta.py',kind='program')
                    job.add('--size',outdir('gene-gene_unique__nuc.txt'),kind='input')
                    job.add('--seqs',outdir('gene-gene_unique__seq.txt'),kind='input')
                    job.add('--threshold',options.limit_bowtie,kind='parameter')
                    job.add('-i',outdir('gene-gene_unique.fa'),kind='input')
                    job.add('-o',outdir('gene-gene_unique_split.fa'),kind='output')
                    job.run()

                    parts = [el.strip() for el in file(outdir('gene-gene_unique_split.fa'),'r').readlines()]
                    for i,part in enumerate(parts):
                        # build index
                        job.add('bowtie-build',kind='program')
                        job.add('-f',kind='parameter')
                        job.add('--quiet',kind='parameter')
#                        job.add('--ntoa',kind='parameter')
                        job.add('--offrate','1',kind='parameter')
                        job.add('--ftabchars','5',kind='parameter')
                        job.add('',part,kind='input',temp_path=temp_flag)
                        job.add('',part+'_dir/',kind='output')
                        job.run()

                        job.add('bowtie',kind='program')
                        job.add('-t',kind='parameter')
                        #job.add('-q',kind='parameter')
                        job.add('-v',options.mismatches,kind='parameter') #options.mismatches
                        job.add('-p',options.processes,kind='parameter',checksum='no')
                        job.add('-k','1',kind='parameter')
                        if os.path.isfile(os.path.join(part+'_dir','.1.ebwtl')):
                            job.add('--large-index',kind='parameter')
                        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                        job.add('--tryhard',kind='parameter')
                        job.add('--suppress','2,3,4,5,6,7,8',kind='parameter')
                        job.add('',part+'_dir/',kind='input',temp_path=temp_flag)
                        job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_final.fq'),kind='input',temp_path=temp_flag)
                        #job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.map'),kind='output')
                        job.add('2>',outdir('log_bowtie_psl_all.txt'),kind='output',checksum='no')
                        #job.add('2>&1',kind='parameter',checksum='no')
                        #job.run()
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('|',kind='parameter')
                        job.add('sort',kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        #job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.map'),kind='input',temp_path = temp_flag) # XXX
                        job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map.')+str(i),kind='output',dest_list='genegeneunique')
                        job.run()

                        job.clean(outdir('log_bowtie_psl_all.txt'),temp_path=temp_flag)

                    job.sink(job.genegeneunique, outdir('genegeneunique.txt'))

                    job.add('concatenate.py',kind='program')
                    job.add('-f',outdir('genegeneunique.txt'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map'),kind='output')
                    job.run()

#                    for ft in job.genegeneunique:
#                        job.clean(ft,temp_path=temp_flag)
                    job.clean(job.genegeneunique,temp_path=temp_flag)

                else:
                    # build index
                    job.add('bowtie-build',kind='program')
                    job.add('-f',kind='parameter')
                    job.add('--quiet',kind='parameter')
#                    job.add('--ntoa',kind='parameter')
                    job.add('--offrate','1',kind='parameter')
                    job.add('--ftabchars','5',kind='parameter')
                    job.add('',outdir('gene-gene_unique.fa'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('gene-gene_index/'),kind='output')
                    job.run()

                    job.add('bowtie',kind='program')
                    job.add('-t',kind='parameter')
                    #job.add('-q',kind='parameter')
                    job.add('-v',options.mismatches,kind='parameter') #options.mismatches
                    job.add('-p',options.processes,kind='parameter',checksum='no')
                    job.add('-k','1',kind='parameter')
                    if os.path.isfile(outdir('gene-gene_index','.1.ebwtl')):
                        job.add('--large-index',kind='parameter')
                    job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                    job.add('--tryhard',kind='parameter')
                    job.add('--suppress','2,3,4,5,6,7,8',kind='parameter')
                    job.add('',outdir('gene-gene_index/'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_final.fq'),kind='input',temp_path=temp_flag)
                    #job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.map'),kind='output')
                    job.add('2>',outdir('log_bowtie_psl_all.txt'),kind='output',checksum='no')
                    #job.add('2>&1',kind='parameter',checksum='no')
                    #job.run()
                    job.add('|',kind='parameter')
                    job.add('uniq',kind='parameter')
                    job.add('|',kind='parameter')
                    job.add('sort',kind='parameter')
                    if sort_buffer:
                        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                    if sort_parallel:
                        job.add('--parallel',options.processes,kind='parameter',checksum='no')
                    if sort_lzop_compress:
                        job.add('--compress-program','lzop',kind='parameter',checksum='no')
                    elif sort_gzip_compress:
                        job.add('--compress-program','gzip',kind='parameter',checksum='no')
                    job.add('-T',tmp_dir,kind='parameter',checksum='no')
                    job.add('|',kind='parameter')
                    job.add('uniq',kind='parameter')
                    #job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all.map'),kind='input',temp_path = temp_flag) # XXX
                    job.add('>',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map'),kind='output')
                    job.run()

                    job.clean(outdir('log_bowtie_psl_all.txt'),temp_path=temp_flag)

                # here I could remove the reads which have pieces mapping very closely to each other (all pieces [more than two] of a read map inside a gene => remove the read)

#                extract the short reads which mapped on gene-gene with Bowtie
#                job.add('extract_short_reads.py',kind='program')
#                job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input',temp_path=temp_flag)
#                job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map'),kind='input',temp_path=temp_flag)
#                job.add('--output',outdir('reads_gene-gene.fq'),kind='output')
#                job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
#                job.run(error_message = ("If this fails due to a memory error then lowering the "+
#                         "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
#                         "of FusionCatcher and running it again might help!"))

                if not options.split_seqtk_subseq:
                    job.add('extract_short_reads.py',kind='program')
                    job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                    job.add('--input',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
                    job.add('--list',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map'),kind='input')
                    job.add('--output','-',kind='parameter',checksum='no')
                else:
                    job.add('seqtk',kind='program')
                    job.add('subseq',kind='parameter')
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),kind='input')
                    job.add('',outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map'),kind='input')
#                job.add('|',kind='parameter')
#                job.add('seqtk',kind='parameter') # convert it to SANGER Qualities scores
#                job.add('seq',kind='parameter')
#                job.add('-Q64',kind='parameter')
#                job.add('-V',kind='parameter')
#                job.add('-',kind='parameter')
                job.add('>',outdir('reads_gene-gene.fq'),kind='output')
                job.run()

                job.clean(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),temp_path=temp_flag)
                job.clean(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl_all_uniq.map'),temp_path=temp_flag)

                job.add('printf',kind='program')
                job.add(('"\nReads Counts after BLAT/STAR/BOWTIE2/BWA prefiltering (and before alignment):\n'+
                         '--------------------------------------------------------------------------------\n"'),kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('cat',kind='program')
                job.add('',outdir('reads_gene-gene.fq'),kind='input')
                job.add('|',kind='parameter')
                job.add("echo $((`wc -l`/4))",kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('printf',kind='program')
                job.add('"\n\n\n"',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()

            else:
                # link
                job.link(outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),
                         outdir('reads_gene-gene.fq'),
                         temp_path = temp_flag)

            if not options.skip_str_filtering:
                # remove STR reads
                job.add('remove_str.py',kind='program')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--threshold','2.1',kind='parameter',checksum='no')
                job.add('--input',outdir('reads_gene-gene.fq'),kind='input',temp_path = temp_flag)
                job.add('--output',outdir('reads_gene-gene_no-str.fq'),kind='output')
                #job.add('--str',outdir('reads_gene-gene_str.fq'),kind='output',temp_path = temp_flag)
                job.add('--log',outdir('log_reads_gene-gene_no-str.txt'),kind='output')
                job.run()
                job.add('cat',kind='program')
                job.add('',outdir('log_reads_gene-gene_no-str.txt'),kind='input',temp_path = temp_flag)
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('printf',kind='program')
                job.add('"\n\n\n"',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
            else:
                job.link(outdir('reads_gene-gene.fq'),
                         outdir('reads_gene-gene_no-str.fq'),
                         temp_path = temp_flag)

            job.add('lengths_reads.py',kind='program')
            job.add('--input',outdir('reads_gene-gene_no-str.fq'),kind='input')
            job.add('--output',outdir('log_lengths_reads_gene-gene_no-str.txt'),kind='output')
            job.add('--counts',outdir('log_counts_reads_gene-gene_no-str.txt'),kind='output')
            job.run()

            # save lengths reads
            info(job,
                 fromfile = outdir('log_lengths_reads_gene-gene_no-str.txt'),
                 tofile = info_file,
                 top = ["Lengths of all reads before BLAT/STAR/BOWTIE2/BWA alignment",
                        "-----------------------------------------------------------"],
                 bottom = "\n\n\n")

            # save lengths reads
            info(job,
                 fromfile = outdir('log_counts_reads_gene-gene_no-str.txt'),
                 tofile = info_file,
                 top = ["Count of all reads before BLAT/STAR/BOWTIE2/BWA alignment",
                        "---------------------------------------------------------"],
                 bottom = "\n\n\n")

################################################################################
# BLAT alignment
################################################################################
            if not options.skip_blat:

                # convert FASTQ to FASTA
#                job.add('fastq2fasta.py',kind='program')
#                job.add('--input',outdir('reads_gene-gene_no-str.fq'),kind='input',temp_path = temp_flag if options.skip_star and options.skip_bowtie2 else 'no')
#                job.add('--output',outdir('reads_gene-gene.fa'),kind='output')
#                job.run()

                job.add('seqtk',kind='program')
                job.add('seq',kind='parameter')
                job.add('-a',kind='parameter')
                job.add('',outdir('reads_gene-gene_no-str.fq'),kind='input',temp_path = temp_flag if options.skip_star and options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('>',outdir('reads_gene-gene.fa'),kind='output')
                job.run()

                # find available memory
                job.add('printf',kind='program')
                job.add('"\n============\nMEMORY (before using BLAT):\n============\n"',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('free',kind='program')
                job.add('-m',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()

                if nucleotides_gg > options.limit_blat:

                    job.add('split-fasta.py',kind='program')
                    job.add('--size',outdir('gene-gene__nuc.txt'),kind='input')
                    job.add('--seqs',outdir('gene-gene__seq.txt'),kind='input')
                    job.add('--threshold',options.limit_blat,kind='parameter')
                    job.add('-i',outdir('gene-gene.fa'),kind='input')
                    job.add('-o',outdir('gene-gene_split_blat.fa'),kind='output')
                    job.add('-x',outdir('gene-gene_split_blat.len'),kind='output')
                    job.run()

                    parts = [el.strip() for el in file(outdir('gene-gene_split_blat.fa'),'r').readlines()]
                    maxlens = [el.strip() for el in file(outdir('gene-gene_split_blat.len'),'r').readlines()]
                    for i,part in enumerate(parts):
                        # file size
#                        job.add('du',kind='program')
#                        job.add('-b',part,kind='input')
#                        job.add('|',kind='parameter')
#                        job.add('cut',kind='parameter')
#                        job.add('-f','1',kind='parameter')
#                        job.add('>',part+'.len',kind='output')
#                        job.run()

                        # convert fasta to 2bit
                        job.add('faToTwoBit',kind='program')
                        job.add('',part,kind='input',temp_path=temp_flag)
                        job.add('',part+'.2bit',kind='output')
                        job.add('-noMask',kind='parameter')
                        job.run()

                        job.add('blat_parallel.py',kind='program')
                        job.add('-noHead',kind='parameter')
                        job.add('-stepSize=','5',kind='parameter',space='no') # 5
                        job.add('-tileSize=','11',kind='parameter',space='no') # 11
                        job.add('-minScore=','30',kind='parameter',space='no') # 30
                        job.add('-t=','DNA',kind='parameter',space='no')
                        job.add('-q=','RNA',kind='parameter',space='no') # orginall was DNA
                        #job.add('-fine',kind='parameter')
                        #job.add('-repMatch=','1000000',kind='parameter',space='no') # 2253 or 1000000?
                        job.add('-repMatch=','2253',kind='parameter',space='no') # 2253 or 1000000?
                        job.add('-minIdentity=','30',kind='parameter',space='no') # 0
                        job.add('-maxIntron=',maxlens[i],kind='parameter',space='no',from_file = 'yes') # default is 750000
                        #job.add('-maxIntron=',outdir('gene-gene_longest.txt'),kind='parameter',space='no',from_file = 'yes') # default is 750000
                        #job.add('-maxIntron=',part+'.len',kind='parameter',space='no',from_file = 'yes',temp_path=temp_flag) # default is 750000
                        job.add('--tmp_dir=',tmp_dir,kind='parameter',checksum='no',space='no')
                        job.add('--cpus=',options.processes,kind='parameter',checksum='no',space='no') # it takes 5 GB per cpu so for here it means 12 * 5 = 60GB
                        job.add('--filter-fusion',kind='parameter',checksum='no') # WARNING: this does a very fast pre-filtering as done in "find_fusion_genes_blat.py"; THIS should be kep in sync with "find_fusion_genes_blat.py"
                        job.add('',part+'.2bit',kind='input',temp_path=temp_flag)
                        job.add('',outdir('reads_gene-gene.fa'),kind='input')
                        job.add('',outdir('reads_blat_mapped_on_fusion_genes.psl.')+str(i),kind='output',dest_list='genegeneblat')
                        job.run()

                    #job.clean(outdir('gene-gene_split_blat.fa'),temp_path=temp_flag)
                    job.clean(outdir('reads_gene-gene.fa'),temp_path=temp_flag)

                    job.sink(job.genegeneblat, outdir('reads_blat_mapped_on_fusion_genes.psl.txt'))

                    job.add('concatenate.py',kind='program')
                    job.add('-f',outdir('reads_blat_mapped_on_fusion_genes.psl.txt'),kind='input',temp_path=temp_flag)
                    #job.add_list('',job.genegene,kind='input',temp_path=temp_flag,command_line='no')
                    job.add('',outdir('reads_blat_mapped_on_fusion_genes.psl'),kind='output')
                    job.run()

#                    for tfile in job.genegeneblat:
#                        job.clean(tfile,temp_path=temp_flag)
                    job.clean(job.genegeneblat,temp_path=temp_flag)
                    
                else: # no splitting of BLAT reference sequences
                    # convert fasta to 2bit
                    job.add('faToTwoBit',kind='program')
                    job.add('',outdir('gene-gene.fa'),kind='input')
                    job.add('',outdir('gene-gene.2bit'),kind='output')
                    job.add('-noMask',kind='parameter')
                    job.run()

                    # align the unmapped reads using BLAT on candidate fusion gene-gene
                    #
                    # web version of blat
                    #  blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 database.2bit query.fa output.psl
                    # from: http://http://genome.ucsc.edu/FAQ/FAQblat.html
                    #
                    # other idea: 	./blat -minIdentity=95 –fine -stepSize=1 –tileSize=6 -repMatch = 1000000
                    # from http://www.gene2drug.com/product/?p=671 by  Sucheta Tripathy
                    job.add('blat_parallel.py',kind='program')
                    job.add('-noHead',kind='parameter')
                    job.add('-stepSize=','5',kind='parameter',space='no') # 5
                    job.add('-tileSize=','11',kind='parameter',space='no') # 11
                    job.add('-minScore=','30',kind='parameter',space='no') # 30
                    job.add('-t=','DNA',kind='parameter',space='no')
                    job.add('-q=','RNA',kind='parameter',space='no') # orginall was DNA
                    #job.add('-fine',kind='parameter')
                    #job.add('-repMatch=','1000000',kind='parameter',space='no') # 2253 or 1000000?
                    job.add('-repMatch=','2253',kind='parameter',space='no') # 2253 or 1000000?
                    job.add('-minIdentity=','30',kind='parameter',space='no') # 0
                    job.add('-maxIntron=',outdir('gene-gene_longest.txt'),kind='parameter',space='no',from_file = 'yes') # default is 750000
                    job.add('--tmp_dir=',tmp_dir,kind='parameter',checksum='no',space='no')
                    job.add('--cpus=',options.processes,kind='parameter',checksum='no',space='no') # it takes 5 GB per cpu so for here it means 12 * 5 = 60GB
                    job.add('--filter-fusion',kind='parameter',checksum='no') # WARNING: this does a very fast pre-filtering as done in "find_fusion_genes_blat.py"; THIS should be kep in sync with "find_fusion_genes_blat.py"
                    job.add('',outdir('gene-gene.2bit'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_gene-gene.fa'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('reads_blat_mapped_on_fusion_genes.psl'),kind='output')
                    job.run()

                # find the best unique alignments of reads
                job.add('psl_best_unique_contigs.py',kind='program')
                job.add('--input',outdir('reads_blat_mapped_on_fusion_genes.psl'),kind='input',temp_path=temp_flag)
                job.add('--output',outdir('reads_best_unique_blat_mapped_on_fusion_genes.psl'),kind='output')
#                if (not empty(outdir('candidate_fusion-genes_further_mark.txt'))) and (not empty(datadir('custom_genes.txt'))):
#                    job.add('--ties',datadir('custom_genes_mark.txt'),kind='input')
                job.add('--ties-overlappings',datadir('ensembl_overlapping_genes.txt'),kind='input')
                job.add('--anchor',length_anchor_blat,kind='parameter') # find_fusion_genes_blat.py --threshold_overlap is enough!
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--tmp_dir',tmp_dir,kind='output',checksum='no')
                job.run()

                # more filtering -- remove the reads from the gene-gene junctions
                # which have the pair read mapping on a totally different gene than
                # those involved in the gene-gene junction
                if not options.all_reads_junction:
                    job.add('remove_reads_exon_exon_psl.py',kind='program')
                    job.add('--input_psl',outdir('reads_best_unique_blat_mapped_on_fusion_genes.psl'),kind='input',temp_path=temp_flag)
                    job.add('--input_transcriptome',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind='input',temp_path=temp_flag if options.skip_star and options.skip_bowtie2 and options.skip_bwa else 'no')
                    job.add('--output_psl',outdir('reads_best_unique_blat_mapped_on_fusion_genes_pairs.psl'),kind='output')
                    job.run()
                else:
                    job.link(outdir('reads_best_unique_blat_mapped_on_fusion_genes.psl'),
                             outdir('reads_best_unique_blat_mapped_on_fusion_genes_pairs.psl'),
                             temp_path=temp_flag)

                job.add('find_fusion_genes_psl.py',kind='program')
                job.add('--input_mappings',outdir('reads_best_unique_blat_mapped_on_fusion_genes_pairs.psl'),kind='input',temp_path=temp_flag)
                job.add('--input_genegene_fasta',outdir('gene-gene.fa'),kind='input',temp_path=temp_flag if options.skip_star and options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('--input_hugo',datadir('genes_symbols.txt'),kind='input')
                job.add('--input_genes_positions',datadir('genes.txt'),kind='input')
                job.add('--threshold_overlap',length_anchor_blat,kind='parameter')
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--output',outdir('candidates_fusion_genes_reads_blat.txt'),kind='output')
                job.run()

                # summary the gene-gene mappings
                job.add('build_report_fusions_psl.py',kind='program')
                job.add('--suporting_unique_reads',spanning_reads_blat,kind='parameter')
                job.add('--anchor2',length_anchor2,kind='parameter')
                job.add('--input_candidate_fusion_genes_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='input',temp_path=temp_flag if options.skip_star and options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('--input_fastq',outdir('original_important.fq.gz'),kind='input',temp_path=temp_flag if options.skip_star and options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('--input_fusion_psl',outdir('candidates_fusion_genes_reads_blat.txt'),kind='input',temp_path=temp_flag)
                job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
                if options.psl_visualization and not empty(datadir('genome.2bit')):
                    job.add('--input_genome_2bit',datadir('genome.2bit'),kind='input')
                    job.add('--psl_alignment_type','web',kind='parameter')
                if options.sam_visualization:
                    job.add('--input_genome_bowtie2',datadir('genome_index2/index'),kind='input')
                    job.add('--sam_alignment','20',kind='parameter')
                    job.add('--threads',options.processes,kind='parameter')
                if options.assembly:
                    job.add('--velvet',kind='parameter')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BLAT.txt'),kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BLAT.zip'),kind='output')
                job.run()

################################################################################
# STAR alignment
################################################################################
            if not options.skip_star:

                # STAR is removing the /1 and /2 from the end of the reads names
                # changing "/1" and "/2" into "-1" "-2" such that STAR does not remove them
                job.add('sed',kind='program')
                job.add("""'s/\/\([1-2]$\)/\-\\1/;n;n;n'""",outdir('reads_gene-gene_no-str.fq'),kind='input')
                job.add('>',outdir('reads_gene-gene_no-str_fixed.fq'),kind='output')
                job.run()


                if job.run():
                    sdjboverhang = int(file(outdir('log_lengths_reads_gene-gene_no-str.txt'),'r').readline().strip()) - 1
                    file(outdir('star_sdjboverhang.txt'),'w').write("%d" % (sdjboverhang,))
                genomesaindexnbases = int(min(14, math.log(nucleotides_gg,2)/(float(2) - 1)))
                genomechrbinnbits = int(min(18, math.log(float(nucleotides_gg)/float(sequences_gg),2)))

                # find available memory
                job.add('printf',kind='program')
                job.add('"\n============\nMEMORY (before using STAR):\n============\n"',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('free',kind='program')
                job.add('-m',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()

                if nucleotides_gg > options.limit_star:

                    job.add('split-fasta.py',kind='program')
                    job.add('--size',outdir('gene-gene__nuc.txt'),kind='input')
                    job.add('--seqs',outdir('gene-gene__seq.txt'),kind='input')
                    job.add('--threshold',options.limit_star,kind='parameter')
                    job.add('-i',outdir('gene-gene.fa'),kind='input')
                    job.add('-o',outdir('gene-gene_split_star.fa'),kind='output')
                    job.add('-x',outdir('gene-gene_split_star.len'),kind='output')
                    job.run()

                    parts = [el.strip() for el in file(outdir('gene-gene_split_star.fa'),'r').readlines()]
                    maxlens = [el.strip() for el in file(outdir('gene-gene_split_star.len'),'r').readlines()]
                    for i,part in enumerate(parts):

                        # get the length of the FASTA file
                        job.add('du',kind='program')
                        job.add('-b',part,kind='input')
                        job.add('|',kind='parameter')
                        job.add('cut',kind='parameter')
                        job.add('-f','1',kind='parameter')
                        job.add('>',part+'.len',kind='output')
                        job.run()


                        t = "[from file: '%s']" % (maxlens[i],)
                        if job.run():
                            t = file(maxlens[i],'r').readline().strip()
                            
                        genomesaindexnbases = int(min(14, math.log(100,2)/(float(2) - 1)))
                        if job.run():
                            ti = file(part+'.len','r').readline().strip()
                            lenparts = len(parts)
                            genomesaindexnbases = int(min(14, math.log(float(ti),2)/(float(2) - 1)))
                            #genomechrbinnbits = int(min(18, math.log(float(t)/(math.ceil(float(sequences_gg)/float(lenparts)))+2,2)))
                        job.clean(part+'.len',temp_path=temp_flag)



                        # build the STAR index
                        gd = "%s_star/" % (part,)
                        gdr = "%s_star-results/" % (part,)
                        job.add('STAR',kind='program')
                        job.add('--genomeChrBinNbits',genomechrbinnbits,kind='parameter')
                        job.add('--genomeSAindexNbases',genomesaindexnbases,kind='parameter')
                        job.add('--runMode','genomeGenerate',kind='parameter')
                        job.add('--runThreadN',options.processes,kind='parameter',checksum='no')
                        job.add('--genomeDir',gd,kind='output')
                        job.add('--genomeFastaFiles',part,kind='input')
                        job.add('--outFileNamePrefix',gd,kind='output')
                        job.run()

                        # align the unmapped reads using STAR on candidate fusion gene-gene
                        job.add('STAR',kind='program')
                        job.add('--twopass1readsN',outdir('log_counts_reads_gene-gene_no-str.txt'),kind='parameter',from_file='yes')
                        job.add('--twopassMode','Basic',kind='parameter')
                        job.add('--genomeSAindexNbases',genomesaindexnbases,kind='parameter')
                        job.add('--sjdbOverhang',outdir('star_sdjboverhang.txt'),kind='parameter',from_file='yes')
                        #job.add('--alignIntronMax',outdir('gene-gene_longest.txt'),kind='parameter',from_file = 'yes')
                        job.add('--alignIntronMax',maxlens[i],kind='parameter',from_file = 'yes')
                        if options.skip_star_bowtie:
                            job.add('--outFilterMatchNmin',int(float(min_len_reads)*0.90),kind='parameter')
                            job.add('--outFilterMatchNminOverLread','0.90',kind='parameter')
                            job.add('--outFilterScoreMinOverLread','0.90',kind='parameter')  # NEW in v0.99.4b
                            job.add('--alignSplicedMateMapLminOverLmate','0.90',kind='parameter') # NEW in v0.99.4b
                        else:
                            job.add('--outFilterMatchNmin',int(float(min_len_reads)*0.49),kind='parameter')
                            job.add('--outFilterMatchNminOverLread','0.49',kind='parameter')
                            job.add('--outFilterScoreMinOverLread','0.49',kind='parameter')  # NEW in v0.99.4b
                            job.add('--alignSplicedMateMapLminOverLmate','0.49',kind='parameter') # NEW in v0.99.4b
                        job.add('--genomeDir',gd,kind='input',temp_path=temp_flag)
                        job.add('--runThreadN',options.processes,kind='parameter',checksum='no')
                        job.add('--seedSearchStartLmax',length_anchor_star,kind='parameter') # 20 # default is: 50
                        job.add('--alignSJoverhangMin',length_anchor_star-1,kind='parameter') # 9 # default is 5? # NEW in v0.99.4b
                        job.add('--outSJfilterOverhangMin','10 10 10 10',kind='parameter')# default is: 30 12 12 12 ("non-canonical motifs","GT/AG"motif,"GC/AG"motif,"AT/AC"motif)
                        job.add('--outSJfilterCountUniqueMin','1 1 1 1',kind='parameter')# default is: 3 1 1 1
                        job.add('--outSJfilterCountTotalMin','1 1 1 1',kind='parameter')# default is: 3 1 1 1
                        job.add('--outSJfilterDistToOtherSJmin','0 0 0 0',kind='parameter')# default is: 10 0 5 10
                        job.add('--outSJfilterIntronMaxVsReadN','%s %s %s' % (t,t,t),kind='parameter')# default is: 50000 100000 200000
                        job.add('--limitOutSAMoneReadBytes','100000000',kind='parameter')
                        job.add('--scoreGapNoncan','-4',kind='parameter')
                        job.add('--scoreGapATAC','-4',kind='parameter')
                        job.add('--readFilesIn',outdir('reads_gene-gene_no-str_fixed.fq'),kind='input')
                        job.add('--outFileNamePrefix',gdr,kind='output')
                        job.add('--outFileNamePrefix',os.path.join(gdr,'Aligned.out.sam'),kind='output',command_line = 'no')
                        job.run()

                        job.add('sed',kind='program')
                        job.add("""'s/\-\([1-2]\\t\)/\/\\1/'""",os.path.join(gdr,'Aligned.out.sam'),kind='input')
                        job.add("",gdr,kind='input',temp_path=temp_flag,command_line='no')
                        job.add('|',kind='parameter')
                        job.add('sam2psl.py',kind='parameter')
                        job.add('--input','-',kind='parameter')
                        job.add('--output','-',kind='output')
                        #job.add('--output',outdir('gene-gene-star.psl.')+str(i),kind='output',dest_list='genegenestar')
                        #job.run()
                        job.add('|',kind='parameter')
                        job.add('LC_ALL=C',kind='parameter')
                        job.add('sort',kind='parameter')
                        job.add('-k','10,10',kind='parameter')
                        job.add('-k','14,14',kind='parameter')
                        job.add('-k','12,12n',kind='parameter')
                        job.add('-k','13,13n',kind='parameter')
                        job.add('-t',"'\t'",kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('>',outdir('gene-gene-star.psl.')+str(i),kind='output')
                        job.run()

                        if options.skip_star_bowtie:
                            job.clean(part,temp_path=temp_flag)
                            job.link(outdir('gene-gene-star.psl.')+str(i),
                                     outdir('gene-gene-star_more.psl.')+str(i),
                                     temp_path=temp_flag,
                                     dest_list='genegenestar')
                        else:
                            job.add('analyze_splits_sam.py',kind='program')
                            job.add('--input',outdir('gene-gene-star.psl.')+str(i),kind='input')
                            job.add('--output',outdir('gene-gene-star_final.psl.')+str(i),kind='output',temp_path=temp_flag)
                            job.add('--clipped-reads-ids',outdir('reads-ids_clip_psl_star.txt.')+str(i),kind='output')
                            job.add('--clipped-reads-refs',outdir('reads-refs_clip_psl_star.txt.')+str(i),kind='output')
                            job.add('--clip-min',length_anchor_star,kind='parameter')
                            job.run()

                            if job.iff(empty(outdir('reads-ids_clip_psl_star.txt.')+str(i)),id = "#reads-ids-clip-psl-star."+str(i)+"#"):
                                job.clean(outdir('reads-ids_clip_psl_star.txt.')+str(i),temp_path=temp_flag)
                                job.clean(outdir('reads-refs_clip_psl_star.txt.')+str(i),temp_path=temp_flag)
                                job.clean(part,temp_path=temp_flag)
                            else:
                                job.add('LC_ALL=C',kind='program')
                                job.add('sort',kind='parameter')
                                if sort_buffer:
                                    job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                                if sort_parallel:
                                    job.add('--parallel',options.processes,kind='parameter',checksum='no')
                                if sort_lzop_compress:
                                    job.add('--compress-program','lzop',kind='parameter',checksum='no')
                                elif sort_gzip_compress:
                                    job.add('--compress-program','gzip',kind='parameter',checksum='no')
                                job.add('-T',tmp_dir,kind='parameter',checksum='no')
                                job.add('',outdir('reads-ids_clip_psl_star.txt.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('|',kind='parameter')
                                job.add('uniq',kind='parameter')
                                job.add('>',outdir('reads-ids_clip_star_psl_uniq.txt.')+str(i),kind='output')
                                job.run()

                                job.add('cut',kind='program')
                                job.add('-f1',outdir('reads-ids_clip_star_psl_uniq.txt.')+str(i),kind='input')
                                job.add('|',kind='parameter')
                                job.add('uniq',kind='parameter')
                                job.add('|',kind='parameter')
                                job.add('LC_ALL=C',kind='parameter')
                                job.add('sort',kind='parameter')
                                if sort_buffer:
                                    job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                                if sort_parallel:
                                    job.add('--parallel',options.processes,kind='parameter',checksum='no')
                                if sort_lzop_compress:
                                    job.add('--compress-program','lzop',kind='parameter',checksum='no')
                                elif sort_gzip_compress:
                                    job.add('--compress-program','gzip',kind='parameter',checksum='no')
                                job.add('-T',tmp_dir,kind='parameter',checksum='no')
                                job.add('|',kind='parameter')
                                job.add('uniq',kind='parameter')
                                job.add('|',kind='parameter')
                                job.add('seqtk',kind='parameter')
                                job.add('subseq',kind='parameter')
                                job.add('',outdir('reads_gene-gene_no-str.fq'),kind='input')
                                job.add('-',kind='parameter')
                                job.add('>',outdir('reads-ids_clip_star_psl.fq.')+str(i),kind='output')
                                job.run()

                                job.add('split-reads.py',kind='program')
                                job.add('--input',outdir('reads-ids_clip_star_psl.fq.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('--list',outdir('reads-ids_clip_star_psl_uniq.txt.')+str(i),kind='input',temp_path=temp_flag)
                                #job.add('--wiggle-size','1',kind='parameter')
                                job.add('--output-1',outdir('reads-ids_clip_star_psl_r1.fq.')+str(i),kind='output')
                                job.add('--output-2',outdir('reads-ids_clip_star_psl_r2.fq.')+str(i),kind='output')
                                job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                                job.run(error_message = ("If this fails due to a memory error then lowering the "+
                                                         "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
                                                         "of FusionCatcher and running it again might help!"))


                                job.add('LC_ALL=C',kind='program')
                                job.add('sort',kind='parameter')
                                if sort_buffer:
                                    job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                                if sort_parallel:
                                    job.add('--parallel',options.processes,kind='parameter',checksum='no')
                                if sort_lzop_compress:
                                    job.add('--compress-program','lzop',kind='parameter',checksum='no')
                                elif sort_gzip_compress:
                                    job.add('--compress-program','gzip',kind='parameter',checksum='no')
                                job.add('-T',tmp_dir,kind='parameter',checksum='no')
                                job.add('',outdir('reads-refs_clip_psl_star.txt.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('|',kind='parameter')
                                job.add('uniq',kind='parameter')
                                job.add('>',outdir('reads-refs_clip_star_psl_uniq.txt.')+str(i),kind='output')
                                job.run()

                                gda = "%s_bowtie_star.fa" % (part,)
                                job.add('seqtk',kind='program')
                                job.add('subseq',kind='parameter')
                                job.add('',part,kind='input',temp_path=temp_flag)
                                job.add('',outdir('reads-refs_clip_star_psl_uniq.txt.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('>',gda,kind='output')
                                job.run()

                                gdb = "%s_bowtie_star/" % (part,)
                                job.add('bowtie-build',kind='program')
                                job.add('-f',kind='parameter')
                                job.add('--quiet',kind='parameter')
    #                            job.add('--ntoa',kind='parameter')
                                job.add('--offrate','1',kind='parameter')
                                job.add('--ftabchars','5',kind='parameter')
                                #job.add('',part,kind='input',temp_path=temp_flag)
                                job.add('',gda,kind='input',temp_path=temp_flag)
                                job.add('',gdb,kind='output',checksum='no')
                                job.add('',gdb,kind='output',command_line='no')
                                job.run()

                                # map using bowtie
                                job.add('bowtie',kind='program')
                                job.add('-t',kind='parameter')
                                #job.add('-q',kind='parameter')
                                job.add('-a',kind='parameter')
                                job.add('-v',options.mismatches,kind='parameter')
                                job.add('-p',options.processes,kind='parameter',checksum='no')
                                if os.path.isfile(os.path.join(gdb,'.1.ebwtl')):
                                    job.add('--large-index',kind='parameter')
                                job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                                job.add('--tryhard',kind='parameter')
                                job.add('--best',kind='parameter')
                                job.add('--strata',kind='parameter')
                                job.add('--sam',kind='parameter')
                                job.add('--ff',kind='parameter')
                                #job.add('-X',outdir('gene-gene_longest.txt'),kind='parameter',from_file="yes")
                                job.add('-X',maxlens[i],kind='parameter',from_file="yes")
                                job.add('',gdb,kind='input',temp_path=temp_flag)
                                job.add('-1',outdir('reads-ids_clip_star_psl_r1.fq.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('-2',outdir('reads-ids_clip_star_psl_r2.fq.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('2>',outdir('log_bowtie_reads_mapped-gene-gene-star.stdout.txt.')+str(i),kind='output',checksum='no',temp_path=temp_flag)
                                job.add('|',kind='parameter')
                                job.add('awk',"""'$3 == "*" { next } { print }'""",kind='parameter')
                                job.add('>',outdir('split_gene-gene_star.sam.')+str(i),kind='output')
                                job.run()

                                job.add('merge-sam.py',kind='program')
                                job.add('--input',outdir('split_gene-gene_star.sam.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('--output',outdir('split_gene-gene_star_patch.sam.')+str(i),kind='output')
                                job.run()

                                job.add('sam2psl.py',kind='program')
                                job.add('--input',outdir('split_gene-gene_star_patch.sam.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('--output',outdir('split_gene-gene_star_patch.psl.')+str(i),kind='output')
                                job.run()

                                job.add('analyze_splits_sam.py',kind='program')
                                job.add('--input',outdir('split_gene-gene_star_patch.psl.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('--output',outdir('split_gene-gene_star_final.psl.')+str(i),kind='output')
                                job.add('--remove-extra',kind='parameter')
                                job.run()

                            if job.iff(empty(outdir('split_gene-gene_star_final.psl.')+str(i)),id = "#split_gene-gene_star_final."+str(i)+"#"):
                                job.link(outdir('gene-gene-star.psl.')+str(i),
                                         outdir('gene-gene-star_more.psl.')+str(i),
                                         temp_path=temp_flag,
                                         dest_list='genegenestar')
                                job.clean(outdir('split_gene-gene_star_final.psl.')+str(i),temp_path=temp_flag)
                            else:
                                job.add('cat',kind='program')
                                job.add('',outdir('split_gene-gene_star_final.psl.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('',outdir('gene-gene-star.psl.')+str(i),kind='input',temp_path=temp_flag)
                                job.add('>',outdir('gene-gene-star_more.psl.')+str(i),kind='output',dest_list='genegenestar')
                                job.run()


                    #job.clean(outdir('gene-gene_split_star.fa'),temp_path=temp_flag)
                    job.clean(outdir('reads_gene-gene_no-str_fixed.fq'),temp_path=temp_flag if options.skip_bwa else 'no')
                    job.clean(outdir('reads_gene-gene_no-str.fq'),temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                    job.sink(job.genegenestar, outdir('gene-gene-star_more.psl.txt'))

                    job.add('concatenate.py',kind='program')
                    job.add('-f',outdir('gene-gene-star_more.psl.txt'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('gene-gene-star_more.psl'),kind='output')
                    job.run()

                    #for tfile in job.genegenestar:
                    #    job.clean(tfile,temp_path=temp_flag)
                    job.clean(job.genegenestar,temp_path=temp_flag)
                    
                else:
                    # build the STAR index
                    job.add('STAR',kind='program')
                    job.add('--genomeChrBinNbits',genomechrbinnbits,kind='parameter')
                    job.add('--genomeSAindexNbases',genomesaindexnbases,kind='parameter')
                    job.add('--runMode','genomeGenerate',kind='parameter')
                    job.add('--runThreadN',options.processes,kind='parameter',checksum='no')
                    job.add('--genomeDir',outdir('gene-gene-star/'),kind='output')
                    job.add('--genomeFastaFiles',outdir('gene-gene.fa'),kind='input')
                    job.add('--outFileNamePrefix',outdir('gene-gene-star-results/'),kind='output')
                    job.run()

                    # align the unmapped reads using STAR on candidate fusion gene-gene
                    t = "[from file: '%s']" % (outdir('gene-gene_longest.txt'),)
                    if job.run():
                        t = file(outdir('gene-gene_longest.txt'),'r').readline().strip()
                    job.add('STAR',kind='program')
                    job.add('--twopass1readsN',outdir('log_counts_reads_gene-gene_no-str.txt'),kind='parameter',from_file='yes')
                    job.add('--twopassMode','Basic',kind='parameter')
                    job.add('--genomeSAindexNbases',genomesaindexnbases,kind='parameter')
                    job.add('--sjdbOverhang',outdir('star_sdjboverhang.txt'),kind='parameter',from_file='yes')
                    job.add('--alignIntronMax',outdir('gene-gene_longest.txt'),kind='parameter',from_file = 'yes')
                    if options.skip_star_bowtie:
                        job.add('--outFilterMatchNmin',int(float(min_len_reads)*0.90),kind='parameter')
                        job.add('--outFilterMatchNminOverLread','0.90',kind='parameter')
                        job.add('--outFilterScoreMinOverLread','0.90',kind='parameter')  # NEW in v0.99.4b
                        job.add('--alignSplicedMateMapLminOverLmate','0.90',kind='parameter') # NEW in v0.99.4b
                    else:
                        job.add('--outFilterMatchNmin',int(float(min_len_reads)*0.49),kind='parameter')
                        job.add('--outFilterMatchNminOverLread','0.49',kind='parameter')
                        job.add('--outFilterScoreMinOverLread','0.49',kind='parameter')  # NEW in v0.99.4b
                        job.add('--alignSplicedMateMapLminOverLmate','0.49',kind='parameter') # NEW in v0.99.4b
                    job.add('--genomeDir',outdir('gene-gene-star/'),kind='input',temp_path=temp_flag)
                    job.add('--runThreadN',options.processes,kind='parameter',checksum='no')
                    job.add('--seedSearchStartLmax',length_anchor_star,kind='parameter')# default is: 50
                    job.add('--alignSJoverhangMin',length_anchor_star-1,kind='parameter') #9 # default is 5? # NEW in v0.99.4b
                    job.add('--outSJfilterOverhangMin','10 10 10 10',kind='parameter')# default is: 30 12 12 12 ("non-canonical motifs","GT/AG"motif,"GC/AG"motif,"AT/AC"motif)
                    job.add('--outSJfilterCountUniqueMin','1 1 1 1',kind='parameter')# default is: 3 1 1 1
                    job.add('--outSJfilterCountTotalMin','1 1 1 1',kind='parameter')# default is: 3 1 1 1
                    job.add('--outSJfilterDistToOtherSJmin','0 0 0 0',kind='parameter')# default is: 10 0 5 10
                    job.add('--outSJfilterIntronMaxVsReadN','%s %s %s' % (t,t,t),kind='parameter')# default is: 50000 100000 200000
                    job.add('--limitOutSAMoneReadBytes','100000000',kind='parameter')
                    job.add('--scoreGapNoncan','-4',kind='parameter')
                    job.add('--scoreGapATAC','-4',kind='parameter')
#                    job.add('--outFilterMultimapScoreRange','10',kind='parameter')
#                    job.add('--outFilterMultimapNmax','10000',kind='parameter')
#                    job.add('--chimScoreJunctionNonGTAG','0',kind='parameter')
#                    job.add('--chimScoreDropMax','10000',kind='parameter')# default is: 0
#                    job.add('--chimScoreMin','0',kind='parameter')# default is: 0
#                    job.add('--chimScoreSeparation','10',kind='parameter')# default is: 0
#                    job.add('--chimSegmentMin',outdir('gene-gene_longest.txt'),kind='parameter',from_file = 'yes')
#                    job.add('--chimJunctionOverhangMin',outdir('gene-gene_longest.txt'),kind='parameter',from_file = 'yes')
                    job.add('--readFilesIn',outdir('reads_gene-gene_no-str_fixed.fq'),kind='input',temp_path=temp_flag if options.skip_bwa else 'no')
                    job.add('--outFileNamePrefix',outdir('gene-gene-star-results/'),kind='output')
                    job.run()

# this works!
#
# STAR \
# --alignIntronMax [ARGUMENT from file '/apps/fusioncatcher/output/snu16/gene-gene_longest.txt'] \
# --outFilterMatchNmin 30 \
# --outFilterMatchNminOverLread 0.90 \
# --genomeDir /apps/fusioncatcher/output/snu16/gene-gene-star/ \
# --runThreadN 4 \
# --seedSearchStartLmax 20 \
# --alignSJoverhangMin 9 \
# --outSJfilterOverhangMin 10 10 10 10 \
# --outSJfilterCountUniqueMin 1 1 1 1 \
# --outSJfilterCountTotalMin 1 1 1 1 \
# --outSJfilterDistToOtherSJmin 0 0 0 0 \
# --outSJfilterIntronMaxVsReadN [from file: '/apps/fusioncatcher/output/snu16/gene-gene_longest.txt'] [from file: '/apps/fusioncatcher/output/snu16/gene-gene_longest.txt'] [from file: '/apps/fusioncatcher/output/snu16/gene-gene_longest.txt'] \
# --readFilesIn /apps/fusioncatcher/output/snu16/reads_gene-gene_no-str_fixed.fq \
# --outFileNamePrefix /apps/fusioncatcher/output/snu16/gene-gene-star-results/


                    job.add('sed',kind='program')
                    job.add("""'s/\-\([1-2]\\t\)/\/\\1/'""",outdir('gene-gene-star-results','Aligned.out.sam'),kind='input')
                    job.add("",outdir('gene-gene-star-results/'),kind='input',temp_path=temp_flag,command_line='no')
                    job.add('|',kind='parameter')
                    job.add('sam2psl.py',kind='parameter')
                    job.add('--input','-',kind='parameter')
                    #job.add('--output',outdir('gene-gene-star.psl'),kind='output')
                    #job.run()
                    job.add('--output','-',kind='output')
                    job.add('|',kind='parameter')
                    job.add('LC_ALL=C',kind='parameter')
                    job.add('sort',kind='parameter')
                    job.add('-k','10,10',kind='parameter')
                    job.add('-k','14,14',kind='parameter')
                    job.add('-k','12,12n',kind='parameter')
                    job.add('-k','13,13n',kind='parameter')
                    job.add('-t',"'\t'",kind='parameter')
                    if sort_buffer:
                        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                    if sort_parallel:
                        job.add('--parallel',options.processes,kind='parameter',checksum='no')
                    if sort_lzop_compress:
                        job.add('--compress-program','lzop',kind='parameter',checksum='no')
                    elif sort_gzip_compress:
                        job.add('--compress-program','gzip',kind='parameter',checksum='no')
                    job.add('-T',tmp_dir,kind='parameter',checksum='no')
                    job.add('>',outdir('gene-gene-star.psl'),kind='output')
                    job.run()


                    if options.skip_star_bowtie:
                        job.link(outdir('gene-gene-star.psl'),
                                 outdir('gene-gene-star_more.psl'),
                                 temp_path=temp_flag)
                        job.clean(outdir('reads_gene-gene_no-str.fq'),temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                    else:
                        job.add('analyze_splits_sam.py',kind='program')
                        job.add('--input',outdir('gene-gene-star.psl'),kind='input')
                        job.add('--output',outdir('gene-gene-star_final.psl'),kind='output',temp_path=temp_flag)
                        job.add('--clipped-reads-ids',outdir('reads-ids_clip_psl_star.txt'),kind='output')
                        job.add('--clipped-reads-refs',outdir('reads-refs_clip_psl_star.txt'),kind='output')
                        job.add('--clip-min',length_anchor_star,kind='parameter')
                        job.run()



    #                    job.add('sed',kind='program')
    #                    job.add("""'s/\-\([1-2]\\t\)/\/\\1/'""",outdir('gene-gene-star-results','Chimeric.out.sam'),kind='input')
    #                    job.add("",outdir('gene-gene-star-results/'),kind='input',temp_path=temp_flag,command_line='no')
    #                    job.add('|',kind='parameter')
    #                    job.add('sam2psl.py',kind='parameter')
    #                    job.add('--input','-',kind='parameter')
    #                    job.add('--output',outdir('gene-gene-star-chimeric.psl'),kind='output')
    #                    job.run()
    #
    #                    job.add('analyze_star_chimeric.py',kind='program')
    #                    job.add('--input',outdir('gene-gene-star-chimeric.psl'),kind='input',temp_path=temp_flag)
    #                    job.add('--output',outdir('gene-gene-star-chimeric_final.psl'),kind='output')
    #                    job.run()
    #
    #                    # group reads which map on transcriptome in one FASTQ file
    #                    job.add('concatenate.py',kind='program')
    #                    job.add('',outdir('gene-gene-star.psl'),kind='input',temp_path=temp_flag)
    #                    job.add('',outdir('gene-gene-star-chimeric_final.psl'),kind='input',temp_path=temp_flag)
    #                    job.add('',outdir('gene-gene-star_all.psl'),kind='output')
    #                    job.run()


                        if job.iff(empty(outdir('reads-ids_clip_psl_star.txt')),id = "#reads-ids-clip-psl-star#"):
                            job.clean(outdir('reads-ids_clip_psl_star.txt'),temp_path=temp_flag)
                            job.clean(outdir('reads-refs_clip_psl_star.txt'),temp_path=temp_flag)
                            job.clean(outdir('reads_gene-gene_no-str.fq'),temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                        else:
                            job.add('LC_ALL=C',kind='program')
                            job.add('sort',kind='parameter')
                            if sort_buffer:
                                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                            if sort_parallel:
                                job.add('--parallel',options.processes,kind='parameter',checksum='no')
                            if sort_lzop_compress:
                                job.add('--compress-program','lzop',kind='parameter',checksum='no')
                            elif sort_gzip_compress:
                                job.add('--compress-program','gzip',kind='parameter',checksum='no')
                            job.add('-T',tmp_dir,kind='parameter',checksum='no')
                            job.add('',outdir('reads-ids_clip_psl_star.txt'),kind='input',temp_path=temp_flag)
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('>',outdir('reads-ids_clip_star_psl_uniq.txt'),kind='output')
                            job.run()

                            job.add('cut',kind='program')
                            job.add('-f1',outdir('reads-ids_clip_star_psl_uniq.txt'),kind='input')
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('|',kind='parameter')
                            job.add('LC_ALL=C',kind='parameter')
                            job.add('sort',kind='parameter')
                            if sort_buffer:
                                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                            if sort_parallel:
                                job.add('--parallel',options.processes,kind='parameter',checksum='no')
                            if sort_lzop_compress:
                                job.add('--compress-program','lzop',kind='parameter',checksum='no')
                            elif sort_gzip_compress:
                                job.add('--compress-program','gzip',kind='parameter',checksum='no')
                            job.add('-T',tmp_dir,kind='parameter',checksum='no')
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('|',kind='parameter')
                            job.add('seqtk',kind='parameter')
                            job.add('subseq',kind='parameter')
                            job.add('',outdir('reads_gene-gene_no-str.fq'),kind='input',temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                            job.add('-',kind='parameter')
                            job.add('>',outdir('reads-ids_clip_star_psl.fq'),kind='output')
                            job.run()

                            job.add('split-reads.py',kind='program')
                            job.add('--input',outdir('reads-ids_clip_star_psl.fq'),kind='input',temp_path=temp_flag)
                            job.add('--list',outdir('reads-ids_clip_star_psl_uniq.txt'),kind='input',temp_path=temp_flag)
                            job.add('--output-1',outdir('reads-ids_clip_star_psl_r1.fq'),kind='output')
                            job.add('--output-2',outdir('reads-ids_clip_star_psl_r2.fq'),kind='output')
                            #job.add('--wiggle-size','1',kind='parameter')
                            job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                            job.run(error_message = ("If this fails due to a memory error then lowering the "+
                                                     "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
                                                     "of FusionCatcher and running it again might help!"))


                            job.add('LC_ALL=C',kind='program')
                            job.add('sort',kind='parameter')
                            if sort_buffer:
                                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                            if sort_parallel:
                                job.add('--parallel',options.processes,kind='parameter',checksum='no')
                            if sort_lzop_compress:
                                job.add('--compress-program','lzop',kind='parameter',checksum='no')
                            elif sort_gzip_compress:
                                job.add('--compress-program','gzip',kind='parameter',checksum='no')
                            job.add('-T',tmp_dir,kind='parameter',checksum='no')
                            job.add('',outdir('reads-refs_clip_psl_star.txt'),kind='input',temp_path=temp_flag)
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('>',outdir('reads-refs_clip_star_psl_uniq.txt'),kind='output')
                            job.run()

                            job.add('seqtk',kind='program')
                            job.add('subseq',kind='parameter')
                            job.add('',outdir('gene-gene.fa'),kind='input')
                            job.add('',outdir('reads-refs_clip_star_psl_uniq.txt'),kind='input',temp_path=temp_flag)
                            job.add('>',outdir('gene-gene-bowtie_star.fa'),kind='output')
                            job.run()

                            job.add('bowtie-build',kind='program')
                            job.add('-f',kind='parameter')
                            job.add('--quiet',kind='parameter')
    #                        job.add('--ntoa',kind='parameter')
                            job.add('--offrate','1',kind='parameter')
                            job.add('--ftabchars','5',kind='parameter')
                            #job.add('',outdir('gene-gene.fa'),kind='input')
                            job.add('',outdir('gene-gene-bowtie_star.fa'),kind='input',temp_path=temp_flag)
                            job.add('',outdir('gene-gene-bowtie_star/'),kind='output',checksum='no')
                            job.add('',outdir('gene-gene-bowtie_star/'),kind='output',command_line='no')
                            job.run()

                            # map using bowtie
                            job.add('bowtie',kind='program')
                            job.add('-t',kind='parameter')
                            #job.add('-q',kind='parameter')
                            job.add('-a',kind='parameter')
                            job.add('-v',options.mismatches,kind='parameter')
                            job.add('-p',options.processes,kind='parameter',checksum='no')
                            if os.path.isfile(os.path.join(outdir('gene-gene-bowtie'),'.1.ebwtl')):
                                job.add('--large-index',kind='parameter')
                            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                            job.add('--tryhard',kind='parameter')
                            job.add('--best',kind='parameter')
                            job.add('--strata',kind='parameter')
                            job.add('--sam',kind='parameter')
                            job.add('--ff',kind='parameter')
                            job.add('-X',outdir('gene-gene_longest.txt'),kind='parameter',from_file="yes")
                            job.add('',outdir('gene-gene-bowtie_star/'),kind='input',temp_path=temp_flag)
                            job.add('-1',outdir('reads-ids_clip_star_psl_r1.fq'),kind='input',temp_path=temp_flag)
                            job.add('-2',outdir('reads-ids_clip_star_psl_r2.fq'),kind='input',temp_path=temp_flag)
                            job.add('2>',outdir('log_bowtie_reads_mapped-gene-gene-star.stdout.txt'),kind='output',checksum='no',temp_path=temp_flag)
                            job.add('|',kind='parameter')
                            job.add('awk',"""'$3 == "*" { next } { print }'""",kind='parameter')
                            job.add('>',outdir('split_gene-gene_star.sam'),kind='output')
                            job.run()

                            job.add('merge-sam.py',kind='program')
                            job.add('--input',outdir('split_gene-gene_star.sam'),kind='input',temp_path=temp_flag)
                            job.add('--output',outdir('split_gene-gene_star_patch.sam'),kind='output')
                            job.run()

                            job.add('sam2psl.py',kind='program')
                            job.add('--input',outdir('split_gene-gene_star_patch.sam'),kind='input',temp_path=temp_flag)
                            job.add('--output',outdir('split_gene-gene_star_patch.psl'),kind='output')
                            job.run()

                            job.add('analyze_splits_sam.py',kind='program')
                            job.add('--input',outdir('split_gene-gene_star_patch.psl'),kind='input',temp_path=temp_flag)
                            job.add('--output',outdir('split_gene-gene_star_final.psl'),kind='output')
                            job.add('--remove-extra',kind='parameter')
                            job.run()


                        if job.iff(empty(outdir('split_gene-gene_star_final.psl')),id = "#split_gene-gene_star_final#"):
                            job.link(outdir('gene-gene-star.psl'),
                                     outdir('gene-gene-star_more.psl'),
                                     temp_path=temp_flag)
                            job.clean(outdir('split_gene-gene_star.psl'),temp_path=temp_flag)
                        else:
                            job.add('cat',kind='program')
                            job.add('',outdir('split_gene-gene_star_final.psl'),kind='input',temp_path=temp_flag)
                            job.add('',outdir('gene-gene-star.psl'),kind='input',temp_path=temp_flag)
                            job.add('>',outdir('gene-gene-star_more.psl'),kind='output')
                            job.run()



                # find the best unique alignments of reads
                job.add('psl_best_unique_contigs.py',kind='program')
                job.add('--input',outdir('gene-gene-star_more.psl'),kind='input',temp_path=temp_flag)
                job.add('--output',outdir('gene-gene-star_best-unique.psl'),kind='output')
#                if (not empty(outdir('candidate_fusion-genes_further_mark.txt'))) and (not empty(datadir('custom_genes.txt'))):
#                    job.add('--ties',datadir('custom_genes_mark.txt'),kind='input')
                job.add('--ties-overlappings',datadir('ensembl_overlapping_genes.txt'),kind='input')
                job.add('--anchor',length_anchor_star,kind='parameter') # find_fusion_genes_blat.py --threshold_overlap is enough!
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--tmp_dir',tmp_dir,kind='output',checksum='no')
                job.run()

                # more filtering -- remove the reads from the gene-gene junctions
                # which have the pair read mapping on a totally different gene than
                # those involved in the gene-gene junction
                if not options.all_reads_junction:
                    job.add('remove_reads_exon_exon_psl.py',kind='program')
                    job.add('--input_psl',outdir('gene-gene-star_best-unique.psl'),kind='input',temp_path=temp_flag)
                    job.add('--input_transcriptome',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind='input',temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                    job.add('--output_psl',outdir('gene-gene-star_best-unique_gene_pairs.psl'),kind='output')
                    job.run()
                else:
                    job.link(outdir('gene-gene-star_best-unique.psl'),
                             outdir('gene-gene-star_best-unique_gene_pairs.psl'),
                             temp_path=temp_flag)

                job.add('find_fusion_genes_psl.py',kind='program')
                job.add('--input_mappings',outdir('gene-gene-star_best-unique_gene_pairs.psl'),kind='input',temp_path=temp_flag)
                job.add('--input_genegene_fasta',outdir('gene-gene.fa'),kind='input',temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('--input_hugo',datadir('genes_symbols.txt'),kind='input')
                job.add('--input_genes_positions',datadir('genes.txt'),kind='input')
                job.add('--threshold_overlap',length_anchor_star,kind='parameter')
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--output',outdir('candidates_fusion_genes_reads_star.txt'),kind='output')
                job.run()

                # summary the gene-gene mappings
                job.add('build_report_fusions_psl.py',kind='program')
                job.add('--suporting_unique_reads',spanning_reads_star,kind='parameter')
                job.add('--anchor2',length_anchor2,kind='parameter')
                job.add('--input_candidate_fusion_genes_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='input',temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('--input_fastq',outdir('original_important.fq.gz'),kind='input',temp_path=temp_flag if options.skip_bowtie2 and options.skip_bwa else 'no')
                job.add('--input_fusion_psl',outdir('candidates_fusion_genes_reads_star.txt'),kind='input',temp_path=temp_flag)
                job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
                if options.psl_visualization and not empty(datadir('genome.2bit')):
                    job.add('--input_genome_2bit',datadir('genome.2bit'),kind='input')
                    job.add('--psl_alignment_type','web',kind='parameter')
                if options.sam_visualization:
                    job.add('--input_genome_bowtie2',datadir('genome_index2/index'),kind='input')
                    job.add('--sam_alignment','20',kind='parameter')
                    job.add('--threads',options.processes,kind='parameter')
                if options.assembly:
                    job.add('--velvet',kind='parameter')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_STAR.txt'),kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_STAR.zip'),kind='output')
                job.run()

################################################################################
# Bowtie2 alignment
################################################################################
            if not options.skip_bowtie2:

                # find available memory
                job.add('printf',kind='program')
                job.add('"\n============\nMEMORY (before using Bowtie2):\n============\n"',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('free',kind='program')
                job.add('-m',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()

                if nucleotides_gg > options.limit_bowtie2:

                    job.add('split-fasta.py',kind='program')
                    job.add('--size',outdir('gene-gene__nuc.txt'),kind='input')
                    job.add('--seqs',outdir('gene-gene__seq.txt'),kind='input')
                    job.add('--threshold',options.limit_bowtie2,kind='parameter')
                    job.add('-i',outdir('gene-gene.fa'),kind='input')
                    job.add('-o',outdir('gene-gene_split_bowtie2.fa'),kind='output')
                    job.add('-x',outdir('gene-gene_split_bowtie2.len'),kind='output')
                    job.run()

                    parts = [el.strip() for el in file(outdir('gene-gene_split_bowtie2.fa'),'r').readlines()]
                    maxlens = [el.strip() for el in file(outdir('gene-gene_split_bowtie2.len'),'r').readlines()]
                    for i,part in enumerate(parts):

                        gd = "%s_bowtie2/" % (part,)
                        gdi = "%s_bowtie2/index" % (part,)
                        # build the BOWTIE2 index
                        job.add('bowtie2-build',kind='program')
                        job.add('-f',kind='parameter')
                        job.add('--quiet',kind='parameter')
                        job.add('--offrate','1',kind='parameter')
                        job.add('--ftabchars','7',kind='parameter')
                        job.add('',part,kind='input')
                        job.add('',gdi,kind='output',checksum='no')
                        job.add('',gd,kind='output',command_line='no')
                        job.run()


                        # align the unmapped reads using BOWTIE2 on candidate fusion gene-gene
                        job.add('bowtie2',kind='program')
                        job.add('-p',options.processes,kind='parameter',checksum='no')
                        job.add('--phred33',kind='parameter')
                        job.add('--no-unal',kind='parameter')
                        job.add('--local',kind='parameter')
                        job.add('-N','1',kind='parameter') # new
                        job.add('-R','3',kind='parameter') # new
                        job.add('-D','20',kind='parameter') # new
                        job.add('-k','5',kind='parameter')
                        job.add('-L','20',kind='parameter')
                        job.add('-x',gdi,kind='input',checksum='no')
                        job.add('-x',gd,kind='input',command_line='no',temp_path=temp_flag)
                        job.add('-U',outdir('reads_gene-gene_no-str.fq'),kind='input')
                        job.add('-S',outdir('gene-gene-bowtie2.sam.')+str(i),kind='output')
                        job.add('2>',outdir('log_bowtie2_reads-gene-gene.stdout.txt.')+str(i),kind='output',checksum='no')
                        job.run()

                        job.clean(outdir('log_bowtie2_reads-gene-gene.stdout.txt.')+str(i),temp_path=temp_flag)

                        job.add('sam2psl.py',kind='program')
                        job.add('--input',outdir('gene-gene-bowtie2.sam.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--output','-',kind='output')
                        job.add('|',kind='parameter')
                        job.add('LC_ALL=C',kind='parameter')
                        job.add('sort',kind='parameter')
                        job.add('-k','10,10',kind='parameter')
                        job.add('-k','14,14',kind='parameter')
                        job.add('-k','12,12n',kind='parameter')
                        job.add('-k','13,13n',kind='parameter')
                        job.add('-t',"'\t'",kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('>',outdir('gene-gene-bowtie2.psl.')+str(i),kind='output')
                        job.run()

                        job.add('analyze_splits_sam.py',kind='program')
                        job.add('--input',outdir('gene-gene-bowtie2.psl.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('gene-gene-bowtie2_final.psl.')+str(i),kind='output')
                        job.add('--clipped-reads-ids',outdir('reads-ids_clip_psl_bowtie2.txt.')+str(i),kind='output')
                        job.add('--clipped-reads-refs',outdir('reads-refs_clip_psl_bowtie2.txt.')+str(i),kind='output')
                        job.add('--clip-min',length_anchor_bowtie2,kind='parameter')
                        job.run()

                        if job.iff(empty(outdir('reads-ids_clip_psl_bowtie2.txt.')+str(i)),id = "#reads-ids-clip-psl-bowtie2."+str(i)+"#"):
                            job.clean(outdir('reads-ids_clip_psl_bowtie2.txt.')+str(i),temp_path=temp_flag)
                            job.clean(outdir('reads-refs_clip_psl_bowtie2.txt.')+str(i),temp_path=temp_flag)
                            job.clean(part,temp_path=temp_flag)
                        else:
                            job.add('LC_ALL=C',kind='program')
                            job.add('sort',kind='parameter')
                            if sort_buffer:
                                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                            if sort_parallel:
                                job.add('--parallel',options.processes,kind='parameter',checksum='no')
                            if sort_lzop_compress:
                                job.add('--compress-program','lzop',kind='parameter',checksum='no')
                            elif sort_gzip_compress:
                                job.add('--compress-program','gzip',kind='parameter',checksum='no')
                            job.add('-T',tmp_dir,kind='parameter',checksum='no')
                            job.add('',outdir('reads-ids_clip_psl_bowtie2.txt.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('>',outdir('reads-ids_clip_bowtie2_psl_uniq.txt.')+str(i),kind='output')
                            job.run()

                            job.add('cut',kind='program')
                            job.add('-f1',outdir('reads-ids_clip_bowtie2_psl_uniq.txt.')+str(i),kind='input')
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('|',kind='parameter')
                            job.add('LC_ALL=C',kind='parameter')
                            job.add('sort',kind='parameter')
                            if sort_buffer:
                                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                            if sort_parallel:
                                job.add('--parallel',options.processes,kind='parameter',checksum='no')
                            if sort_lzop_compress:
                                job.add('--compress-program','lzop',kind='parameter',checksum='no')
                            elif sort_gzip_compress:
                                job.add('--compress-program','gzip',kind='parameter',checksum='no')
                            job.add('-T',tmp_dir,kind='parameter',checksum='no')
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('|',kind='parameter')
                            job.add('seqtk',kind='parameter')
                            job.add('subseq',kind='parameter')
                            job.add('',outdir('reads_gene-gene_no-str.fq'),kind='input')
                            job.add('-',kind='parameter')
                            job.add('>',outdir('reads-ids_clip_bowtie2_psl.fq.')+str(i),kind='output')
                            job.run()

                            job.add('split-reads.py',kind='program')
                            job.add('--input',outdir('reads-ids_clip_bowtie2_psl.fq.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('--list',outdir('reads-ids_clip_bowtie2_psl_uniq.txt.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('--output-1',outdir('reads-ids_clip_bowtie2_psl_r1.fq.')+str(i),kind='output')
                            job.add('--output-2',outdir('reads-ids_clip_bowtie2_psl_r2.fq.')+str(i),kind='output')
                            job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                            job.run(error_message = ("If this fails due to a memory error then lowering the "+
                                                     "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
                                                     "of FusionCatcher and running it again might help!"))


                            job.add('LC_ALL=C',kind='program')
                            job.add('sort',kind='parameter')
                            if sort_buffer:
                                job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                            if sort_parallel:
                                job.add('--parallel',options.processes,kind='parameter',checksum='no')
                            if sort_lzop_compress:
                                job.add('--compress-program','lzop',kind='parameter',checksum='no')
                            elif sort_gzip_compress:
                                job.add('--compress-program','gzip',kind='parameter',checksum='no')
                            job.add('-T',tmp_dir,kind='parameter',checksum='no')
                            job.add('',outdir('reads-refs_clip_psl_bowtie2.txt.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('|',kind='parameter')
                            job.add('uniq',kind='parameter')
                            job.add('>',outdir('reads-refs_clip_bowtie2_psl_uniq.txt.')+str(i),kind='output')
                            job.run()

                            gda = "%s_bowtie_bowtie2.fa" % (part,)
                            job.add('seqtk',kind='program')
                            job.add('subseq',kind='parameter')
                            job.add('',part,kind='input',temp_path=temp_flag)
                            job.add('',outdir('reads-refs_clip_bowtie2_psl_uniq.txt.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('>',gda,kind='output')
                            job.run()


                            gdb = "%s_bowtie_bowtie2/" % (part,)
                            job.add('bowtie-build',kind='program')
                            job.add('-f',kind='parameter')
                            job.add('--quiet',kind='parameter')
#                            job.add('--ntoa',kind='parameter')
                            job.add('--offrate','1',kind='parameter')
                            job.add('--ftabchars','5',kind='parameter')
                            #job.add('',part,kind='input',temp_path=temp_flag)
                            job.add('',gda,kind='input',temp_path=temp_flag)
                            job.add('',gdb,kind='output',checksum='no')
                            job.add('',gdb,kind='output',command_line='no')
                            job.run()

                            # map using bowtie
                            job.add('bowtie',kind='program')
                            job.add('-t',kind='parameter')
                            #job.add('-q',kind='parameter')
                            job.add('-a',kind='parameter')
                            job.add('-v',options.mismatches,kind='parameter')
                            job.add('-p',options.processes,kind='parameter',checksum='no')
                            if os.path.isfile(os.path.join(gdb,'.1.ebwtl')):
                                job.add('--large-index',kind='parameter')
                            job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                            job.add('--tryhard',kind='parameter')
                            job.add('--best',kind='parameter')
                            job.add('--strata',kind='parameter')
                            job.add('--sam',kind='parameter')
                            job.add('--ff',kind='parameter')
                            #job.add('-X',outdir('gene-gene_longest.txt'),kind='parameter',from_file="yes")
                            job.add('-X',maxlens[i],kind='parameter',from_file="yes")
                            job.add('',gdb,kind='input',temp_path=temp_flag)
                            job.add('-1',outdir('reads-ids_clip_bowtie2_psl_r1.fq.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('-2',outdir('reads-ids_clip_bowtie2_psl_r2.fq.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('2>',outdir('log_bowtie_reads_mapped-gene-gene-bowtie2.stdout.txt.')+str(i),kind='output',checksum='no',temp_path=temp_flag)
                            job.add('|',kind='parameter')
                            job.add('awk',"""'$3 == "*" { next } { print }'""",kind='parameter')
                            job.add('>',outdir('split_gene-gene_bowtie2.sam.')+str(i),kind='output')
                            job.run()

                            job.add('merge-sam.py',kind='program')
                            job.add('--input',outdir('split_gene-gene_bowtie2.sam.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('--output',outdir('split_gene-gene_bowtie2_patch.sam.')+str(i),kind='output')
                            job.run()

                            job.add('sam2psl.py',kind='program')
                            job.add('--input',outdir('split_gene-gene_bowtie2_patch.sam.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('--output',outdir('split_gene-gene_bowtie2_patch.psl.')+str(i),kind='output')
                            job.run()

                            job.add('analyze_splits_sam.py',kind='program')
                            job.add('--input',outdir('split_gene-gene_bowtie2_patch.psl.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('--output',outdir('split_gene-gene_bowtie2_final.psl.')+str(i),kind='output')
                            job.add('--remove-extra',kind='parameter')
                            job.run()

                        if job.iff(empty(outdir('split_gene-gene_bowtie2_final.psl.')+str(i)),id = "#split_gene-gene_bowtie2_final."+str(i)+"#"):
                            job.link(outdir('gene-gene-bowtie2_final.psl.')+str(i),
                                     outdir('gene-gene-bowtie2_final_more.psl.')+str(i),
                                     temp_path=temp_flag,
                                     dest_list='genegenebowtie2')
                            job.clean(outdir('split_gene-gene_bowtie2_final.psl.')+str(i),temp_path=temp_flag)
                        else:
                            job.add('cat',kind='program')
                            job.add('',outdir('split_gene-gene_bowtie2_final.psl.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('',outdir('gene-gene-bowtie2_final.psl.')+str(i),kind='input',temp_path=temp_flag)
                            job.add('>',outdir('gene-gene-bowtie2_final_more.psl.')+str(i),kind='output',dest_list='genegenebowtie2')
                            job.run()


                    #job.clean(outdir('gene-gene_split_bowtie2.fa'),temp_path=temp_flag)
                    job.clean(outdir('reads_gene-gene_no-str.fq'),temp_path=temp_flag if options.skip_bwa else 'no')
                    job.sink(job.genegenebowtie2, outdir('gene-gene-bowtie2_final_more.psl.txt'))

                    job.add('concatenate.py',kind='program')
                    job.add('-f',outdir('gene-gene-bowtie2_final_more.psl.txt'),kind='input',temp_path=temp_flag)
                    job.add('',outdir('gene-gene-bowtie2_final_more.psl'),kind='output')
                    job.run()

#                    for tfile in job.genegenebowtie2:
#                        job.clean(tfile,temp_path=temp_flag)
                    job.clean(job.genegenebowtie2,temp_path=temp_flag)
                    
                else:
                    # build the BOWTIE2 index
                    job.add('bowtie2-build',kind='program')
                    job.add('-f',kind='parameter')
                    job.add('--quiet',kind='parameter')
                    job.add('--offrate','1',kind='parameter')
                    job.add('--ftabchars','7',kind='parameter')
                    job.add('',outdir('gene-gene.fa'),kind='input')
                    job.add('',outdir('gene-gene-bowtie2/index'),kind='output',checksum='no')
                    job.add('',outdir('gene-gene-bowtie2/'),kind='output',command_line='no')
                    job.run()


                    # align the unmapped reads using BOWTIE2 on candidate fusion gene-gene
                    job.add('bowtie2',kind='program')
                    job.add('-p',options.processes,kind='parameter',checksum='no')
                    job.add('--phred33',kind='parameter')
                    job.add('--no-unal',kind='parameter')
                    job.add('--local',kind='parameter')
                    job.add('-N','1',kind='parameter') # new
                    job.add('-R','3',kind='parameter') # new
                    job.add('-D','20',kind='parameter') # new
                    job.add('-k','5',kind='parameter')
                    job.add('-L','20',kind='parameter')
                    job.add('-x',outdir('gene-gene-bowtie2/index'),kind='input',checksum='no')
                    job.add('-x',outdir('gene-gene-bowtie2/'),kind='input',command_line='no',temp_path=temp_flag)
                    job.add('-U',outdir('reads_gene-gene_no-str.fq'),kind='input')
                    job.add('-S',outdir('gene-gene-bowtie2.sam'),kind='output')
                    job.add('2>',outdir('log_bowtie2_reads-gene-gene.stdout.txt'),kind='output',checksum='no')
                    job.run()
                    # -D 20 -R 3 -N 1 -L 20 => almost like bwa-mem
                    # tried  -D 20 -R 3 -N 0 -i 'S,1,0.5' but it was slow

                    job.clean(outdir('log_bowtie2_reads-gene-gene.stdout.txt'),temp_path=temp_flag)

                    job.add('sam2psl.py',kind='program')
                    job.add('--input',outdir('gene-gene-bowtie2.sam'),kind='input',temp_path=temp_flag)
                    job.add('--output','-',kind='output')
                    job.add('|',kind='parameter')
                    job.add('LC_ALL=C',kind='parameter')
                    job.add('sort',kind='parameter')
                    job.add('-k','10,10',kind='parameter')
                    job.add('-k','14,14',kind='parameter')
                    job.add('-k','12,12n',kind='parameter')
                    job.add('-k','13,13n',kind='parameter')
                    job.add('-t',"'\t'",kind='parameter')
                    if sort_buffer:
                        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                    if sort_parallel:
                        job.add('--parallel',options.processes,kind='parameter',checksum='no')
                    if sort_lzop_compress:
                        job.add('--compress-program','lzop',kind='parameter',checksum='no')
                    elif sort_gzip_compress:
                        job.add('--compress-program','gzip',kind='parameter',checksum='no')
                    job.add('-T',tmp_dir,kind='parameter',checksum='no')
                    job.add('>',outdir('gene-gene-bowtie2.psl'),kind='output')
                    job.run()

                    job.add('analyze_splits_sam.py',kind='program')
                    job.add('--input',outdir('gene-gene-bowtie2.psl'),kind='input',temp_path=temp_flag)
                    job.add('--output',outdir('gene-gene-bowtie2_final.psl'),kind='output')
                    job.add('--clipped-reads-ids',outdir('reads-ids_clip_psl_bowtie2.txt'),kind='output')
                    job.add('--clipped-reads-refs',outdir('reads-refs_clip_psl_bowtie2.txt'),kind='output')
                    job.add('--clip-min',length_anchor_bowtie2,kind='parameter')
                    job.run()

                    if job.iff(empty(outdir('reads-ids_clip_psl_bowtie2.txt')),id = "#reads-ids-clip-psl-bowtie2#"):
                        job.clean(outdir('reads-ids_clip_psl_bowtie2.txt'),temp_path=temp_flag)
                        job.clean(outdir('reads-refs_clip_psl_bowtie2.txt'),temp_path=temp_flag)
                    else:
                        job.add('LC_ALL=C',kind='program')
                        job.add('sort',kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('',outdir('reads-ids_clip_psl_bowtie2.txt'),kind='input',temp_path=temp_flag)
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('>',outdir('reads-ids_clip_bowtie2_psl_uniq.txt'),kind='output')
                        job.run()

                        job.add('cut',kind='program')
                        job.add('-f1',outdir('reads-ids_clip_bowtie2_psl_uniq.txt'),kind='input')
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('|',kind='parameter')
                        job.add('LC_ALL=C',kind='parameter')
                        job.add('sort',kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('|',kind='parameter')
                        job.add('seqtk',kind='parameter')
                        job.add('subseq',kind='parameter')
                        job.add('',outdir('reads_gene-gene_no-str.fq'),kind='input',temp_path=temp_flag if options.skip_bwa else 'no')
                        job.add('-',kind='parameter')
                        job.add('>',outdir('reads-ids_clip_bowtie2_psl.fq'),kind='output')
                        job.run()

                        job.add('split-reads.py',kind='program')
                        job.add('--input',outdir('reads-ids_clip_bowtie2_psl.fq'),kind='input',temp_path=temp_flag)
                        job.add('--list',outdir('reads-ids_clip_bowtie2_psl_uniq.txt'),kind='input',temp_path=temp_flag)
                        job.add('--output-1',outdir('reads-ids_clip_bowtie2_psl_r1.fq'),kind='output')
                        job.add('--output-2',outdir('reads-ids_clip_bowtie2_psl_r2.fq'),kind='output')
                        job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                        job.run(error_message = ("If this fails due to a memory error then lowering the "+
                                                 "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
                                                 "of FusionCatcher and running it again might help!"))


                        job.add('LC_ALL=C',kind='program')
                        job.add('sort',kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('',outdir('reads-refs_clip_psl_bowtie2.txt'),kind='input',temp_path=temp_flag)
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('>',outdir('reads-refs_clip_bowtie2_psl_uniq.txt'),kind='output')
                        job.run()

                        job.add('seqtk',kind='program')
                        job.add('subseq',kind='parameter')
                        job.add('',outdir('gene-gene.fa'),kind='input')
                        job.add('',outdir('reads-refs_clip_bowtie2_psl_uniq.txt'),kind='input',temp_path=temp_flag)
                        job.add('>',outdir('gene-gene-bowtie_bowtie2.fa'),kind='output')
                        job.run()

                        job.add('bowtie-build',kind='program')
                        job.add('-f',kind='parameter')
                        job.add('--quiet',kind='parameter')
#                        job.add('--ntoa',kind='parameter')
                        job.add('--offrate','1',kind='parameter')
                        job.add('--ftabchars','5',kind='parameter')
                        job.add('',outdir('gene-gene-bowtie_bowtie2.fa'),kind='input',temp_path=temp_flag)
                        job.add('',outdir('gene-gene-bowtie_bowtie2/'),kind='output',checksum='no')
                        job.add('',outdir('gene-gene-bowtie_bowtie2/'),kind='output',command_line='no')
                        job.run()

                        # map using bowtie
                        job.add('bowtie',kind='program')
                        job.add('-t',kind='parameter')
                        #job.add('-q',kind='parameter')
                        job.add('-a',kind='parameter')
                        job.add('-v',options.mismatches,kind='parameter')
                        job.add('-p',options.processes,kind='parameter',checksum='no')
                        if os.path.isfile(os.path.join(outdir('gene-gene-bowtie'),'.1.ebwtl')):
                            job.add('--large-index',kind='parameter')
                        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                        job.add('--tryhard',kind='parameter')
                        job.add('--best',kind='parameter')
                        job.add('--strata',kind='parameter')
                        job.add('--sam',kind='parameter')
                        job.add('--ff',kind='parameter')
                        job.add('-X',outdir('gene-gene_longest.txt'),kind='parameter',from_file="yes")
                        job.add('',outdir('gene-gene-bowtie_bowtie2/'),kind='input',temp_path=temp_flag)
                        job.add('-1',outdir('reads-ids_clip_bowtie2_psl_r1.fq'),kind='input',temp_path=temp_flag)
                        job.add('-2',outdir('reads-ids_clip_bowtie2_psl_r2.fq'),kind='input',temp_path=temp_flag)
                        job.add('2>',outdir('log_bowtie_reads_mapped-gene-gene-bowtie2.stdout.txt'),kind='output',checksum='no',temp_path=temp_flag)
                        job.add('|',kind='parameter')
                        job.add('awk',"""'$3 == "*" { next } { print }'""",kind='parameter')
                        job.add('>',outdir('split_gene-gene_bowtie2.sam'),kind='output')
                        job.run()

                        job.add('merge-sam.py',kind='program')
                        job.add('--input',outdir('split_gene-gene_bowtie2.sam'),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('split_gene-gene_bowtie2_patch.sam'),kind='output')
                        job.run()

                        job.add('sam2psl.py',kind='program')
                        job.add('--input',outdir('split_gene-gene_bowtie2_patch.sam'),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('split_gene-gene_bowtie2_patch.psl'),kind='output')
                        job.run()

                        job.add('analyze_splits_sam.py',kind='program')
                        job.add('--input',outdir('split_gene-gene_bowtie2_patch.psl'),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('split_gene-gene_bowtie2_final.psl'),kind='output')
                        job.add('--remove-extra',kind='parameter')
                        job.run()


                    if job.iff(empty(outdir('split_gene-gene_bowtie2_final.psl')),id = "#split_gene-gene_bowtie2_final#"):
                        job.link(outdir('gene-gene-bowtie2_final.psl'),
                                 outdir('gene-gene-bowtie2_final_more.psl'),
                                 temp_path=temp_flag)
                        job.clean(outdir('split_gene-gene_bowtie2_final.psl'),temp_path=temp_flag)
                    else:
                        job.add('cat',kind='program')
                        job.add('',outdir('split_gene-gene_bowtie2_final.psl'),kind='input',temp_path=temp_flag)
                        job.add('',outdir('gene-gene-bowtie2_final.psl'),kind='input',temp_path=temp_flag)
                        job.add('>',outdir('gene-gene-bowtie2_final_more.psl'),kind='output')
                        job.run()

                # find the best unique alignments of reads
                job.add('psl_best_unique_contigs.py',kind='program')
                job.add('--input',outdir('gene-gene-bowtie2_final_more.psl'),kind='input',temp_path=temp_flag)
                job.add('--output',outdir('gene-gene-bowtie2_best-unique.psl'),kind='output')
#                if (not empty(outdir('candidate_fusion-genes_further_mark.txt'))) and (not empty(datadir('custom_genes.txt'))):
#                    job.add('--ties',datadir('custom_genes_mark.txt'),kind='output')
                job.add('--ties-overlappings',datadir('ensembl_overlapping_genes.txt'),kind='input')
                job.add('--anchor',length_anchor_bowtie2,kind='parameter') # find_fusion_genes_blat.py --threshold_overlap is enough!
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--tmp_dir',tmp_dir,kind='output',checksum='no')
                job.run()

                # more filtering -- remove the reads from the gene-gene junctions
                # which have the pair read mapping on a totally different gene than
                # those involved in the gene-gene junction
                if not options.all_reads_junction:
                    job.add('remove_reads_exon_exon_psl.py',kind='program')
                    job.add('--input_psl',outdir('gene-gene-bowtie2_best-unique.psl'),kind='input',temp_path=temp_flag)
                    job.add('--input_transcriptome',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind='input',temp_path=temp_flag if options.skip_bwa else 'no')
                    job.add('--output_psl',outdir('gene-gene-bowtie2_best-unique_gene_pairs.psl'),kind='output')
                    job.run()
                else:
                    job.link(outdir('gene-gene-bowtie2_best-unique.psl'),
                             outdir('gene-gene-bowtie2_best-unique_gene_pairs.psl'),
                             temp_path=temp_flag)

                job.add('find_fusion_genes_psl.py',kind='program')
                job.add('--input_mappings',outdir('gene-gene-bowtie2_best-unique_gene_pairs.psl'),kind='input',temp_path=temp_flag)
                job.add('--input_genegene_fasta',outdir('gene-gene.fa'),kind='input',temp_path=temp_flag if options.skip_bwa else 'no')
                job.add('--input_hugo',datadir('genes_symbols.txt'),kind='input')
                job.add('--input_genes_positions',datadir('genes.txt'),kind='input')
                job.add('--threshold_overlap',length_anchor_bowtie2,kind='parameter')
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--output',outdir('candidates_fusion_genes_reads_bowtie2.txt'),kind='output')
                job.run()

                # summary the gene-gene mappings
                job.add('build_report_fusions_psl.py',kind='program')
                job.add('--suporting_unique_reads',spanning_reads_bowtie2,kind='parameter')
                job.add('--anchor2',length_anchor2,kind='parameter')
                job.add('--input_candidate_fusion_genes_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='input',temp_path=temp_flag if options.skip_bwa else 'no')
                job.add('--input_fastq',outdir('original_important.fq.gz'),kind='input',temp_path=temp_flag if options.skip_bwa else 'no')
                job.add('--input_fusion_psl',outdir('candidates_fusion_genes_reads_bowtie2.txt'),kind='input',temp_path=temp_flag)
                job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
                if options.psl_visualization and not empty(datadir('genome.2bit')):
                    job.add('--input_genome_2bit',datadir('genome.2bit'),kind='input')
                    job.add('--psl_alignment_type','web',kind='parameter')
                if options.sam_visualization:
                    job.add('--input_genome_bowtie2',datadir('genome_index2/index'),kind='input')
                    job.add('--sam_alignment','20',kind='parameter')
                    job.add('--threads',options.processes,kind='parameter')
                if options.assembly:
                    job.add('--velvet',kind='parameter')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BOWTIE2.txt'),kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BOWTIE2.zip'),kind='output')
                job.run()


################################################################################
# BWA alignment
################################################################################
            if not options.skip_bwa:

                if options.skip_star:
                    # BWA is removing the /1 and /2 from the end of the reads names
                    # changing "/1" and "/2" into "-1" "-2" such that BWA does not remove them
                    job.add('sed',kind='program')
                    job.add("""'s/\/\([1-2]$\)/\-\\1/;n;n;n'""",outdir('reads_gene-gene_no-str.fq'),kind='input',temp_path=temp_flag)
                    job.add('>',outdir('reads_gene-gene_no-str_fixed.fq'),kind='output')
                    job.run()
                else:
                    job.clean(outdir('reads_gene-gene_no-str.fq'),temp_path=temp_flag)

                job.add('split-fasta.py',kind='program')
                job.add('--size',outdir('gene-gene__nuc.txt'),kind='input')
                job.add('--seqs',outdir('gene-gene__seq.txt'),kind='input')
                job.add('--threshold',options.limit_bwa,kind='parameter')
                job.add('-i',outdir('gene-gene.fa'),kind='input')
                job.add('-o',outdir('gene-gene_split_bwa.fa'),kind='output')
                job.add('-x',outdir('gene-gene_split_bwa.len'),kind='output')
                job.run()

                # find available memory
                job.add('printf',kind='program')
                job.add('"\n============\nMEMORY (before using BWA):\n============\n"',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()
                job.add('free',kind='program')
                job.add('-m',kind='parameter')
                job.add('>>',info_file,kind='output')
                job.run()

                parts = [el.strip() for el in file(outdir('gene-gene_split_bwa.fa'),'r').readlines()]
                maxlens = [el.strip() for el in file(outdir('gene-gene_split_bwa.len'),'r').readlines()]
                for i,part in enumerate(parts):


                    job.add('bwa',kind='program')
                    job.add('index',kind='parameter')
                    job.add('-a',kind='parameter')
                    job.add('bwtsw',kind='parameter')
                    job.add('',part,kind='input')
                    job.add('',part+'.amb',kind='output',command_line='no')
                    job.run()

                    job.add('bwa',kind='program')
                    job.add('mem',kind='parameter')
                    job.add('-k',length_anchor_bwa,kind='parameter')
                    job.add('-a',kind='parameter')
                    job.add('-T','1',kind='parameter')
                    job.add('-t',options.processes,kind='parameter',checksum='no')
                    job.add('',part,kind='input')
                    job.add('',part+'.amb',kind='input',command_line='no',temp_path=temp_flag)
                    job.add('',part+'.ann',kind='input',command_line='no',temp_path=temp_flag)
                    job.add('',part+'.bwt',kind='input',command_line='no',temp_path=temp_flag)
                    job.add('',part+'.pac',kind='input',command_line='no',temp_path=temp_flag)
                    job.add('',part+'.sa',kind='input',command_line='no',temp_path=temp_flag)
                    job.add('',outdir('reads_gene-gene_no-str_fixed.fq'),kind='input')
#                    job.add('>',outdir('gene-gene-bwa.sam.')+str(i),kind='output')
#                    job.add('2>',outdir('log_bwa_reads-gene-gene.stdout.txt.')+str(i),kind='output',checksum='no')
#                    job.run()
#                    job.clean(outdir('log_bwa_reads-gene-gene.stdout.txt.')+str(i),temp_path=temp_flag)

                    job.add('|',kind='parameter')
                    job.add('sam2psl.py',kind='parameter')
                    job.add('--input','-',kind='parameter')
                    job.add('--output','-',kind='output')
                    job.add('|',kind='parameter')
                    job.add('sed',kind='parameter')
#                    job.add("""'s/\-\([1-2]\\t\)/\/\\1/'""",outdir('gene-gene-bwa.sam.')+str(i),kind='input',temp_path=temp_flag)
                    job.add("""'s/\-\([1-2]\\t\)/\/\\1/'""",kind='parameter')
                    job.add('|',kind='parameter')
                    job.add('LC_ALL=C',kind='parameter')
                    job.add('sort',kind='parameter')
                    job.add('-k','10,10',kind='parameter')
                    job.add('-k','14,14',kind='parameter')
                    job.add('-k','12,12n',kind='parameter')
                    job.add('-k','13,13n',kind='parameter')
                    job.add('-t',"'\t'",kind='parameter')
                    if sort_buffer:
                        job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                    if sort_parallel:
                        job.add('--parallel',options.processes,kind='parameter',checksum='no')
                    if sort_lzop_compress:
                        job.add('--compress-program','lzop',kind='parameter',checksum='no')
                    elif sort_gzip_compress:
                        job.add('--compress-program','gzip',kind='parameter',checksum='no')
                    job.add('-T',tmp_dir,kind='parameter',checksum='no')
                    job.add('>',outdir('gene-gene-bwa.psl.')+str(i),kind='output')
                    job.run()

                    job.add('analyze_splits_sam.py',kind='program')
                    job.add('--input',outdir('gene-gene-bwa.psl.')+str(i),kind='input',temp_path=temp_flag)
                    job.add('--output',outdir('gene-gene-bwa_final.psl.')+str(i),kind='output')
                    job.add('--clipped-reads-ids',outdir('reads-ids_clip_psl_bwa.txt.')+str(i),kind='output')
                    job.add('--clip-min',length_anchor_bwa,kind='parameter')
                    job.run()


                    if job.iff(empty(outdir('reads-ids_clip_psl_bwa.txt.')+str(i)),id = "#reads-ids-clip-psl-bwa."+str(i)+"#"):
                        job.clean(outdir('reads-ids_clip_psl_bwa.txt.')+str(i),temp_path=temp_flag)
                        job.clean(part,temp_path=temp_flag)
                    else:
                        job.add('LC_ALL=C',kind='program')
                        job.add('sort',kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('',outdir('reads-ids_clip_psl_bwa.txt.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('>',outdir('reads-ids_clip_bwa_psl_uniq.txt.')+str(i),kind='output')
                        job.run()

                        job.add('cut',kind='program')
                        job.add('-f1',outdir('reads-ids_clip_bwa_psl_uniq.txt.')+str(i),kind='input')
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('|',kind='parameter')
                        job.add('LC_ALL=C',kind='parameter')
                        job.add('sort',kind='parameter')
                        if sort_buffer:
                            job.add('--buffer-size',sort_buffer,kind='parameter',checksum='no')
                        if sort_parallel:
                            job.add('--parallel',options.processes,kind='parameter',checksum='no')
                        if sort_lzop_compress:
                            job.add('--compress-program','lzop',kind='parameter',checksum='no')
                        elif sort_gzip_compress:
                            job.add('--compress-program','gzip',kind='parameter',checksum='no')
                        job.add('-T',tmp_dir,kind='parameter',checksum='no')
                        job.add('|',kind='parameter')
                        job.add('uniq',kind='parameter')
                        job.add('|',kind='parameter')
                        job.add('seqtk',kind='parameter')
                        job.add('subseq',kind='parameter')
                        job.add('',outdir('reads_gene-gene_no-str_fixed.fq'),kind='input')
                        job.add('-',kind='parameter')
                        job.add('>',outdir('reads-ids_clip_bwa_psl.fq.')+str(i),kind='output')
                        job.run()

                        job.add('split-reads.py',kind='program')
                        job.add('--input',outdir('reads-ids_clip_bwa_psl.fq.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--list',outdir('reads-ids_clip_bwa_psl_uniq.txt.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--output-1',outdir('reads-ids_clip_bwa_psl_r1.fq.')+str(i),kind='output')
                        job.add('--output-2',outdir('reads-ids_clip_bwa_psl_r2.fq.')+str(i),kind='output')
                        job.add('--buffer-size',options.extract_buffer_size,kind='parameter',checksum='no')
                        job.run(error_message = ("If this fails due to a memory error then lowering the "+
                                                 "buffer size (to 50% or 25%) using the command line option --extra-buffer-size "+
                                                 "of FusionCatcher and running it again might help!"))

                        gdb = "%s_bowtie_bwa/" % (part,)
                        job.add('bowtie-build',kind='program')
                        job.add('-f',kind='parameter')
                        job.add('--quiet',kind='parameter')
#                        job.add('--ntoa',kind='parameter')
                        job.add('--offrate','1',kind='parameter')
                        job.add('--ftabchars','5',kind='parameter')
                        job.add('',part,kind='input',temp_path=temp_flag)
                        job.add('',gdb,kind='output',checksum='no')
                        job.add('',gdb,kind='output',command_line='no')
                        job.run()

                        # map using bowtie
                        job.add('bowtie',kind='program')
                        job.add('-t',kind='parameter')
                        #job.add('-q',kind='parameter')
                        job.add('-a',kind='parameter')
                        job.add('-v',options.mismatches,kind='parameter')
                        job.add('-p',options.processes,kind='parameter',checksum='no')
                        if os.path.isfile(os.path.join(gdb,'.1.ebwtl')):
                            job.add('--large-index',kind='parameter')
                        job.add('--chunkmbs',options.chunkmbs,kind='parameter',checksum='no')
                        job.add('--tryhard',kind='parameter')
                        job.add('--best',kind='parameter')
                        job.add('--strata',kind='parameter')
                        job.add('--sam',kind='parameter')
                        job.add('--ff',kind='parameter')
                        #job.add('-X',outdir('gene-gene_longest.txt'),kind='parameter',from_file="yes")
                        job.add('-X',maxlens[i],kind='parameter',from_file="yes")
                        job.add('',gdb,kind='input',temp_path=temp_flag)
                        job.add('-1',outdir('reads-ids_clip_bwa_psl_r1.fq.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('-2',outdir('reads-ids_clip_bwa_psl_r2.fq.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('2>',outdir('log_bowtie_reads_mapped-gene-gene-bwa.stdout.txt.')+str(i),kind='output',checksum='no',temp_path=temp_flag)
                        job.add('|',kind='parameter')
                        job.add('awk',"""'$3 == "*" { next } { print }'""",kind='parameter')
                        job.add('>',outdir('split_gene-gene_bwa.sam.')+str(i),kind='output')
                        job.run()

                        job.add('merge-sam.py',kind='program')
                        job.add('--input',outdir('split_gene-gene_bwa.sam.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('split_gene-gene_bwa_patch.sam.')+str(i),kind='output')
                        job.run()

                        job.add('sam2psl.py',kind='program')
                        job.add('--input',outdir('split_gene-gene_bwa_patch.sam.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('split_gene-gene_bwa_patch.psl.')+str(i),kind='output')
                        job.run()

                        job.add('analyze_splits_sam.py',kind='program')
                        job.add('--input',outdir('split_gene-gene_bwa_patch.psl.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('--output',outdir('split_gene-gene_bwa_final.psl.')+str(i),kind='output')
                        job.add('--remove-extra',kind='parameter')
                        job.run()

                    if job.iff(empty(outdir('split_gene-gene_bwa_final.psl.')+str(i)),id = "#split_gene-gene_bwa_final."+str(i)+"#"):
                        job.link(outdir('gene-gene-bwa_final.psl.')+str(i),
                                 outdir('gene-gene-bwa_final_more.psl.')+str(i),
                                 temp_path=temp_flag,
                                 dest_list='genegenebwa')
                        job.clean(outdir('split_gene-gene_bwa_final.psl.')+str(i),temp_path=temp_flag)
                    else:
                        job.add('cat',kind='program')
                        job.add('',outdir('split_gene-gene_bwa_final.psl.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('',outdir('gene-gene-bwa_final.psl.')+str(i),kind='input',temp_path=temp_flag)
                        job.add('>',outdir('gene-gene-bwa_final_more.psl.')+str(i),kind='output',dest_list='genegenebwa')
                        job.run()


                #job.clean(outdir('gene-gene_split_bwa.fa'),temp_path=temp_flag)
                job.clean(outdir('reads_gene-gene_no-str_fixed.fq'),temp_path=temp_flag)
                job.sink(job.genegenebwa, outdir('gene-gene-bwa_final_more.psl.txt'))

                job.add('concatenate.py',kind='program')
                job.add('-f',outdir('gene-gene-bwa_final_more.psl.txt'),kind='input',temp_path=temp_flag)
                job.add('',outdir('gene-gene-bwa_final_more.psl'),kind='output')
                job.run()

#                for tfile in job.genegenebwa:
#                    job.clean(tfile,temp_path=temp_flag)
                job.clean(job.genegenebwa,temp_path=temp_flag)

                # find the best unique alignments of reads
                job.add('psl_best_unique_contigs.py',kind='program')
                job.add('--input',outdir('gene-gene-bwa_final_more.psl'),kind='input',temp_path=temp_flag)
                job.add('--output',outdir('gene-gene-bwa_best-unique.psl'),kind='output')
#                if (not empty(outdir('candidate_fusion-genes_further_mark.txt'))) and (not empty(datadir('custom_genes.txt'))):
#                    job.add('--ties',datadir('custom_genes_mark.txt'),kind='input')
                job.add('--ties-overlappings',datadir('ensembl_overlapping_genes.txt'),kind='input')
                job.add('--anchor',length_anchor_bwa,kind='parameter') # find_fusion_genes_blat.py --threshold_overlap is enough!
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--processes',options.processes,kind='parameter',checksum='no')
                job.add('--tmp_dir',tmp_dir,kind='output',checksum='no')
                job.run()


                # more filtering -- remove the reads from the gene-gene junctions
                # which have the pair read mapping on a totally different gene than
                # those involved in the gene-gene junction
                if not options.all_reads_junction:
                    job.add('remove_reads_exon_exon_psl.py',kind='program')
                    job.add('--input_psl',outdir('gene-gene-bwa_best-unique.psl'),kind='input',temp_path=temp_flag)
                    job.add('--input_transcriptome',outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),kind='input',temp_path=temp_flag)
                    job.add('--output_psl',outdir('gene-gene-bwa_best-unique_gene_pairs.psl'),kind='output')
                    job.run()
                else:
                    job.link(outdir('gene-gene-bwa_best-unique.psl'),
                             outdir('gene-gene-bwa_best-unique_gene_pairs.psl'),
                             temp_path=temp_flag)

                job.add('find_fusion_genes_psl.py',kind='program')
                job.add('--input_mappings',outdir('gene-gene-bwa_best-unique_gene_pairs.psl'),kind='input',temp_path=temp_flag)
                job.add('--input_genegene_fasta',outdir('gene-gene.fa'),kind='input',temp_path=temp_flag)
                job.add('--input_hugo',datadir('genes_symbols.txt'),kind='input')
                job.add('--input_genes_positions',datadir('genes.txt'),kind='input')
                job.add('--threshold_overlap',length_anchor_bwa,kind='parameter')
                job.add('--mismatches',options.mismatches_psl,kind='parameter')
                job.add('--output',outdir('candidates_fusion_genes_reads_bwa.txt'),kind='output')
                job.run()

                # summary the gene-gene mappings
                job.add('build_report_fusions_psl.py',kind='program')
                job.add('--suporting_unique_reads',spanning_reads_bwa,kind='parameter')
                job.add('--anchor2',length_anchor2,kind='parameter')
                job.add('--input_candidate_fusion_genes_reads',outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),kind='input',temp_path=temp_flag)
                job.add('--input_fastq',outdir('original_important.fq.gz'),kind='input',temp_path=temp_flag)
                job.add('--input_fusion_psl',outdir('candidates_fusion_genes_reads_bwa.txt'),kind='input',temp_path=temp_flag)
                job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
                if options.psl_visualization and not empty(datadir('genome.2bit')):
                    job.add('--input_genome_2bit',datadir('genome.2bit'),kind='input')
                    job.add('--psl_alignment_type','web',kind='parameter')
                if options.sam_visualization:
                    job.add('--input_genome_bowtie2',datadir('genome_index2/index'),kind='input')
                    job.add('--sam_alignment','20',kind='parameter')
                    job.add('--threads',options.processes,kind='parameter')
                if options.assembly:
                    job.add('--velvet',kind='parameter')
                job.add('--output_super_summary',outdir('candidate_fusion_genes_summary_BWA.txt'),kind='output')
                job.add('--output_zip_fasta',outdir('supporting-reads_gene-fusions_BWA.zip'),kind='output')
                job.run()

    #
    # merge all reports
    #
    job.add('merge_reports.py',kind='program')
    job.add('--input_ambiguous',outdir('all_ambiguous_genes.txt'),kind='input',temp_path=temp_flag)
    job.add('--input_bowtie',outdir('candidate_fusion_genes_summary_BOWTIE.txt'),kind='input',temp_path=temp_flag)
    job.add('--input_candidate_fusion_genes',outdir('candidate_fusion-genes_further.txt'),kind='input',temp_path=temp_flag)
    job.add('--anchor2',length_anchor2,kind='parameter')
    if (not options.skip_blat) and (not empty(outdir('candidate_fusion_genes_summary_BLAT.txt'))):
        job.add('--supporting_pairs_blat',spanning_pairs_blat,kind='parameter')
        job.add('--supporting_reads_blat',spanning_reads_blat,kind='parameter')
        job.add('--input_blat',outdir('candidate_fusion_genes_summary_BLAT.txt'),kind='input',temp_path=temp_flag)
    if (not options.skip_star) and (not empty(outdir('candidate_fusion_genes_summary_STAR.txt'))):
        job.add('--supporting_pairs_star',spanning_pairs_star,kind='parameter')
        job.add('--supporting_reads_star',spanning_reads_star,kind='parameter')
        job.add('--input_star',outdir('candidate_fusion_genes_summary_STAR.txt'),kind='input',temp_path=temp_flag)
    if (not options.skip_bowtie2) and (not empty(outdir('candidate_fusion_genes_summary_BOWTIE2.txt'))):
        job.add('--supporting_pairs_bowtie2',spanning_pairs_bowtie2,kind='parameter')
        job.add('--supporting_reads_bowtie2',spanning_reads_bowtie2,kind='parameter')
        job.add('--input_bowtie2',outdir('candidate_fusion_genes_summary_BOWTIE2.txt'),kind='input',temp_path=temp_flag)
    if (not options.skip_bwa) and (not empty(outdir('candidate_fusion_genes_summary_BWA.txt'))):
        job.add('--supporting_pairs_bwa',spanning_pairs_bwa,kind='parameter')
        job.add('--supporting_reads_bwa',spanning_reads_bwa,kind='parameter')
        job.add('--input_bwa',outdir('candidate_fusion_genes_summary_BWA.txt'),kind='input',temp_path=temp_flag)
    if not options.long_report:
        job.add('--squish-report',kind='parameter')
    job.add('--output',outdir('final-list_candidate-fusion-genes-temp.txt'),kind='output')
    job.run()

    # predict effect of fusion
    job.add('label_found_fusions.py',kind='program')
    job.add('--data',datadir('readthroughs.txt'),kind='input')
    job.add('--input',outdir('final-list_candidate-fusion-genes-temp.txt'),kind='input',temp_path=temp_flag)
    job.add('--output',outdir('final-list_candidate-fusion-genes-t2.txt'),kind='output')
    job.add('--data-not-commutative',kind='parameter')
    job.add('--label','readthrough',kind='parameter')
    job.run()

    # predict effect of fusion
    job.add('predict_frame.py',kind='program')
    job.add('--gtf',datadir('organism.gtf'),kind='input')
    job.add('--transcripts',datadir('transcripts.fa'),kind='input')
    job.add('--input',outdir('final-list_candidate-fusion-genes-t2.txt'),kind='input',temp_path=temp_flag)
    job.add('--output',outdir('final-list_candidate-fusion-genes.txt'),kind='output')
    job.run()

    # predict effect of fusion
    job.add('build_summary.py',kind='program')
    job.add('--input',outdir('final-list_candidate-fusion-genes.txt'),kind='input')
    job.add('--output',outdir('summary_candidate_fusions.txt'),kind='output')
    job.run()

    if not options.skip_conversion_grch37:
        human = [line for line in file(datadir('version.txt'),'r').readlines() if line.lower().find('homo sapiens')!=-1 or line.lower().find('grch38')!=-1]
        #print "Human Genome GRCh38 found!"
        if human and len(human) == 2:
            # predict effect of fusion
            job.add('liftover.py',kind='program')
            job.add('--input',outdir('final-list_candidate-fusion-genes.txt'),kind='input')
            job.add('--chain',datadir('hg38ToHg19.over.chain.gz'),kind='input')
            job.add('--output',outdir('final-list_candidate-fusion-genes.GRCh37.txt'),kind='output')
            job.add('--tmp_dir',tmp_dir,kind='parameter',checksum='no')
            job.run()

    #
    # CLEANING
    #
    to_delete = [
        outdir('gene-gene_split_blat.len'),
        outdir('gene-gene_split_bowtie2.len'),
        outdir('gene-gene_split_bwa.len'),
        outdir('gene-gene_split_star.len'),
        outdir('pre-fusion')]
    job.clean(to_delete,list_file='yes',temp_path=temp_flag)
    
    to_delete = [
        outdir('exon-exon_junction_cut_split.fa'),
        outdir('gene-gene_split_blat.fa'),
        outdir('gene-gene_split_bowtie2.fa'),
        outdir('gene-gene_split_bwa.fa'),
        outdir('gene-gene_split_star.fa'),
        outdir('all_ambiguous_genes.txt'),
        outdir('candidate_fusion-genes_exon-exon.txt'),
        outdir('candidate_fusion-genes_further.txt'),
        outdir('candidate_fusion-genes_missing_mates.txt'),
        outdir('candidate_fusion-genes_no-offending-reads_supporting_paired-reads.txt'),
        outdir('count_reads_left_after_filtering.txt'),
        outdir('gene-gene.2bit'),
        outdir('gene-gene.fa'),
        outdir('gene-gene-star/'),
        outdir('gene-gene-bowtie2/'),
        outdir('gene-gene-star-results/'),
        outdir('gene-gene_longest.txt'),
        outdir('gene-gene_unique.fa'),
        outdir('original.fq'),
        outdir('original.fq.gz'),
        outdir('originala.fq'),
        outdir('originala.fq.gz'),
        outdir('original_important.fq.gz'),
        outdir('original_important.txt'),
        outdir('reads_filtered_mapped-transcriptome.fq'),
        outdir('reads-filtered_multiple-mappings-genome.fq'),
        outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-p.fq'),
        outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_psl-pp.fq'),
        outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_ex-ex_final.fq'),
        outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final.fq'),
        outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome.fq'),
        outdir('reads_filtered_transcriptome_sorted-read_end.map'),
        outdir('reads_filtered_transcriptome_sorted-read_end_important.map'),
        outdir('reads_gene-gene.fq'),
        outdir('exon-exon_junction_cut__seq.txt'),
        outdir('exon-exon_junction_cut__nuc.txt'),
        outdir('gene-gene_unique__nuc.txt'),
        outdir('gene-gene__nuc.txt'),
        outdir('gene-gene__seq.txt'),
        outdir('gene-gene_unique__seq.txt'),
        outdir('log_lengths_original_reads.txt'),
        outdir('log_lengths_original_reads_plus.txt'),
        outdir('log_lengths_original_reads_final.txt'),
        outdir('log_bowtie_reads_filtered_all-possible-mappings-transcriptome_map.stdout.txt'),
        outdir('log_bowtie_reads_filtered_all-possible-mappings-transcriptome_multiple_map.stdout.txt'),
        outdir('log_bowtie_reads-filtered-out.stdout.txt'),
        outdir('log_bowtie_reads_mapped-exon-exon-fusion-genes_map.stdout.txt'),
        outdir('log_bowtie_reads_mapped-genome.stdout.txt'),
        outdir('log_bowtie_reads_not-mapped-genome_but_mapped-transcriptome.stdout.txt'),
        outdir('log_bowtie_reads_unique-mapped-genome_mapped-transcriptome.stdout.txt'),
        outdir('log_bowtie_reads-unmapped-filtered-out-genome-blat_last.stdout.txt'),
        outdir('log_bowtie_reads-unmapped-filtered-out-genome_last.stdout.txt'),
        outdir('log_bowtie_reads-unmapped-filtered-out-genome.stdout.txt'),
        outdir('log_bowtie_reads-unmapped-filtered-out-transcriptome.stdout.txt'),
        outdir('log_lengths_reads.txt'),
        outdir('log_minimum_length_short_read.txt'),
        outdir('log_number_of_reads_processed.txt'),
        outdir('log_overlaps_error.txt'),
        outdir('log_lengths_reads_gene-gene_no-str.txt'),
        outdir('log_counts_reads_gene-gene_no-str.txt'),
        outdir('candidate_fusion-genes_further_paired-reads.txt'),
        outdir('reads_filtered_not-mapped-genome_not-mapped-transcriptome_final2.txt'),
        outdir('log_removed_single_reads1.txt'),
        outdir('candidate_fusion-genes_further_mark.txt'),
        outdir('star_sdjboverhang.txt'),
        outdir('split_gene-gene_star_final.psl'),
        outdir('split_gene-gene_bowtie2_final.psl'),
        outdir('split_gene-gene_bwa_final.psl'),
        outdir('single/'),
        outdir('restart.sh')]
    job.clean(to_delete,temp_path=temp_flag)
    
    if urls:
        job.unprotect(urls)
        job.clean(new_input_output,temp_path = temp_flag)


    job.clean(tmp_dir,temp_path='yes')


    job.close()

    #
    #The End!
    #
