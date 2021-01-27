#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It produces a report with the summary of the fusion genes found. Also
FASTQ and FASTA files containing the supporting reads corresponding to each
fusion gene is generated.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2021 Daniel Nicorici

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

This file is running/executing/using BLAT.

"""
import sys
import os
import optparse
import gc
import string
import zipfile
import Bio.SeqIO
import datetime
import tempfile
import shutil
import gzip

empty_zip_data = 'PK\x05\x06\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'

ttable = string.maketrans("ACGTN","TGCAA") # global
#ttable = string.maketrans("ACGTYRSWKMBDHV-","TGCARYSWMKVHDB-")

mapping_solexa2sanger = "".join([chr(0) for ascii in range(0, 59)]
                        + [chr(33 + int(round(Bio.SeqIO.QualityIO.phred_quality_from_solexa(q)))) for q in range(-5, 62 + 1)]
                        + [chr(0) for ascii in range(127, 256)])

mapping_illumina2sanger = "".join([chr(0) for ascii in range(0, 64)]
                          + [chr(33 + q) for q in range(0, 62 + 1)]
                          + [chr(0) for ascii in range(127, 256)])

def solexa2sanger(qual):
    return qual.translate(mapping_solexa2sanger)

def illumina2sanger(qual):
    return qual.translate(mapping_illumina2sanger)

def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

def myorder(a,b):
    return (a,b) if a <= b else (b,a)

def dnaReverseComplement(seq):
    seq = seq.upper()
    seq = seq.translate(ttable)
    return seq[::-1]


def reads_from_fastq_file(file_name, size_read_buffer = 10**8):
    fid = None
    if file_name == '-':
        fid = sys.stdin
    elif file_name.lower().endswith('.gz'):
        fid = gzip.open(file_name,'r')
    else:
        fid = open(file_name,'r')
    piece = [None,None,None,None]
    ij = 0
    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for line in lines:
            ij = ij + 1
            piece[ij-1] = line
            if ij == 4:
                bucket = (piece[0].rstrip('\r\n')[1:],
                          piece[1].rstrip('\r\n'),
                          piece[3].rstrip('\r\n'))
                yield bucket
                piece = [None,None,None,None]
                ij = 0
    fid.close()

def delete_file(some_file):
    if os.path.isfile(some_file) or os.path.islink(some_file):
        os.remove(some_file)
    elif os.path.isdir(some_file):
        shutil.rmtree(some_file)

def give_me_psl(fasta, twobit, blat_dir = None, tmp_dir = None, align_type = 'web'):
    # give as input a file as a list of strings it runs BLAT and it returns
    # the PSL output as a list of strings
    fasta_file = give_me_temp_filename(tmp_dir = tmp_dir)
    psl_file = give_me_temp_filename(tmp_dir = tmp_dir)
    file(fasta_file,'w').writelines(fasta)

    # web version of blat
    #  blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 database.2bit query.fa output.psl
    # from: http://http://genome.ucsc.edu/FAQ/FAQblat.html
    #
    # other idea: 	./blat -minIdentity=95 –fine -stepSize=1 –tileSize=6 -repMatch = 1000000
    # from http://www.gene2drug.com/product/?p=671 by  Sucheta Tripathy
    # BLAT stands for Blast like Alignment Tool and was designed by Jim Kent.
    # It is relatively easy to install the software and is an excellent way of
    # mapping assembly files generated from abyss into reference genome for
    # finding new transcripts and new intron exon junctions. BLAT has recently
    # added few fine tuning options for short read mapping. Setting the stepSize
    # and tileSize parameters for mapping reads of length n, where
    # n = 2 * stepSize + tileSize – 1. While tileSize can range from 6 to 15,
    # stepSize can be from 1 to tileSize. So, in other words, reads as short
    # as 7 bases can be mapped into reference [2 * 1 + 6 – 1 = 7]. Also other
    # commandline options can be used to make the mapping more sensitive such
    # as –fine and -repMatch = 1000000. –fastmap option and –ooc option should
    # be avoided for mapping short reads. In addition –minIdentity may be set
    # to 95%.
    
    _BT_ = ""
    if blat_dir and blat_dir.strip():
        _BT_ = blat_dir.rstrip("/")+"/"
    
    cmd = None
    if align_type == 'web':
        cmd = [_BT_+'blat',
               '-stepSize=5', # 5
               '-repMatch=2253', # 2253
               '-minScore=0',  # 0
               '-minIdentity=0', # 0
               twobit,
               fasta_file,
               psl_file]
    elif align_type == 'sensitive':
        cmd = [_BT_+'blat',
               '-stepSize=5', # 5
               '-repMatch=2253', # 2253
               '-minScore=0',  # 0
               '-minIdentity=95', # 0
               '-fine',
               twobit,
               fasta_file,
               psl_file]
    else:
        print "ERROR: Not known type of BLAT search!"
        sys.exit(1)
    cmd = ' '.join(cmd)
    proc = os.system(cmd)
    if proc:
        print >>sys.stderr, "ERROR while executing '%s'" % (cmd,)
        sys.exit(1)

    psl = file(psl_file,'r').readlines()
    # add chr to the column number 14 (index = 1) so that can be loaded into UCSC
    # genome browser
    chr_psl = []
    for line in psl:
        li = line.split('\t')
        if li and len(li) > 14:
            l = li[13]
            li[13] = 'chr' + l if not l.startswith('chr') else l
            if li[13] == 'chrMT':
                li[13] = 'chrM'
            chr_psl.append('\t'.join(li))
        else:
            chr_psl.append(line)

    delete_file(fasta_file)
    delete_file(psl_file)
    return chr_psl



def give_me_sam(fastq, anchor, bowtie2index, bowtie2_dir = None, tmp_dir = None, cpus = 1):
    # give as input a file as a list of strings it runs BOWTIE2 and it returns
    # the SAM output as a list of strings
    fastq_file = give_me_temp_filename(tmp_dir = tmp_dir)
    sam_file = give_me_temp_filename(tmp_dir = tmp_dir)
    file(fastq_file,'w').writelines(fastq)

    _B2_ = ""
    if bowtie2_dir and bowtie2_dir.strip():
        _B2_ = bowtie2_dir.rstrip("/")+"/"

    cmd = [_B2_+'bowtie2',
           '-p',str(cpus),
           '--local',
           '-k','10',
           '-L',str(anchor),
           '-x',bowtie2index,
           '-U',fastq_file,
           '-S',sam_file]

    sam = []
    cmd = ' '.join(cmd)
    proc = os.system(cmd)
    if proc:
        print >>sys.stderr, "WARNING: unable to execute: '%s'" % (cmd,)
    else:
        sam = file(sam_file,'r').readlines()


    delete_file(fastq_file)
    delete_file(sam_file)
    return sam

def trimleft(t,s):
    #t="105856140,105856216,105863813,105863848,105863892,105863898"
    #s = 40
    d = map(int,t.split(","))
    n = len(d)
    p = 0
    f = []
    # go left
    for i in xrange(n-1,0,-2):
        #print d[i-1],d[i]
        l = abs(d[i]-d[i-1])
        if p + l > s:
            f.append(d[i])
            break
        else:
            p = p + l
            f.append(d[i])
            f.append(d[i-1])
    f.reverse()
    z = ",".join(map(str,f))
    return z
    
def trimright(t,s):
    #t="105856140,105856216,105863813,105863848,105863892,105863898"
    #s = 40
    d = map(int,t.split(","))
    n = len(d)
    p = 0
    f = []
    # go left
    for i in xrange(0,n,2):
        #print d[i],d[i+1]
        l = abs(d[i+1]-d[i])
        if p + l > s:
            f.append(d[i])
            break
        else:
            p = p + l
            f.append(d[i])
            f.append(d[i+1])
    z = ",".join(map(str,f))
    return z

def trimmycols(t,s):
    left = trimleft(t[4],s)
    right = trimright(t[11],s)
    z = t[0:4] + [left] + t[6:10] + [right]
    return z

def strtrimmycols(t,s):
    return '\t'.join(trimmycols(t,s))

def mycols(t):
    return t[0:4] + t[5:11]

def strmycols(t):
    return '\t'.join(mycols(t))

def mygenes(t):
    return tuple(t[0],t[6])

def ordmygenes(t):
    return myorder(mygenes(t))

def give_me_assembly(fasta, kmer = 31, velvet_dir = None, tmp_dir = None):
    # use Velvet to assembl the supporting reads
    #
    # velveth /tmp/velvet-unmapped-reads/ 17 -fasta -short myfasta.fa
    # velvetg /tmp/velvet-unmapped-reads/
    if fasta:
        fasta_file = give_me_temp_filename(tmp_dir = tmp_dir)
        ase_dir = give_me_temp_filename(tmp_dir = tmp_dir)
        if os.path.isfile(ase_dir) or os.path.islink(ase_dir):
            os.remove(ase_dir)
        elif os.path.isdir(ase_dir):
            os.rmtree(ase_dir)
        os.makedirs(ase_dir)
        file(fasta_file,'w').writelines(fasta)
        
        _VT_ = ""
        if velvet_dir and velvet_dir.strip():
            _VT_ = velvet_dir.rstrip("/")+"/"
        
        cmd = [_VT_+'velveth',
               ase_dir,
               str(kmer),
               '-fasta',
               '-short',
               fasta_file,
               ';',
               _VT_+'velvetg',
               ase_dir,
               '>',
               '/dev/null',
               '2>&1'
               ]
        cmd = ' '.join(cmd)
        proc = os.system(cmd)
        if proc:
            print >>sys.stderr, "ERROR while executing '%s'" % (cmd,)
            sys.exit(1)
    else:
        return []

    ase = file(os.path.join(ase_dir,'contigs.fa'),'r').readlines()
    delete_file(fasta_file)
    shutil.rmtree(ase_dir)
    return ase

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It analyzes the mappings of reads on exon-exon junctions."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input_fastq",
                      action = "store",
                      type = "string",
                      dest = "input_fastq_filename",
                      help = """The input FASTQ file containing all the reads.""")

    parser.add_option("--input_fusion_psl",
                      action = "store",
                      type = "string",
                      dest = "input_fusion_psl_filename",
                      help = """The input PSL file containing the candidate fusion genes.""")

    parser.add_option("--input_candidate_fusion_genes_reads",
                      action = "store",
                      type = "string",
                      dest = "input_candidate_fusion_genes_reads_filename",
                      help = """The input list of candidate fusion genes and ids of the supporting reads, for example 'candidate_fusion-genes_not-filtered_supporting_paired-reads.txt'. This is processed even further.""")


    parser.add_option("--input_unmapped_reads",
                      action = "store",
                      type = "string",
                      dest = "input_unmapped_reads_filename",
                      help = """The input list of ids of reads that are unmapped (that are mapping over the fusion junction).""")


    parser.add_option("--output_super_summary",
                      action = "store",
                      type = "string",
                      dest = "output_super_summary_filename",
                      help = """The output super summary report for candidate fusion genes.""")

    parser.add_option("--output_zip_fasta",
                      action = "store",
                      type = "string",
                      dest = "output_zip_fasta_filename",
                      help = """The ouput FASTQ file containing the reads which support each candidate fusion gene.""")

    parser.add_option("--suporting_unique_reads",
                      action = "store",
                      type = "int",
                      dest = "supporting_unique_reads",
                      default = 1,
                      help = """The minimum number of unique reads which overlap over an exon-exon junction. Default is %default.""")


    parser.add_option("--anchor2",
                      action = "store",
                      type = "int",
                      dest = "anchor2",
                      default = 150,
                      help = """For anchors longer (or equal) with this value it is enough to have only one supporting read. Default is '%default'.""")

    parser.add_option("--trim-complex",
                      action = "store",
                      type = "int",
                      dest = "complex",
                      default = 40,
                      help = """This value will be used for trimming at the fusion junction left adn right in order to be able to match the PSL mappings when they contain several exon-introns. Default is '%default'.""")

    parser.add_option("--input_genome_2bit",
                      action = "store",
                      type = "string",
                      dest = "input_genome_2bit",
                      help = """Path to the genome in 2bit format (generated with faToTwoBit) which will be used for aligning using BLAT the supporting reads and their alignment in PSL format is added to file specified with '--output_zip_fasta'.""")

    parser.add_option("--input_genome_bowtie2",
                      action = "store",
                      type = "string",
                      dest = "input_genome_bowtie2",
                      help = """Path to the genome in BOWTIE2 index format which will be used for aligning using BOWTIE2 the supporting reads and their alignment in PSL format is added to file specified with '--output_zip_fasta'.""")


    choices = ('web','sensitive')
    parser.add_option("--psl_alignment_type",
                      action = "store",
                      type = "choice",
                      choices = choices,
                      dest = "psl_search_type",
                      default = "web",
                      help = "The type of BLAT alignment to be used for aligning "+
                             "the supporting reads when BLAT is chosen. The choices "+
                             "are ['"+"','".join(choices)+"']. "+
                             "Default is '%default'.")

    parser.add_option("--blat-dir",
                      action = "store",
                      type = "string",
                      dest = "blat_directory",
                      help = """Path to Blat's executable.""")


    parser.add_option("--sam_alignment",
                      action = "store",
                      type = "int",
                      dest = "sam_alignment",
                      default = 10,
                      help = """If set then a SAM file will be generated using BOWTIE2. Default is '%default'.""")

    parser.add_option("--bowtie2-dir",
                      action = "store",
                      type = "string",
                      dest = "bowtie2_directory",
                      help = """Path to Bowtie2's executable.""")

    parser.add_option("--mismatches",
                      action = "store",
                      type = "int",
                      dest = "mismatches",
                      default = 3,
                      help = """The minimum number of mismatches accepted in the alignment. Default is '%default'.""")

    parser.add_option("--mismatches-gap",
                      action = "store",
                      type = "int",
                      dest = "mismatches_gap",
                      default = 30,
                      help = """The minimum number of mismatches accepted in the gap alignment. Default is '%default'.""")


    parser.add_option("--junction",
                      action = "store_true",
                      dest = "junction",
                      default = False,
                      help = """If used then the junction sequence is added to the FASTA file with the supporting reads. Default is '%default'.""")

    parser.add_option("--threads","-p",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 1,
                      help = "Number or processes to be used for running Bowtie2. "+
                             "Default is '%default'. ")

    parser.add_option("--tmp_dir",'-t',
                      action = "store",
                      type = "string",
                      dest = "tmp_directory",
                      default = None,
                      help = "The directory which should be used as temporary directory. By default is the OS temporary directory.")

    parser.add_option("--velvet",
                      action = "store_true",
                      dest = "velvet",
                      default = False,
                      help = """If used then the supporting reads from the FASTA file are assembled using VELVET. Default is '%default'.""")

    parser.add_option("--velvet-dir",
                      action = "store",
                      type = "string",
                      dest = "velvet_directory",
                      help = """Path to Velvet's executable.""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_fusion_psl_filename and
            options.output_super_summary_filename
            ):
        parser.print_help()
        sys.exit(1)

    #
    # HEADER PSL file
    #
    #header = ['gene-5end',                                           # 0
    #          'gene-5end_symbol',                                    # 1
    #          'chromosome_gene-5end',                                # 2
    #          'strand_gene-5end',                                    # 3
    #          'start_chromosome_part-1-of-read-mapped-gene-5end',    # 4
    #          'end_chromosome_part-1-read-mapped-gene-5end',         # 5
    #          'gene-3end',                                           # 6
    #          'gene-3end_symbol',                                    # 7
    #          'chromosome_gene-3end',                                # 8
    #          'strand_gene-3end',                                    # 9
    #          'start_chromosome_part-2-of-read-mapped-gene-3end',    # 10
    #          'end_chromosome_part-2-read-mapped-gene-3end',         # 11
    #          'short_read',                                          # 12
    #          'mismatches',                                          # 13
    #          'length_short_read',                                   # 14
    #          'start_part-1-read_on_gene-5end',                      # 15
    #          'end_part-1-read_on_gene-5end',                        # 16
    #          'start_part-2-read_on_gene-3end',                      # 17
    #          'end_part-2-read_on_gene-3end',                        # 18
    #          'anchor_length',                                       # 19
    #          'fusion_sequence'                                      # 20
    #          ]


    unmapped_reads = set()
    if options.input_unmapped_reads_filename:
        print "Reading...",options.input_unmapped_reads_filename
        unmapped_reads = set([line.rstrip('\r\n') for line in file(options.input_unmapped_reads_filename,'r').readlines()])


    print "Reading...",options.input_fusion_psl_filename
    data = [line.rstrip('\r\n').split('\t') for line in file(options.input_fusion_psl_filename,'r') if line.rstrip('\r\n')]
    header = data.pop(0)


    # filter for mismatches
    #data = [line for line in data if int(line[13])<=options.mismatches]
    dudu = []
    for line in data:
        if int(line[13])<=options.mismatches_gap:
            dudu.append(line)
            # here I have gaps in alignment (i.e. IGH fusions) because there is "*NNNN"
    #data = [line for line in data if int(line[13])<=options.mismatches] # ORIGINAL
    data = dudu

    comp = options.complex

    # find unique reads

    data_uniq = list(set(['\t'.join(line[:12]) for line in data]))
#    print "1 -> ",data_uniq
    data_uniq = [line.split('\t') for line in data_uniq]
    # find same splicing sites = remove cols 4 and 11
    data_uniq = [strtrimmycols(line,comp) for line in data_uniq]
    # counts the unique reads for unique splicing sites
    data_dict = dict()
    for line in data_uniq:
        data_dict[line] = data_dict.get(line,0) + 1
    # sort the counts

    dd = sorted(data_dict.items(),key = lambda x: -x[1])

    # filter those fusion with too few counts
    #dd = [(k,v) for (k,v) in dd if v >= options.supporting_unique_reads]
    dd = [(k,v) for (k,v) in dd if v >= 1] # in order to allow the use of options.anchor2


    # find those reads and the fusion sequence for the unique fusion points
    summary = []
    summary_reads = []
    summary.append("%s\tcounts\tlongest_anchor\tfusion_sequence\n"%(strmycols(header),))
    summary_reads.append('header')
    singles = set()
    ggenes_e = list()
    ggenes_s = list()
    ggenes_p = list()
    ggenes_e.append('header')
    ggenes_s.append('header')
    ggenes_p.append('header')
    fast_gg = dict()
    for (k,v) in dd:
        if v >= options.supporting_unique_reads:
            r = []
            fs = None
            gg = None
            anchor_max = 0
            kk = ""
            for li in data:
                if strtrimmycols(li,comp) == k:
                    r.append(li[12])
                    fs = li[20]
                    anchor = int(li[19])
                    if anchor > anchor_max:
                        anchor_max = anchor
                    gg_e = (li[0],li[6])
                    gg_s = (li[1],li[7])
                    gg_p = (li[5],li[10])
                    if not kk:
                        kk = strmycols(li)
            summary.append("%s\t%d\t%d\t%s\n"%(kk,v,anchor_max,fs))
            #print k
            r = set(r)
            rr = set(el[:-1] + '1' if el.endswith('/2') else el[:-1]+'2' for el in r)
            summary_reads.append(list(r)+list(rr))
            singles.update(r)
            singles.update(rr)
            ggenes_e.append(gg_e)
            ggenes_s.append(gg_s)
            ggenes_p.append(gg_p)
            fast_gg[myorder(*gg_e)] = None
        elif options.supporting_unique_reads > 1:
            r = []
            fs = None
            gg = None
            anchor_max = 0
            kk = ""
            for li in data:
                if strtrimmycols(li,comp) == k:
                    r.append(li[12])
                    fs = li[20]
                    anchor = int(li[19])
                    if anchor > anchor_max:
                        anchor_max = anchor
                    gg_e = (li[0],li[6])
                    gg_s = (li[1],li[7])
                    gg_p = (li[5],li[10])
                    if not kk:
                        kk = strmycols(li)
            if anchor_max >= options.anchor2:
                summary.append("%s\t%d\t%d\t%s\n"%(kk,v,anchor_max,fs))
                r = set(r)
                rr = set(el[:-1] + '1' if el.endswith('/2') else el[:-1]+'2' for el in r)
                summary_reads.append(list(r)+list(rr))
                singles.update(r)
                singles.update(rr)
                ggenes_e.append(gg_e)
                ggenes_s.append(gg_s)
                ggenes_p.append(gg_p)
                fast_gg[myorder(*gg_e)] = None
    print "Writing the summary file...", options.output_super_summary_filename
    file(options.output_super_summary_filename,'w').writelines(summary)

#    print "--------------------------"
#    print summary
#    print "--------------------------"

    print "Reading...",options.input_candidate_fusion_genes_reads_filename
    # 0  - Fusion_gene_symbol_1
    # 1  - Fusion_gene_symbol_2
    # 2  - Fusion_gene_1
    # 3  - Fusion_gene_2
    # 4  - Count_paired-end_reads
    # 5  - Supporting_paired-read_ids  ==> separated by commas, e.g.  F1000018349733,F1000033997513,F1000046358541,F1000034322437,...
    candidate_fusions_reads = [line.rstrip('\r\n').split('\t') for line in file(options.input_candidate_fusion_genes_reads_filename,'r').readlines()]
    candidate_fusions_reads.pop(0) # remove the header
    candidate_fusions_reads = dict([(myorder(el[2],el[3]),el[5].split(',')) for el in candidate_fusions_reads if fast_gg.has_key(myorder(el[2],el[3]))])

    #
    # for each candidate fusion genes build a FASTA file containing the reads which support it
    #
    print "Processing the supporting reads..."
    fasta = dict()
    fastq = dict()
    pairs = dict()
    for (k,v) in candidate_fusions_reads.items():
        pairs[k] = []
        for vv in v:
            s1 = '%s/1' % (vv,)
            s2 = '%s/2' % (vv,)
            fasta[s1] = None
            fasta[s2] = None
            pairs[k].append(s1)
            pairs[k].append(s2)

    for k in singles:
        fasta[k] = None


    print "Scanning the FASTQ file...",options.input_fastq_filename
    for a_read in reads_from_fastq_file(options.input_fastq_filename):
        if fasta.has_key(a_read[0]):
            ev = a_read[0]
            w = a_read[1]
            q = a_read[2]
            if ev.endswith("/1"):
                fasta[ev] = w
            elif ev.endswith("/2"):
                fasta[ev] = dnaReverseComplement(w)
            fastq[ev] = (w,q)

    # create a ZIP FASTA file where is a file for each candidate fusion gene
    print "Writing the FASTA/FASTQ files containing the supporting reads...",options.output_zip_fasta_filename
    archive = zipfile.ZipFile(options.output_zip_fasta_filename, 'w', zipfile.ZIP_STORED, allowZip64 = True)
    for i in xrange(len(summary)):
        if i == 0: # skip header
            continue
        # for each candidate fusion
        #gg = "%s:%s_%s:%s_%s:%s" % (ggenes_s[i][0],ggenes_s[i][1],ggenes_e[i][0],ggenes_e[i][1],ggenes_p[i][0],ggenes_p[i][1])
        gg = "%s--%s__%s--%s" % (ggenes_s[i][0],ggenes_s[i][1],ggenes_p[i][0],ggenes_p[i][1])
        gg_e = myorder(ggenes_e[i][0],ggenes_e[i][1])

        da = []
        if options.junction:
            u = summary[i].rstrip('\n').split('\t')
            da = ['>JUNCTION__%s\n%s\n' % ('_'.join(u[:-1]),u[-1])]
        # write the junction sequence
        #archive.writestr("%s_junction.fa" % (gg,), da)
        # PSL
        #if options.input_genome_2bit:
        #    psl = give_me_psl(da,
        #                      options.input_genome_2bit,
        #                      tmp_dir = options.tmp_directory,
        #                      align_type = options.blat_search_type)
        #    archive.writestr("%s_junction.psl" % (gg,), ''.join(psl))


        # write the reads in FASTA file
        #da = []
        for v in sorted(summary_reads[i]):
            da.append(">%s_supports_fusion_junction\n"%(v,))
            da.append("%s\n"%(fasta[v],))
        for v in sorted(pairs[gg_e]):
            da.append(">%s_supports_fusion_pair\n"%(v,))
            da.append("%s\n"%(fasta[v],))
        archive.writestr("%s_reads.fa" % (gg,), ''.join(da))
        # PSL
        if options.input_genome_2bit:
            psl = give_me_psl(da,
                              options.input_genome_2bit,
                              blat_dir = options.blat_directory,
                              tmp_dir = options.tmp_directory,
                              align_type = options.psl_search_type)
            archive.writestr("%s_reads.psl" % (gg,), ''.join(psl))
        # VELVET
        if options.velvet:
            ase = give_me_assembly(da,
                                   17,
                                   velvet_dir = options.velvet_directory,
                                   tmp_dir = options.tmp_directory)
            archive.writestr("%s_assembly.fa" % (gg,), ''.join(ase))
        # write the reads in FASTQ file
        da = []
        for v in sorted(summary_reads[i]):
            da.append("@%s_supports_fusion_junction%s\n"%(v[:-2],v[-2:]))
            sq = fastq[v]
            da.append("%s\n+\n%s\n"%(sq[0],sq[1]))
        for v in sorted(pairs[gg_e]):
            da.append("@%s_supports_fusion_pair%s\n"%(v[:-2],v[-2:]))
            sq = fastq[v]
            #da.append("%s\n+\n%s\n"%(sq[0],illumina2sanger(sq[1])))
            da.append("%s\n+\n%s\n"%(sq[0],sq[1]))
        archive.writestr("%s_reads.fq" % (gg,), ''.join(da))
        if options.input_genome_bowtie2:
            sam = give_me_sam(da,
                              options.sam_alignment,
                              options.input_genome_bowtie2,
                              bowtie2_dir = options.bowtie2_directory,
                              tmp_dir = options.tmp_directory,
                              cpus = options.processes)
            archive.writestr("%s_reads.sam" % (gg,), ''.join(sam))


        # Ensembl ids of genes
        archive.writestr("%s_ensembl_ids.txt" % (gg,), '%s\n%s\n' % (ggenes_e[i][0],ggenes_e[i][1]))

    archive.close()              
