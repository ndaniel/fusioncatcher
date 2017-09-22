#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It produces a report with the summary of the fusion genes found. Also
FASTQ and FASTA files containing the supporting reads corresponding to each
fusion gene is generated.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2017 Daniel Nicorici

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

This file is not running/executing/using BLAT.

"""
import sys
import os
import optparse
import gc
import string
import zipfile
import Bio.SeqIO
import Bio.SeqIO.QualityIO
import datetime
import tempfile
import shutil
import gzip

ttable = string.maketrans("ACGTYRSWKMBDHV-","TGCARYSWMKVHDB-")

mapping_solexa2sanger = "".join([chr(0) for ascii in range(0, 59)]
                        + [chr(33 + int(round(Bio.SeqIO.QualityIO.phred_quality_from_solexa(q)))) for q in range(-5, 62 + 1)]
                        + [chr(0) for ascii in range(127, 256)])

mapping_illumina2sanger = "".join([chr(0) for ascii in range(0, 64)]
                          + [chr(33 + q) for q in range(0, 62 + 1)]
                          + [chr(0) for ascii in range(127, 256)])

empty_zip_data = 'PK\x05\x06\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'

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
    psl = []
    cmd = ' '.join(cmd)
    proc = os.system(cmd)
    if proc:
        print >>sys.stderr, "WARNING: unable to execute: '%s'" % (cmd,)
    else:
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


def pos_junction(t):
    su = [elx.split('junction=')[1] for elx in t.lower().rstrip('\r\n').split(';') if elx.startswith('junction=')]
    return int(su[0])

def star_junction(t,po):
    return "%s*%s" % (t[:po],t[po:])

def give_me_assembly(fasta, kmer = 31, velvet_dir = None ,tmp_dir = None):
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
            shutil.rmtree(ase_dir)
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
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input_fastq",
                      action = "store",
                      type = "string",
                      dest = "input_fastq_filename",
                      help = """The input FASTQ file containing all the reads (can be given as gzipped file too).""")

    parser.add_option("--input_fasta_juncs",
                      action = "store",
                      type = "string",
                      dest = "input_fasta_juncs_filename",
                      help = """The input FASTA file containing all the sequences of the exon-exon junctions.""")

    parser.add_option("--input_fusion_summary",
                      action = "store",
                      type = "string",
                      dest = "input_fusion_summary_filename",
                      help = """The input summary for fusion genes which is processed further, for example 'candidate_fusion-genes_exon-exon-junctions_summary.txt'.""")

    parser.add_option("--input_fusion_summary_more",
                      action = "store",
                      type = "string",
                      dest = "input_fusion_summary_more_filename",
                      help = """The input summary for fusion genes with even more information, for example 'candidate_fusion-genes_exon-exon-junctions_reads-positions.txt'. This is processed even further.""")

    parser.add_option("--input_candidate_fusion_genes",
                      action = "store",
                      type = "string",
                      dest = "input_candidate_fusion_genes_filename",
                      help = """The input list of candidate fusion genes, for example 'candidate_fusion-genes_no-offending-reads_label-no-proteins-paralogs-readthrough-similar-pseudogenes_further.txt'. This is processed even further.""")

    parser.add_option("--input_candidate_fusion_genes_reads",
                      action = "store",
                      type = "string",
                      dest = "input_candidate_fusion_genes_reads_filename",
                      help = """The input list of candidate fusion genes and ids of the supporting reads, for example 'candidate_fusion-genes_not-filtered_supporting_paired-reads.txt'. This is processed even further.""")

    parser.add_option("--input_candidate_fusions_missing_mates",
                      action = "store",
                      type = "string",
                      dest = "input_candidate_fusions_missing_mates_filename",
                      help = """The input list mate reads from pairs of reads together with their mappings on the genes, for example 'candidate_fusion-genes_missing_mates.txt'.""")

    parser.add_option("--input_exons",
                      action="store",
                      type="string",
                      dest="input_exons_filename",
                      help="""Database with exons position on chromosomes, e.g. 'more_exons_ensembl.txt'. This is used for filtering the UTRs extensions by removing any extension which overlaps with any exons from the database. This is optional.""")

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
                      help = """The minimum number of unique reads which overlap over an exon-exon junction. Default is '%default'.""")

    parser.add_option("--anchor2",
                      action = "store",
                      type = "int",
                      dest = "anchor2",
                      default = 40,
                      help = """For anchors longer (or equal) with this value it is enough to have only one supporting read. Default is '%default'.""")


    parser.add_option("--input_genome_2bit",
                      action = "store",
                      type = "string",
                      dest = "input_genome_2bit",
                      help = """Path to the genome in 2bit format (generated with faToTwoBit) which will be used for aligning using BLAT the supporting reads and their alignment in PSL format is added to file specified with '--output_zip_fasta'.""")

    parser.add_option("--blat-dir",
                      action = "store",
                      type = "string",
                      dest = "blat_directory",
                      help = """Path to Blat's executable.""")

    parser.add_option("--input_genome_bowtie2",
                      action = "store",
                      type = "string",
                      dest = "input_genome_bowtie2",
                      help = """Path to the genome BOWTIE2 index which is used to generate the alignments in SAM format which is added to file specified with '--output_zip_fasta'.""")

    parser.add_option("--bowtie2-dir",
                      action = "store",
                      type = "string",
                      dest = "bowtie2_directory",
                      help = """Path to Bowtie2's executable.""")

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

    parser.add_option("--sam_alignment",
                      action = "store",
                      type = "int",
                      dest = "sam_alignment",
                      default = 10,
                      help = """If set then a SAM file will be generated using BOWTIE2. Default is '%default'.""")


    parser.add_option("--junction",
                      action = "store_true",
                      dest = "junction",
                      default = False,
                      help = """If used then the junction sequence is added to the FASTA file with the supporting reads. Default is '%default'.""")

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

    parser.add_option("--output_all_candidate_fusion_genes_reads",
                      action = "store",
                      type = "string",
                      dest = "output_all_candidate_fusion_genes_reads_filename",
                      help = """The output list of candidate fusion genes and the supporting reads.""")



    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_fastq_filename and
            options.input_fasta_juncs_filename and
            options.input_fusion_summary_filename and
            options.input_fusion_summary_more_filename and
            options.input_candidate_fusion_genes_filename and
            options.output_super_summary_filename and
            options.input_candidate_fusion_genes_reads_filename and
            options.output_zip_fasta_filename and
            options.input_exons_filename and
            options.input_candidate_fusions_missing_mates_filename and
            options.output_super_summary_filename
            ):
        if options.output_super_summary_filename:
            # print only the header
            line = [
                "Fusion_gene_symbol_1(5end_partner)",            # 0
                "Fusion_gene_symbol_2(3end_partner)",            # 1
                "Fusion_gene_1(5end_partner)",                   # 2
                "Fusion_gene_2(3end_partner)",                   # 3
                "Fusion_transcript_1(5end_partner)",             # 4
                "Fusion_transcript_2(3end_partner)",             # 5
                "Fusion_exon_1(5end_partner)",                   # 6
                "Fusion_exon_2(3end_partner)",                   # 7
                "Fusion_exon_1_exon-location",                   # 8
                "Fusion_exon_2_exon-location",                   # 9
                "Fusion_gene_1_position(5end_partner)",          # 10
                "Fusion_gene_2_position(3end_partner)",          # 11
                "Spanning_pairs",                                # 12
                "Spanning_unique_reads",                         # 13
                "Longest_spanning_read",                         # 14
                "Fusion_sequence",                               # 15
                "Aligner(s)"                                     # 16
                    ]
            da = [(line[0],
                   line[1],
                   line[2],
                   line[3],
                   line[6],
                   line[7],
                   line[10],
                   line[11],
                   line[12],
                   line[13],
                   line[14],
                   line[16],
                   line[15])]
            file(options.output_super_summary_filename,"w").writelines(['\t'.join(map(str,line))+'\n' for line in da])

            if options.output_zip_fasta_filename:
                file(options.output_zip_fasta_filename,'wb').write(empty_zip_data)

            sys.exit(0)
        else:
            parser.print_help()
            parser.error("One of the options has not been specified.")
            sys.exit(1)


    print "Reading known exons...",options.input_exons_filename
    # 0  - ensembl_peptide_id
    # 1  - ensembl_gene_id
    # 2  - ensembl_transcript_id
    # 3  - ensembl_exon_id
    # 4  - exon_chrom_start
    # 5  - exon_chrom_end
    # 6  - rank
    # 7  - start_position
    # 8  - end_position
    # 9  - transcript_start
    # 10 - transcript_end
    # 11 - strand
    # 12 - chromosome_name
    # 13 - cds_start
    # 14 - cds_end
    # 15 - 5_utr_start
    # 16 - 5_utr_end
    # 17 - 3_utr_start
    # 18 - 3_utr_end
    gc.disable()
    exon = [line.rstrip('\r\n').split('\t') for line in file(options.input_exons_filename,'r').readlines() if line.rstrip('\r\n')]
    exon = dict([ (line[3], ("%s:%s-%s:%s" % (line[12],line[4],line[5],"+" if line[11]=="1" else "-"), # start & end positions
                             "%s:%s:%s" % (line[12],line[4] if line[11]=="1" else line[5],"+" if line[11]=="1" else "-"), # 5 end position of exon
                             "%s:%s:%s" % (line[12],line[5] if line[11]=="1" else line[4],"+" if line[11]=="1" else "-") # 3 end position of exons
                            )) for line in exon]) # chrom:start-end:strand
    gc.enable()


    print "Reading...",options.input_fusion_summary_filename
    # 0  - 5end_gene
    # 1  - 3end_gene
    # 2  - 5end_transcript
    # 3  - 3end_transcript
    # 4  - 5end_exon
    # 5  - 3end_exon
    # 6  - 5end_rank_exon_transcript
    # 7  - 3end_rank_exon_transcript
    # 8  - counts_reads
    # 9  - counts_unique_reads
    # 10 - mean_lengths_of_shorter_overlaping_parts_of_reads_over_junction
    # 11 - longest_of_shorter_overlaping_parts_of_reads_over_junction
    # 12 - lengths_of_shorter_overlaping_parts_of_reads_over_junction
    # 13 - average_number_mismatches_per_read
    # 14 - 5end_gene_symbol
    # 15 - 3end_gene_symbol
    # 16 - .
    # 17 - .
    # 18 - .
    # 19 - length_exon-exon_junction
    # 20 - position_junction
    # 21 - index_sequence_exon-exon_junction
    # 22 - .
    # 23 - .
    # 24 - .
    # 25 - id_sequence_in_fasta_file
    # 26 --> added here the 'fusion position': 5-partner-gene--3end-exon-position <-> 3-partner-gene--5end-exon-position
    juncs_fasta = dict()
    fusion_summary = [line.rstrip('\r\n').split('\t') for line in file(options.input_fusion_summary_filename,'r').readlines()]
    fusion_summary.pop(0) # remove the header
    fusion_summary = [line + [ "%s===%s" % (exon[line[4]][2], exon[line[5]][1])] for line in fusion_summary]
    fusion_summary = sorted(fusion_summary, key = lambda x: (-int(x[9]), -float(x[11]), float(x[13]) ) ) # order by the unique counts
    fusion_summary = [line for line in fusion_summary if ((int(line[9]) >= options.supporting_unique_reads) or (int(line[11]) >= options.anchor2))]
    temp = list()
    uniq = set()
    # remove duplicate lines which represent the same fusion (based on fusion point with chromosomal coordinates)
    for line in fusion_summary:
        juncs_fasta[line[25]] = ''
        if line[26] not in uniq:
            uniq.add(line[26])
            temp.append(line)
    fusion_summary = temp


    print "Reading...",options.input_fasta_juncs_filename
    input_seq_iterator = Bio.SeqIO.parse(open(options.input_fasta_juncs_filename, "rU"), "fasta")
    for record in input_seq_iterator:
        if juncs_fasta.has_key(record.id):
            juncs_fasta[record.id] = str(record.seq)
    # now column 27 contains the sequence of the junction
    #fusion_summary = [line + [juncs_fasta[line[25]]] for line in fusion_summary]
    fusion_summary = [line + [juncs_fasta[line[25]]] for line in fusion_summary]
    # now column 28 contains the nucleotide before the fusion
    # now column 29 contains the nucleotide after the fusion
    #fusion_summary = [line + [line[27][len(line[27])/2-1], line[27][len(line[27])/2]] for line in fusion_summary]



    print "Reading...",options.input_fusion_summary_more_filename
    # 0  - ENSG00000014138
    # 1  - ENSG00000149798
    # 2  - ENST00000265465
    # 3  - ENST00000279249
    # 4  - ENSE00002479086
    # 5  - ENSE00000992560
    # 6  - 17-66,19-68,27-76 <- positions suporting reads
    # 7  - S000035555693/2,S000049560905/2,S000099140369/2 <- supporting read ids
    fusion_summary_more = [line.rstrip('\r\n').split('\t') for line in file(options.input_fusion_summary_more_filename,'r').readlines()]
    fusion_summary_more = dict([( (el[0],el[1],el[2],el[3],el[4],el[5]),
                                   {'loc':el[6].split(','),
                                    'ids':el[7].split(',')
                                   }
                                ) for el in fusion_summary_more])


    print "Reading...",options.input_candidate_fusion_genes_filename
    # 0  - Fusion_gene_1
    # 1  - Fusion_gene_2
    # 2  - Count_paired-end_reads
    # 3  - Fusion_gene_symbol_1
    # 4  - Fusion_gene_symbol_2
    # 5  - Information_fusion_genes
    # 6  - Analysis_status -> further or skipped
    candidate_fusions = [line.rstrip('\r\n').split('\t') for line in file(options.input_candidate_fusion_genes_filename,'r').readlines()]
    candidate_fusions.pop(0) # remove the header
    candidate_fusions = dict([( tuple(myorder(line[0], line[1])), line[2]) for line in candidate_fusions if line[6].lower() == 'further_analysis'])


    print "Reading...",options.input_candidate_fusion_genes_reads_filename
    # 0  - Fusion_gene_symbol_1
    # 1  - Fusion_gene_symbol_2
    # 2  - Fusion_gene_1
    # 3  - Fusion_gene_2
    # 4  - Count_paired-end_reads
    # 5  - Supporting_paired-read_ids  ==> separated by commas, e.g.  F1000018349733,F1000033997513,F1000046358541,F1000034322437,...
    candidate_fusions_reads = [line.rstrip('\r\n').split('\t') for line in file(options.input_candidate_fusion_genes_reads_filename,'r').readlines()]
    candidate_fusions_reads.pop(0) # remove the header
    if options.output_all_candidate_fusion_genes_reads_filename:
        all_candidate_fusions_reads = [(el[2]+"--"+el[3]+"__"+el[0]+"--"+el[1],el[5].split(',')) for el in candidate_fusions_reads]
        xlist = []
        for xfusion,xreads in all_candidate_fusions_reads:
            f = "%s.%s.txt" % (options.output_all_candidate_fusion_genes_reads_filename,xfusion)
            xlist.append(f)
            xr = []
            for xrs in xreads:
                xr.append("%s/1" %(xrs,))
                xr.append("%s/2" %(xrs,))
            xr = sorted(set(xr))
            file(f,'w').writelines([el+'\n' for el in xr])
        file(options.output_all_candidate_fusion_genes_reads_filename,'w').writelines([el+'\n' for el in xlist])
        all_candidate_fusions_reads = []
    candidate_fusions_reads = dict([(tuple(myorder(el[2],el[3])),el[5].split(',')) for el in candidate_fusions_reads if myorder(el[2],el[3]) in candidate_fusions])
    
    print "Reading...",options.input_candidate_fusions_missing_mates_filename
    # 0  - missing_mate
    # 1  - found_mate
    # 2  - genes_where_the_found_mate_maps
    missing = [line.rstrip('\r\n').split('\t') for line in file(options.input_candidate_fusions_missing_mates_filename,'r') if line.find(",")==-1]
    missing.pop(0) # remove the header
    missing = dict([(el[0],el[2]) for el in missing])


    ############################################################################
    # build the summary
    ############################################################################

    # 0  - Fusion_gene_symbol_1
    # 1  - Fusion_gene_symbol_2
    # 2  - Fusion_gene_1
    # 3  - Fusion_gene_2
    # 4  - Fusion_transcript_1
    # 5  - Fusion_transcript_2
    # 6  - Fusion_exon_1
    # 7  - Fusion_exon_2
    # 8  - Count_paired-end_reads
    # 9  - Count_unique_reads_exon-exon_junction
    # 10 - Max_overlapping_exon-exon_junction
    support = dict()
    support_pair = set()
    juncs_fasta = dict()
    data = []
    extra = dict()
    for line in fusion_summary:

        # search for info about this gene-gene in exon-exon junctions
        g1 = line[0] # gene id
        g2 = line[1] # gene id
        s1 = line[14] # gene symbol
        s2 = line[15] # gene symbol
        t1 = line[2] # transcript id
        t2 = line[3] # transcript id
        e1 = line[4] # exon id
        e2 = line[5] # exon id
        ep1 = exon[e1][0] # exon position
        ep2 = exon[e2][0] # exon position
        f1 = exon[e1][2]
        f2 = exon[e2][1]
        cp = int(candidate_fusions[myorder(line[0],line[1])]) # count pairs
        cu = int(line[9]) # count unique reads
        mo = int(line[11]) # longest_of_shorter_overlaping_parts_of_reads_over_junction


        # find and save the supporting pairs
        #key = "%s:%s_%s:%s_%s:%s" % (s1,s2,g1,g2,e1,e2)
        key =  "%s--%s__%s--%s" % (s1,s2,f1.split(':')[1],f2.split(':')[1])
        extra[key] = (e1,e2)
        if not support.has_key(key):
            support[key] = set()
            juncs_fasta[key] = set()
        pairs = candidate_fusions_reads.get((g1,g2),None)
        if not pairs:
            pairs = candidate_fusions_reads.get((g2,g1),None)
            if not pairs:
                print g1
                print g2
                print pairs
                print "Error: something wrong!"
                sys.exit(1)
        pairs = set(pairs)

        for pai in pairs:
            support[key].add("%s/1"%(pai,))
            support[key].add("%s/2"%(pai,))
            support_pair.add("%s/1"%(pai,))
            support_pair.add("%s/2"%(pai,))
        star = star_junction(line[27],pos_junction(line[25]))
        juncs_fasta[key].add(">JUNCTION__%s\n%s\n" % (line[25],star))
        # supporting junction-reads (and corresponding mate if the the mate maps
        # on one of the fusion genes to which the junction belongs)


        # find reads which overlap over exon-exon junction
        #we have a junction hit!
        juncs = fusion_summary_more.get((g1,g2,t1,t2,e1,e2),None)
        juncs_ids = set(juncs["ids"])
        support[key].update(juncs_ids)
        # add their mates too
        mate = 0
        for el in juncs_ids:
            support[key].add( el[:-1] + '1' if el.endswith('/2') else el[:-1]+'2')

            # count also their mate read if it maps on the this gene fusion
            res = missing.get(el,None)
            if res and ((res == g1) or (res == g2)):
                # found it => add it
                mate = mate + 1


        #
        temp = [s1,
                s2,
                g1,
                g2,
                t1,
                t2,
                e1,
                e2,
                ep1,
                ep2,
                f1,
                f2,
                cp,
                cu,
                mo,
                star,
                'BOWTIE'
                ]
        #
        data.append(temp)



    print "Writing report",options.output_super_summary_filename
    data = sorted(data, key = lambda x:( -x[13], -x[14], -x[12] ))
    data.insert(0, [
        "Gene_1_symbol(5end_fusion_partner)",            # 0
        "Gene_2_symbol(3end_fusion_partner)",          # 1
        "Gene_1_id(5end_fusion_partner)",                # 2
        "Gene_2_id(3end_fusion_partner)",                # 3
        "Transcript_1_id(5end_fusion_partner)",          # 4
        "Transcript_2_id(3end_fusion_partner)",          # 5
        "Exon_1_id(5end_fusion_partner)",                # 6
        "Exon_2_id(3end_fusion_partner)",                # 7
        "Exon_1_genomic_location(5end_fusion_partner)",  # 8
        "Exon_2_genomic_location(3end_fusion_partner)",  # 9
        "Fusion_point_for_gene_1(5end_fusion_partner)",  # 10
        "Fusion_point_for_gene_2(3end_fusion_partner)",  # 11
        "Spanning_pairs",                                # 12
        "Spanning_unique_reads",                         # 13
        "Longest_anchor_found",                          # 14
        "Fusion_sequence",                               # 15
        "Fusion_finding_method"                          # 16
    ])
    da = [(line[0],
           line[1],
           line[2],
           line[3],
           line[6],
           line[7],
           line[10],
           line[11],
           line[12],
           line[13],
           line[14],
           line[16],
           line[15]) for line in data]
    file(options.output_super_summary_filename,"w").writelines(['\t'.join(map(str,line))+'\n' for line in da])


    #
    # for each candidate fusion genes build a FASTA file containing the reads which support it
    #
    print "Processing the supporting reads..."
    # write the list of supporting reads (pair-reads, junction-reads)
    fasta = dict()
    fastq = dict()
    for (k,v) in support.items():
        for ev in v:
            fasta[ev] = None
            fastq[ev] = None


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
    for (gg,vv) in support.items():
        # for each candidate fusion

        # write the junction sequence
        da = []
        if options.junction:
            for v in juncs_fasta[gg]:
                da.append(v)
        #archive.writestr("%s_junction.fa" % (gg,), ''.join(da))
        # PSL
        #if options.input_genome_2bit:
        #    psl = give_me_psl(da,
        #                      options.input_genome_2bit,
        #                      tmp_dir = options.tmp_directory,
        #                      align_type = options.blat_search_type)
        #    archive.writestr("%s_junction.psl" % (gg,), ''.join(psl))


        # write the reads in FASTA file
        #da = []
        vvv = sorted(vv)
        for v in vvv:
            if v in support_pair:
                da.append(">%s_supports_fusion_pair\n"%(v,))
            else:
                da.append(">%s_supports_fusion_junction\n"%(v,))
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
        vvv = sorted(vv)
        for v in vvv:
            if v in support_pair:
                da.append("@%s_supports_fusion_pair%s\n"%(v[:-2],v[-2:]))
            else:
                da.append("@%s_supports_fusion_junction%s\n"%(v[:-2],v[-2:]))
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
        u = extra[gg]
        archive.writestr("%s_ensembl_ids.txt" % (gg,), '%s\n%s\n' % (u[0],u[1]))


    archive.close()
    #
