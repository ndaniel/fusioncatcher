#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Takes a list of genes with two genes per line and generates a FASTA with the two genes' sequences concatenated.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

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

import os
import sys
import optparse
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet
import itertools
import gc
import shutil
import random

def divergent(list_seq):
    """
    Given a list it splits the input list into several list such that each gene appears only once in one output list.
    """
    n_list_seq = len(list_seq)
    flag = [True] * n_list_seq
    now = set()
    while any(flag):
        block = []
        now = set()
        for i in xrange(n_list_seq):
            if flag[i]:
                ngg = list_seq[i].id.split('|')
                if (ngg[0] not in now) and (ngg[1] not in now):
                    block.append(list_seq[i])
                    flag[i] = False
                    now.add(ngg[0])
                    now.add(ngg[1])
        if block:
            yield block


def comb_all(some_id,
         gene_a,
         gene_b,
         some_strand):
    """
    It generates candidate junctions for TopHat, by doing all exon-exon combinations of the two genes (even inside the same gene)!

    Format of junctions for TopHat is:

    #<chrom> <left> <right> <+/->

    left and right are zero-based coordinates, and specify the last character of the left sequenced to be spliced to the first character of the right sequence, inclusive.
    """
    all_juncs = list(gene_a['exons'])
    r = gene_a['length']
    for (ja,jb) in gene_b['exons']:
        all_juncs.append((r+ja,r+jb))
    nj = len(all_juncs)
    result = set()
    for i in xrange(nj-1):
        for j in xrange(i+1,nj):
            a = all_juncs[i][0]
            b = all_juncs[i][1]
            c = all_juncs[j][0]
            d = all_juncs[j][1]
            if a > b or c > d:
                print "Error: unexpected junction position!"
                sys.exit(1)
            if b < c:
                result.add((some_id,b-1,c-1,some_strand))
    result = list(result)
    return result


def comb_fusion(some_id,
         gene_a,
         gene_b,
         some_strand):
    """
    It generates candidate junctions for TopHat, by doing all exon-exon combinations of the two genes (not inside the same gene)!

    Format of junctions for TopHat is:

    #<chrom> <left> <right> <+/->

    left and right are zero-based coordinates, and specify the last character of the left sequenced to be spliced to the first character of the right sequence, inclusive.
    """
    juncs_a = list(gene_a['exons'])
    r = gene_a['length']
    juncs_b = []
    for (ja,jb) in gene_b['exons']:
        juncs_b.append((r+ja,r+jb))
    nja = len(juncs_a)
    njb = len(juncs_b)
    result = set()
    for i in xrange(nja):
        for j in xrange(njb):
            a = juncs_a[i][0]
            b = juncs_a[i][1]
            c = juncs_b[j][0]
            d = juncs_b[j][1]
            if a > b or c > d:
                print "Error: unexpected junction position!"
                sys.exit(1)
            if b < c:
                result.add((some_id,b-1,c-1,some_strand))
    result = list(result)
    return result


def int2str(x,n=10):
    x = str(x)
    return '0' * int(n - len(x)) + x


#########################

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""Takes a list of genes with two genes per line and generates a FASTA with the two genes' sequences concatenated."""
    version="%prog 0.12 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The text separated file containing on each line two names of genes.""")

    parser.add_option("--input_database",
                      action="store",
                      type="string",
                      dest="input_database_filename",
                      help="""The FASTA file containg the sequences of all genes.""")

    parser.add_option("--input_exons",
                      action="store",
                      type="string",
                      dest="exons_filename",
                      help="""Database with exons position on chromosomes, e.g. 'more_exons_ensembl.txt'. This is used for filtering the UTRs extensions by removing any extension which overlaps with any exons from the database. This is optional.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""A FASTA file containing the sequences of the two genes per line concatenated.""")

    parser.add_option("--output_genes",
                      action="store",
                      type="string",
                      dest="output_genes_filename",
                      help="""A FASTA file containing the sequences of all genes from the input (a gene will appear only once in the output).""")

    parser.add_option("--output_dir",
                      action="store",
                      type="string",
                      dest="output_dir",
                      help="""An output directory containing FASTA files containing the sequences of the two genes per line concatenated. One Fasta file contains one sequence.""")

    parser.add_option("--output_tophat_juncs",
                      action="store",
                      type="string",
                      dest="output_tophat_juncs_filename",
                      help="""A junctions file which can be used further as input to TopHat.""")

    parser.add_option("--longest",
                      action="store",
                      type="string",
                      dest="output_longest",
                      help="""A text file where it will be written the size of the longest sequence of two genes which have been concatenated.""")

    parser.add_option("--reverse",
                      action="store_true",
                      dest="reverse",
                      default=False,
                      help="""If this is True then for a given set of two genes A and B two sequences will be generated for A+B and B+A. Default is '%default'.""")

    parser.add_option("--output_genes_count_nuc",
                      action="store",
                      type='string',
                      dest="output_genes_count_nuc_filename",
                      help="""If used then the number of nucleotides of all sequences from the output FASTA file (i.e. --output_genes) will be reported.""")

    parser.add_option("--output_genes_count_seq",
                      action="store",
                      type='string',
                      dest="output_genes_count_seq_filename",
                      help="""If used then the number of sequences from the output FASTA file (i.e. --output_genes) will be reported.""")

    parser.add_option("--padding",
                      action="store",
                      type="int",
                      dest="padding",
                      default=0,
                      help="""If this is larger than 0 than 'N' characters will be added at the end of each sequence from the output (it might help for avoiding cross-scaffolding). Default is '%default'.""")



    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.input_database_filename and
            (options.output_filename or options.output_dir)
            ):
        parser.print_help()
        parser.error('Missing arguments!')


    print "Reading the list of gene-gene...",options.input_filename
    genes = [line.rstrip('\r\n').split('\t') for line in file(options.input_filename,'r').readlines() if line.rstrip('\r\n')]
    set_genes = set(list(itertools.chain(*genes)))

    if options.exons_filename:
        print "Processing the exons database..."
        # ensembl_peptide_id             0
        # ensembl_gene_id                1
        # ensembl_transcript_id          2
        # ensembl_exon_id                3
        # exon_chrom_start               4
        # exon_chrom_end                 5
        # rank                           6
        # start_position                 7
        # end_position                   8
        # transcript_start               9
        # transcript_end                 10
        # strand                         11
        # chromosome_name                12
        # cds_start                      13
        # cds_end                        14
        # 5_utr_start                    15
        # 5_utr_end                      16
        # 3_utr_start                    17
        # 3_utr_end                      18
        exons = [line.rstrip('\r\n').split('\t') for line in file(options.exons_filename,'r').readlines() if line.rstrip('\r\n')]
        exons = [(line[1], # gene_id              0
                  line[7], # gene_start           1
                  line[8], # gene_end             2
                  line[4], # exon_start           3
                  line[5], # exon_end             4
                  line[11] # strand               5
                  ) for line in exons if line[1] in set_genes]
        database = {}
        for line in exons:
            gn = line[0]
            gs = int(line[1])
            ge = int(line[2])
            es = int(line[3])
            ee = int(line[4])
            st = int(line[5])
            if gs > ge:
                (gs,ge) = (ge,gs)
            if es > ee:
                (es,ee) = (ee,es)
            if gn not in database:
                database[gn] = {'length':0,
                                'strand':0,
                                'exons':set()}
            ps = 0
            pe = 0
            if st == 1:
                ps = es - gs + 1
                pe = ee - gs + 1
            elif st == -1:
                ps = ge - ee + 1
                pe = ge - es + 1
            else:
                print "Unknown strand!"
                sys.exit(1)
            database[gn]['length'] = ge - gs + 1
            database[gn]['strand'] = st
            database[gn]['exons'].add((ps,pe))
        print "Generating list of junctions..."
        #<chrom> <left> <right> <+/->      left and right are zero-based coordinates, and specify the last character of the left sequenced to be spliced to the first character of the right sequence, inclusive.
        junctions = []
        for line in genes:
            gene_1 = line[0]
            gene_2 = line[1]
            id = gene_1 + '|' + gene_2 + '|' + str(database[gene_1]['length'])
            tj = comb_fusion(id,database[gene_1],database[gene_2],'+')
            junctions.extend(['\t'.join(map(str,li))+'\n' for li in tj])
            if options.reverse:
                id = gene_2 + '|' + gene_1 + '|' + str(database[gene_2]['length'])
                tj = comb_fusion(id,database[gene_2],database[gene_1],'+')
                junctions.extend(['\t'.join(map(str,li))+'\n' for li in tj])
        if options.output_tophat_juncs_filename:
            print 'Writing the ',len(junctions),'junctions for TopHat...'
            file(options.output_tophat_juncs_filename,'w').writelines(junctions)

    print "Reading all genes' sequences...",options.input_database_filename
    input_seq_iterator = Bio.SeqIO.parse(open(options.input_database_filename, "rU"), "fasta")
    seq_dict = dict((record.id,str(record.seq)) for record in input_seq_iterator if record.id in set_genes)

    print "Writing..."
    random.seed(274876858367)
    sequences = []
    max_len = 0
    padding = options.padding
    i_shift = 0
    output_handle = open(options.output_filename, "w") if options.output_filename else None
    if options.output_dir:
        if not os.path.exists(options.output_dir):
            os.makedirs(options.output_dir)
        mylist = []
    for line in genes:
        gene_1 = line[0]
        gene_2 = line[1]
        if not seq_dict.has_key(gene_1):
            print "The sequence of gene",gene_1,"has not been found!"
            continue
        if not seq_dict.has_key(gene_2):
            print "The sequence of gene",gene_2,"has not been found!"
            continue
        seq_1 = seq_dict[gene_1]
        seq_2 = seq_dict[gene_2]
        id = gene_1 + '|' + gene_2 + '|' + str(len(seq_1))
        seq = seq_1 + seq_2
        if len(seq) > max_len:
            max_len = len(seq)
        if padding:
            seq = seq + 'A'*padding
        gc.disable()
        sequences.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq,Bio.Alphabet.IUPAC.ambiguous_dna),id=id,name="",description=""))
        gc.enable()
        if options.reverse:
            id = gene_2 + '|' + gene_1 + '|' + str(len(seq_2))
            seq = seq_2 + seq_1
            if padding:
                seq = seq + 'A'*padding
            gc.disable()
            sequences.insert((len(sequences)-1)/2,Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq,Bio.Alphabet.IUPAC.ambiguous_dna),id=id,name="",description=""))
            gc.disable()

        if options.output_filename and len(sequences) > 10**5:
            # I do this in order to avoid the A+B and B+A (or A+C) combination of two genes to be near each other
            sequences = tuple(sequences)
            ls = len(sequences)
            if ls > 4:
                sequences = tuple(random.sample(sequences,ls))
            Bio.SeqIO.write(sequences, output_handle, "fasta")
            sequences = []

        if options.output_dir:
            for seqs in divergent(sequences):
                file_name = os.path.join(options.output_dir,'gene-gene_'+int2str(i_shift)+'.fa')
                mylist.append(file_name)
                another_output_handle = open(file_name, "w")
                Bio.SeqIO.write(seqs, another_output_handle, "fasta")
                another_output_handle.close()
                i_shift = i_shift + 1

    if options.output_filename:
        if sequences:
            sequences = tuple(sequences)
            ls = len(sequences)
            if ls > 4:
                sequences = tuple(random.sample(sequences,ls))
            Bio.SeqIO.write(sequences, output_handle, "fasta")
            sequences = []
        output_handle.close()

    if options.output_dir:
        file(os.path.join(options.output_dir,'list_divergent_fasta_gene-gene.txt'),'w').writelines([line+'\n' for line in mylist])

    print "Length of longest sequence = ",max_len
    if options.output_longest:
        file(options.output_longest,'w').write(str(max_len))

    count_nuc = 0
    count_seq = 0
    if options.output_genes_filename:
        print "Writing the sequences of all genes from the input in a separate file (one gene will be written only once)..."
        sequences = []
        for k,v in seq_dict.iteritems():
            sequences.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(v,Bio.Alphabet.IUPAC.ambiguous_dna),id=k,name="",description=""))
            count_nuc = count_nuc + len(v)
        count_seq = len(sequences)
        output_handle = open(options.output_genes_filename, "w")
        Bio.SeqIO.write(sequences, output_handle, "fasta")
        output_handle.close()

    if options.output_genes_count_seq_filename:
        file(options.output_genes_count_seq_filename,"w").write("%d" % (count_seq,))

    if options.output_genes_count_nuc_filename:
        file(options.output_genes_count_nuc_filename,"w").write("%d" % (count_nuc,))


    print "The end."
    #
