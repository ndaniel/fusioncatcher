#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Given a PSL format file with alignments of contigs on genome it gives
a list of candidate fusion genes.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2022 Daniel Nicorici

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

#
"""
More about PSL format is here: http://genome.ucsc.edu/FAQ/FAQformat#format2


PSL format

PSL lines represent alignments, and are typically taken from files generated by
BLAT or psLayout. See the BLAT documentation for more details. All of the
following fields are required on each data line within a PSL file:

   1. matches - Number of bases that match that aren't repeats
   2. misMatches - Number of bases that don't match
   3. repMatches - Number of bases that match but are part of repeats
   4. nCount - Number of 'N' bases
   5. qNumInsert - Number of inserts in query
   6. qBaseInsert - Number of bases inserted in query
   7. tNumInsert - Number of inserts in target
   8. tBaseInsert - Number of bases inserted in target
   9. strand - '+' or '-' for query strand. For translated alignments,
      second '+' or '-' is for genomic strand
  10. qName - Query sequence name
  11. qSize - Query sequence size
  12. qStart - Alignment start position in query
  13. qEnd - Alignment end position in query
  14. tName - Target sequence name
  15. tSize - Target sequence size
  16. tStart - Alignment start position in target
  17. tEnd - Alignment end position in target
  18. blockCount - Number of blocks in the alignment (a block contains no gaps)
  19. blockSizes - Comma-separated list of sizes of each block
  20. qStarts - Comma-separated list of starting positions of each block in query
  21. tStarts - Comma-separated list of starting positions of each block in target

Example:
Here is an example of an annotation track in PSL format. Note that line breaks
have been inserted into the PSL lines in this example for documentation display
purposes. Click here for a copy of this example that can be pasted into the
browser without editing.

track name=fishBlats description="Fish BLAT" useScore=1
59 9 0 0 1 823 1 96 +- FS_CONTIG_48080_1 1955 171 1062 chr22
    47748585 13073589 13073753 2 48,20,  171,1042,  34674832,34674976,
59 7 0 0 1 55 1 55 +- FS_CONTIG_26780_1 2825 2456 2577 chr22
    47748585 13073626 13073747 2 21,45,  2456,2532,  34674838,34674914,
59 7 0 0 1 55 1 55 -+ FS_CONTIG_26780_1 2825 2455 2676 chr22
    47748585 13073727 13073848 2 45,21,  249,349,  13073727,13073827,

Be aware that the coordinates for a negative strand in a PSL line are handled
in a special way. In the qStart and qEnd fields, the coordinates indicate the
position where the query matches from the point of view of the forward strand,
even when the match is on the reverse strand. However, in the qStarts list,
the coordinates are reversed.

Example:
Here is a 30-mer containing 2 blocks that align on the minus strand and 2 blocks
that align on the plus strand (this sometimes can happen in response to
assembly errors):

0         1         2         3 tens position in query
0123456789012345678901234567890 ones position in query
            ++++          +++++ plus strand alignment on query
    --------    ----------      minus strand alignment on query

Plus strand:
     qStart=12
     qEnd=31
     blockSizes=4,5
     qStarts=12,26

Minus strand:
     qStart=4
     qEnd=26
     blockSizes=10,8
     qStarts=5,19

Essentially, the minus strand blockSizes and qStarts are what you would get if
you reverse-complemented the query. However, the qStart and qEnd are not
reversed. To convert one to the other:

     qStart = qSize - revQEnd
     qEnd = qSize - revQStart
"""

import os
import sys
import optparse
import tempfile
import Bio.SeqIO
import itertools

#########################
def coord_gene2genome(a_gene,positions,a_database):
    r = []
    sta = a_database[a_gene]['strand']
    s = a_database[a_gene]['start']
    e = a_database[a_gene]['end']
    if sta == 1:
        for ps in positions:
            r.append(s + ps - 1)
    elif sta == -1:
        for ps in positions:
            r.append(e - ps + 1)
    else:
        print "Unknown strand!", strand
        sys.exit(1)
    return r



#####################################
#####################################
#####################################
if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""Given a PSL format file with alignments of contigs on genome it gives
a candidate list of fusion genes (where short reads have been aligned using BLAT).
"""
    version="%prog 0.11 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_mappings",
                      action="store",
                      type="string",
                      dest="input_mappings_filename",
                      help="""The input file in PSL format containing the reads/contigs uniquely mapped on fusion genes.""")

    parser.add_option("--input_genes_positions",
                      action="store",
                      type="string",
                      dest="input_genes_positions_filename",
                      help="""A database containing the genes positions on the genome, e.g. 'ensembl/genes_positions_ensembl.txt'.""")

    parser.add_option("--input_genegene_fasta",
                      action="store",
                      type="string",
                      dest="input_genegene_fasta_filename",
                      help="""A FAST file containing the sequences of the gene-gene combinations used for finding fusion genes, e.g. 'gene-gene.fa'.""")

    parser.add_option("--input_hugo",
                      action="store",
                      type="string",
                      dest="input_hugo_filename",
                      help="""The input database used for linking ENSEMBL GENE ID to HUGO gene names, e.g. 'genes_info_ensembl.txt'.""")

    parser.add_option("--threshold_matches",
                      action="store",
                      type="float",
                      dest="threshold_matches",
                      default=0.90,
                      help="""The threshold for matches above which the contigs which align are taking into consideration. Default is '%default'.""")

    parser.add_option("--mismatches",
                      action="store",
                      type="float",
                      dest="mismatches",
                      default=1000,
                      help="""All alignments having strictly more mismatches will be removed. Default is '%default'.""")

    parser.add_option("--threshold_overlap",
                      action="store",
                      type="float",
                      dest="threshold_overlap",
                      default=17,
                      help="""The threshold for the minimum length of the read overlap over the fusion point (i.e. overhang/anchor). Default is '%default'.""")


    parser.add_option("--separator",
                      action="store",
                      type="string",
                      dest="separator",
                      default="*",
                      help="""The separator string to be used for marking the breakpoint in the fusion junction.""")


    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""A text file containg a report regarding new candidate fusion genes.""")


    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_mappings_filename and
            options.input_hugo_filename and
            options.input_genes_positions_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    print "Reading...",options.input_hugo_filename
    # it is assumed that it is in this format
    #    ensembl_gene_id
    #    hgnc_symbol
    #    description
    #    description
    hugo=dict()
    for el in [line.rstrip('\r\n').split('\t')[0:2] for line in file(options.input_hugo_filename,'r').readlines() if line.rstrip('\r\n')]:
        g=el[0]
        h=el[1].replace(' ','_')
        if hugo.has_key(g):
            hugo[g]=hugo[g]+','+h
        else:
            hugo[g]=h

    print "Reading...",options.input_genes_positions_filename
    #file: genes_positions_ensembl.txt
    #
    # columns:
    #
    # ensembl_gene_id
    # end_position
    # start_position
    # strand
    # chromosome_name
    col=dict()
    col['gene']=0
    col['end']=1
    col['start']=2
    col['strand']=3
    col['chr']=4
    gene = [line.rstrip('\r\n').split('\t') for line in file(options.input_genes_positions_filename,'r').readlines()]
    gene_dict = dict()
    for line in gene:
        g = line[col['gene']]
        s = line[col['start']]
        e = line[col['end']]
        c = line[col['chr']]
        t = line[col['strand']]
        if not gene_dict.has_key(g):
            gene_dict[g]={'chr':'','strand':'','start':0,'end':0,'hugo':''}
        else:
            if gene_dict[g]['chr'] == c and gene_dict[g]['start'] == int(s) and gene_dict[g]['end'] == int(e) and gene_dict[g]['strand'] == int(t):
                print >>sys.stderr,"WARNING: gene id %s is not unique!" % (g,)
            else:
                print >>sys.stderr,"ERROR: gene id %s is not unique!" % (g,)
                sys.exit(1)
        gene_dict[g]['chr'] = c
        gene_dict[g]['start'] = int(s)
        gene_dict[g]['end'] = int(e)
        gene_dict[g]['strand'] = int(t)
        if gene_dict[g]['start'] > gene_dict[g]['end']:
            print "Error: bad gene coordinates!"
            print g,gene_dict[g]['start'],gene_dict[g]['end'],gene_dict[g]['strand']
            sys.exit(1)
        gene_dict[g]['hugo'] = hugo[g]
    gene = gene_dict
    del gene_dict
    del hugo

    ###########################################################
    # find the skipped exon parts and exons in transcriptome
    ###########################################################
    print "Reading mappings...",options.input_mappings_filename
    data_mappings = [line.rstrip('\r\n').split('\t') for line in file(options.input_mappings_filename,'r').readlines()]

    data_mappings = [line for line in data_mappings if ((float(line[0])/float(line[10])) >= options.threshold_matches and  # number of matches / length of query sequence
                                                         line[17] == '2' and  # blockCount
                                                         line[6] == '1' and   # number inserts in target
                                                         (line[4] == '0' or (line[4] == '1' and int(line[5]) <= 3))) and  # number inserts in query
                                                         min(map(int,line[18][:-1].split(','))) >= options.threshold_overlap] # blockSizes
# original
#    data_mappings = [line for line in data_mappings if ((float(line[0])/float(line[10])) >= options.threshold_matches and  # number of matches / length of query sequence
#                                                         line[17] == '2' and  # blockCount
#                                                         line[6] == '1' and   # number inserts in target
#                                                         (line[4] == '0' or (line[4] == '1' and int(line[5]) <= 3))) and  # number inserts in query
#                                                         min(map(int,line[18][:-1].split(','))) >= options.threshold_overlap] # blockSizes


    data = [] # remove reads which do not map on two genes
    for line in data_mappings:
        # fusion genes names and border
        border = int(line[13].split('|')[2])
        # strand
        st = line[8]
        #qSize = int(line[10])
        # starts
        #p = map(int,line[20][:-1].split(',')) # start positions on target seq = tStarts
        p = map(int,line[20][:-1].split(',')) # start positions on target seq = tStarts
        #if st == '-':
        #    blockSizes = map(int,line[20][:-1].split(','))
        #    p = [qSize - blockSizes[i] - e + 1 for i,e in enumerate(p)]
        #    p.reverse()
        tStart = int(line[15])
        tEnd = int(line[16])
        if border < p[1] + 1 and border > p[0] + 1 and tStart < border and border < tEnd:
            data.append(line)

    seq_dict = None
    if options.input_genegene_fasta_filename:
        print "Reading all gene-gene's sequences...",options.input_genegene_fasta_filename
        input_seq_iterator = Bio.SeqIO.parse(open(options.input_genegene_fasta_filename, "rU"), "fasta")
        seq_dict = dict((record.id,str(record.seq)) for record in input_seq_iterator)


    flank = 50

    print "Processing..."
    result = []
    for line in data:

        # skip it if there are too many mismatches
        if int(line[1]) > options.mismatches:
            continue

        short_read = line[9]

        temp = line[13].split('|') # target sequence name
        gene_1 = temp[0]
        gene_2 = temp[1]
        border = int(temp[2])
        tsn = line[13]

        strand = line[8] # strand

        s1=0
        s2=0
        e1=0
        e2=0

        #if strand == "+":
        temp = map(int,line[20][:-1].split(',')) # start positions on target seq = tStarts
        start_1 = temp[0] + 1
        start_2 = temp[1] - border + 1
        s1 = temp[0] + 1
        s2 = temp[1] + 1

        temp = map(int,line[18][:-1].split(',')) # blockSizes
        end_1 = start_1 + temp[0] - 1
        end_2 = start_2 + temp[1] - 1
        e1 = temp[0] + s1 - 1
        e2 = temp[1] + s2 - 1

        #elif strand == "-":
        #    ts = int(line[10])
#
        #    temp = map(int,line[18][:-1].split(',')) # blockSizes
        #    bs = temp[:]
        #    temp.reverse()
        #    end_1 = start_1 + temp[0] - 1
        #    end_2 = start_2 + temp[1] - 1
        #    e1 = temp[0] + s1 - 1
        #    e2 = temp[1] + s2 - 1
#
        #    temp = map(int,line[20][:-1].split(',')) # start positions on target seq = tStarts
        #    temp = [ts - bs[i] - e + 1 for i,e in enumerate(temp)]
        #    temp.reverse()
        #    start_1 = temp[0] + 1
        #    start_2 = temp[1] - border + 1
        #    s1 = temp[0] + 1
        #    s2 = temp[1] + 1






        chr_1 = gene[gene_1]['chr']
        chr_2 = gene[gene_2]['chr']

        str_1 = gene[gene_1]['strand']
        str_2 = gene[gene_2]['strand']

        coord_gene_1 = coord_gene2genome(gene_1,[start_1,end_1],gene)
        coord_gene_2 = coord_gene2genome(gene_2,[start_2,end_2],gene)

        anchor_length = e1-s1+1 if e1-s1 < e2-s2 else e2-s2+1

        fs = ''
        if seq_dict:
            myseq = seq_dict.get(tsn,None)
            if myseq:
                e1flank = e1 - flank
                if e1flank < 0:
                    e1flank = 0
                fs = '%s%s%s' % (myseq[e1flank:e1],options.separator,myseq[s2-1:s2-1+flank])

        t = [gene_1,                 # 0
             gene[gene_1]['hugo'],   # 1
             chr_1,                  # 2
             str_1,                  # 3
             coord_gene_1[0],        # 4
             coord_gene_1[1],        # 5
             gene_2,                 # 6
             gene[gene_2]['hugo'],   # 7
             chr_2,                  # 8
             str_2,                  # 9
             coord_gene_2[0],        # 10
             coord_gene_2[1],        # 11
             short_read,             # 12
             int(line[1]),           # 13 - mismatches
             int(line[10]),          # 14
             start_1,                # 15
             end_1,                  # 16
             start_2,                # 17
             end_2,                  # 18
             anchor_length,          # 19
             fs                      # 20
             ]
        result.append(t)





    # sorting
    result = sorted(result,
                    key = lambda x: (x[0], x[6], x[5], x[10])
                    )
    header = ['gene-5end',
              'gene-5end_symbol',
              'chromosome_gene-5end',
              'strand_gene-5end',
              'start_chromosome_part-1-of-read-mapped-gene-5end',
              'end_chromosome_part-1-read-mapped-gene-5end',
              'gene-3end',
              'gene-3end_symbol',
              'chromosome_gene-3end',
              'strand_gene-3end',
              'start_chromosome_part-2-of-read-mapped-gene-3end',
              'end_chromosome_part-2-read-mapped-gene-3end',
              'short_read',
              'mismatches',
              'length_short_read',
              'start_part-1-read_on_gene-5end',
              'end_part-1-read_on_gene-5end',
              'start_part-2-read_on_gene-3end',
              'end_part-2-read_on_gene-3end',
              'anchor_length',
              'fusion_sequence'
              ]

    result.insert(0, header)
    res = result

    res = ['\t'.join(map(str,line))+'\n' for line in res]
    file(options.output_filename,'w').writelines(res)





    print "The end."
