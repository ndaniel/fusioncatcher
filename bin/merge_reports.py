#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It merges the reports from BOWTIE-method with BLAT-method.



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

import os
import sys
import optparse

def myorder(a,b):
    return (a,b) if a <= b else (b,a)


def process_psl(data, psl_filename, tag, threshold = 0, tpairs = 0, anchor2 = 40):
    # reading the BLAT report
    print "Reading...",psl_filename,anchor2
    blat = [line.rstrip('\r\n').split('\t') for line in file(psl_filename,'r') if line.rstrip('\r\n')]
    blat_header = dict([(v.lower(),i) for (i,v) in enumerate(blat.pop(0))])
    # gene-5end
    # gene-5end_symbol
    # chromosome_gene-5end
    # strand_gene-5end
    # end_chromosome_part-1-read-mapped-gene-5end
    # gene-3end
    # gene-3end_symbol
    # chromosome_gene-3end
    # strand_gene-3end
    # start_chromosome_part-2-of-read-mapped-gene-3end
    # counts
    # fusion_sequence

    # filter
    if ((threshold > 0) or (anchor2>0)):
        blat = [row for row in blat if ((int(row[blat_header['counts']]) >= threshold) or (int(row[blat_header['longest_anchor']]) >= anchor2))]


    # merging => adding to BOWTIE list
    print "Processing..."
    for line in blat:
        symbol_1 = line[blat_header['gene-5end_symbol']]
        symbol_2 = line[blat_header['gene-3end_symbol']]
        gene_1 = line[blat_header['gene-5end']]
        gene_2 = line[blat_header['gene-3end']]
        exon_1 = ''
        exon_2 = ''
        pos_1 = "%s:%s:%s" % (line[blat_header['chromosome_gene-5end']],
                              line[blat_header['end_chromosome_part-1-read-mapped-gene-5end']],
                              '+' if line[blat_header['strand_gene-5end']] == '1' else '-')
        pos_2 = "%s:%s:%s" % (line[blat_header['chromosome_gene-3end']],
                              line[blat_header['start_chromosome_part-2-of-read-mapped-gene-3end']],
                              '+' if line[blat_header['strand_gene-3end']] == '1' else '-')
        pairs =  candidate_fusions[myorder(gene_1,gene_2)]
        reads =  line[blat_header['counts']]
        longest = line[blat_header['longest_anchor']]
        aligner = tag
        seq = line[blat_header['fusion_sequence']]
        if int(pairs) >= tpairs:
            data.append([
                symbol_1,
                symbol_2,
                gene_1,
                gene_2,
                exon_1,
                exon_2,
                pos_1,
                pos_2,
                pairs,
                reads,
                longest,
                aligner,
                seq
                        ])



#####################################
#####################################
#####################################
if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It merges the reports from BOWTIE-method with BLAT-method.
"""
    version="%prog 0.10 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_bowtie",
                      action="store",
                      type="string",
                      dest="input_bowtie_filename",
                      help="""The report with candidate fusion genes found using the Bowtie.""")

    parser.add_option("--input_blat",
                      action="store",
                      type="string",
                      dest="input_blat_filename",
                      help="""The report with candidate fusion genes found using the BLAT aligner.""")

    parser.add_option("--input_star",
                      action="store",
                      type="string",
                      dest="input_star_filename",
                      help="""The report with candidate fusion genes found using the STAR aligner.""")

    parser.add_option("--input_bowtie2",
                      action="store",
                      type="string",
                      dest="input_bowtie2_filename",
                      help="""The report with candidate fusion genes found using the BOWTIE2 aligner.""")

    parser.add_option("--input_bwa",
                      action="store",
                      type="string",
                      dest="input_bwa_filename",
                      help="""The report with candidate fusion genes found using the BWA aligner.""")

    parser.add_option("--input_spotlight",
                      action="store",
                      type="string",
                      dest="input_spotlight_filename",
                      help="""The report with candidate fusion genes found using the SPOTLIGHT method.""")


    parser.add_option("--input_candidate_fusion_genes",
                      action = "store",
                      type = "string",
                      dest = "input_candidate_fusion_genes_filename",
                      help = """The input list of candidate fusion genes, for example 'candidate_fusion-genes_no-offending-reads_label-no-proteins-paralogs-readthrough-similar-pseudogenes_further.txt'.""")

    parser.add_option("--input_ambiguous",
                      action = "store",
                      type = "string",
                      dest = "input_ambiguous_filename",
                      help = """File containing the pairs of genes and their corresponding number of reads which map ambiguously on each other.""")

    parser.add_option("--supporting_reads_blat",
                      action="store",
                      type="int",
                      dest="supporting_reads_blat",
                      default=2,
                      help="""The minimum number of supporting reads (found using BLAT aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_pairs_blat",
                      action="store",
                      type="int",
                      dest="supporting_pairs_blat",
                      default=2,
                      help="""The minimum number of supporting pairs (found using BLAT aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")


    parser.add_option("--supporting_reads_star",
                      action="store",
                      type="int",
                      dest="supporting_reads_star",
                      default=2,
                      help="""The minimum number of supporting reads (found using STAR aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_pairs_star",
                      action="store",
                      type="int",
                      dest="supporting_pairs_star",
                      default=2,
                      help="""The minimum number of supporting pairs (found using STAR aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")


    parser.add_option("--supporting_reads_bowtie2",
                      action="store",
                      type="int",
                      dest="supporting_reads_bowtie2",
                      default=2,
                      help="""The minimum number of supporting reads (found using BOWTIE2 aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_pairs_bowtie2",
                      action="store",
                      type="int",
                      dest="supporting_pairs_bowtie2",
                      default=2,
                      help="""The minimum number of supporting pairs (found using BOWTIE2 aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_reads_bwa",
                      action="store",
                      type="int",
                      dest="supporting_reads_bwa",
                      default=2,
                      help="""The minimum number of supporting reads (found using BOWTIE2 aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_pairs_bwa",
                      action="store",
                      type="int",
                      dest="supporting_pairs_bwa",
                      default=2,
                      help="""The minimum number of supporting pairs (found using BWA aligner) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_reads_spotlight",
                      action="store",
                      type="int",
                      dest="supporting_reads_spotlight",
                      default=2,
                      help="""The minimum number of supporting reads (found using SPOTLIGHT method) necessary for considering valid a candidate fusion gene. Default is '%default'.""")

    parser.add_option("--supporting_pairs_spotlight",
                      action="store",
                      type="int",
                      dest="supporting_pairs_spotlight",
                      default=2,
                      help="""The minimum number of supporting pairs (found using SPOTLIGHT method) necessary for considering valid a candidate fusion gene. Default is '%default'.""")



    parser.add_option("--squish-report",
                      action = "store_true",
                      dest = "squish",
                      default = False,
                      help = """If set then the report is squished (i.e. fusion genes with same junction coordinates are listed once even that they are found by severeal methods). Default is '%default'.""")



    parser.add_option("--anchor2",
                      action = "store",
                      type = "int",
                      dest = "anchor2",
                      default = 40,
                      help = """For anchors longer (or equal) with this value it is enough to have only one supporting read. Default is '%default'.""")


    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""Merged report of candidate fusion genes.""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_bowtie_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    print "Reading...",options.input_candidate_fusion_genes_filename
    # 0  - Fusion_gene_1
    # 1  - Fusion_gene_2
    # 2  - Count_paired-end_reads
    # 3  - Fusion_gene_symbol_1
    # 4  - Fusion_gene_symbol_2
    # 5  - Information_fusion_genes
    # 6  - Analysis_status -> further or skipped
    draft_fusions = [line.rstrip('\r\n').split('\t') for line in file(options.input_candidate_fusion_genes_filename,'r').readlines()]
    draft_fusions.pop(0) # remove the header
    candidate_fusions = dict([( tuple(myorder(line[0], line[1])), line[2]) for line in draft_fusions if line[6].lower() == 'further_analysis'])
    label_fusions = dict([( tuple(myorder(line[0], line[1])), line[5]) for line in draft_fusions if line[6].lower() == 'further_analysis'])


    # reading the BOWTIE report
    # this is considered the standard
    print "Reading...",options.input_bowtie_filename
    bowtie = [line.rstrip('\r\n').split('\t') for line in file(options.input_bowtie_filename,'r') if line.rstrip('\r\n')]
    bowtie_header = dict([(v.lower(),i) for (i,v) in enumerate(bowtie[0])])
    #    Fusion_gene_symbol_1(5end_partner)              # 0
    #    Fusion_gene_symbol_2(3end_partner)              # 1
    #    Fusion_gene_1(5end_partner)                     # 2
    #    Fusion_gene_2(3end_partner)                     # 3
    #    Fusion_exon_1(5end_partner)                     # 4
    #    Fusion_exon_2(3end_partner)                     # 5
    #    Fusion_gene_1_position(5end_partner)            # 6
    #    Fusion_gene_2_position(3end_partner)            # 7
    #    Spanning_pairs                                  # 8
    #    Spanning_unique_reads                           # 9
    #    Longest_spanning_read                           # 10
    #    Aligner(s)                                      # 11
    #    Fusion_sequence                                 # 12


    if options.input_blat_filename:
        process_psl(bowtie, options.input_blat_filename, tag="BOWTIE+BLAT", threshold = options.supporting_reads_blat, tpairs = options.supporting_pairs_blat, anchor2=options.anchor2)

    if options.input_star_filename:
        process_psl(bowtie, options.input_star_filename, tag="BOWTIE+STAR", threshold = options.supporting_reads_star, tpairs = options.supporting_pairs_star, anchor2=options.anchor2)

    if options.input_bowtie2_filename:
        process_psl(bowtie, options.input_bowtie2_filename, tag="BOWTIE+BOWTIE2", threshold = options.supporting_reads_bowtie2, tpairs = options.supporting_pairs_bowtie2, anchor2=options.anchor2)

    if options.input_bwa_filename:
        process_psl(bowtie, options.input_bwa_filename, tag="BOWTIE+BWA", threshold = options.supporting_reads_bwa, tpairs = options.supporting_pairs_bwa, anchor2=options.anchor2)

    if options.input_spotlight_filename:
        process_psl(bowtie, options.input_spotlight_filename, tag="BOWTIE+SPOTLIGHT", threshold = options.supporting_reads_spotlight, tpairs = options.supporting_pairs_spotlight, anchor2=options.anchor2)



    # add an extra column to the report with the labels of the fusion genes
    new_bowtie = []
    for i,line in enumerate(bowtie):
        if i == 0:
            new_bowtie.append(line + ['Fusion_description'])
        else:
            gene_1 = line[2]
            gene_2 = line[3]
            new_bowtie.append(line + [label_fusions[myorder(gene_1,gene_2)]])
    bowtie = new_bowtie

    # add an extra column to the report with the counts of ambiguous counts of reads
    ambiguous = dict()
    if options.input_ambiguous_filename and (os.path.isfile(options.input_ambiguous_filename) or os.path.islink(options.input_ambiguous_filename)):
        ambiguous = [line.rstrip('\r\n').split('\t') for line in file(options.input_ambiguous_filename,'r').readlines() if line.rstrip('\r\n')]
        ambiguous = dict([(myorder(line[0],line[1]),line[2]) for line in ambiguous if label_fusions.has_key(myorder(line[0],line[1])) ])
    new_bowtie = []
    for i,line in enumerate(bowtie):
        if i == 0:
            new_bowtie.append(line + ['Counts_of_common_mapping_reads'])
        else:
            gene_1 = line[2]
            gene_2 = line[3]
            new_bowtie.append(line + [ambiguous.get(myorder(gene_1,gene_2),'0')])
    bowtie = new_bowtie

    # here reshuffle the columns so that they look beautiful
    # move column Fusion_description from 13 to 3
    # move column Counts_of_common_mapping_reads from 14 to 4 (after Fusion_description)
    # Spanning_pairs from column 8 to 5
    # Spanning_unique_reads from column 9 to 6
    # Longest_anchor_found from column 10 to 7
    # Fusion_finding_method from column 11 to 8
    new_bowtie = zip(*bowtie)
    new_bowtie = [ new_bowtie[0],
                   new_bowtie[1],
                   new_bowtie[13],
                   new_bowtie[14],
                   new_bowtie[8],
                   new_bowtie[9],
                   new_bowtie[10],
                   new_bowtie[11],
                   new_bowtie[6],
                   new_bowtie[7],
                   new_bowtie[2],
                   new_bowtie[3],
                   new_bowtie[4],
                   new_bowtie[5],
                   new_bowtie[12]
    ]
    new_bowtie = zip(*new_bowtie)
    header = new_bowtie.pop(0)
    bowtie = sorted(new_bowtie, key = lambda x: (-int(x[4]),x[0],x[1],-int(x[5]),-int(x[6]),x[7]))
    bowtie.insert(0,header)

    # squish report
    if options.squish and bowtie and len(bowtie)>1:
        b = bowtie[:]
        h = b.pop(0) # header
        r = [h]
        n = len(b)
        f = [True] * n
        clean_labels = {
            frozenset(['distance200kbp', 'distance100kbp', 'distance10kbp', 'distance1000bp']): 'gap<1K',
            frozenset(['distance200kbp', 'distance100kbp', 'distance10kbp']): '1K<gap<10K',
            frozenset(['distance200kbp', 'distance100kbp']): '10K<gap<100K',
            frozenset(['distance200kbp']): '100K<gap<200K',
            }
        for i in xrange(n):
            if f[i]:
                method = set([b[i][7]])
                exon1 = ''
                exon2 = ''
                ids = []
                for j in xrange(i+1,n):
                    if f[j] and b[i][8] == b[j][8] and b[i][9] == b[j][9] and b[i][10] == b[j][10] and b[i][11] == b[j][11]:
                        method.add(b[j][7])
                        if b[j][12]:
                            exon1 = b[j][12]
                            exon2 = b[j][13]
                        f[j] = False
                line = list(b[i][:])
                if method:
                    line[7] = ';'.join(sorted(method))
                if exon1 and (not line[12]):
                    line[12] = exon1
                    line[13] = exon2
                # convert IGH_LOCUS to IGH@ and IGK_LOCUS to IGK@
                if line[0].lower().find('_locus') != -1:
                    line[0] = line[0].split('_')[0].upper()+'@'
                if line[1].lower().find('_locus') != -1:
                    line[1] = line[1].split('_')[0].upper()+'@'
                # clean the short_distance and distance labels
                labels = line[2].split(',')
                labels = [elk for elk in labels if elk != 'short_distance'] # remove short_distance
                distances = frozenset([elk for elk in labels if elk.startswith('distance')])
                if distances:
                    labels = [elk for elk in labels if not elk.startswith('distance')] # remove the distance labels
                    labels.append(clean_labels.get(distances,''))
                    line[2] = ','.join(labels)
                r.append(line)
        bowtie = r


    # FINAL REPORT
    print "Writing the merged report..."
    file(options.output_filename,'w').writelines(['\t'.join(line)+'\n' for line in bowtie])
    #
