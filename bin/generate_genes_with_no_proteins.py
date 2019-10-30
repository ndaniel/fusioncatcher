#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates list of genes which have no protein product.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2019 Daniel Nicorici

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


#####################################################################
#####################################################################
#####################################################################
# main

if __name__=='__main__':
    print "Starting..."

    #command line parsing

    usage="%prog [options]"
    description="""It generates a text file with the ensembl gene ids which do not have a protein product."""
    version="%prog 0.10 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""A text file containg all the information regarding exons, genes, proteins and their positions. (see 'more_exons_ensembl.txt' file)""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output file containing all Ensembl gene ids which have no protein product.""")

    parser.add_option("--header",
                      action="store",
                      type="string",
                      dest="output_header_filename",
                      help="""The header of the output file.""")


    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    
    print "Reading the exon information...", options.input_filename
    database=[line.strip('\r\n').split('\t') for line in file(options.input_filename,'r').readlines()]
    # the database is expected to have the following columns
    # 1: ensembl_peptide_id
    # 2: ensembl_gene_id
    # 3: ensembl_transcript_id
    # 4: ensembl_exon_id
    # 5: exon_chrom_start
    # 6: exon_chrom_end
    # 7: rank (of exon in the transcript)
    # 8: start_position (of gene)
    # 9: end_position (of gene)
    # 10: transcript_start
    # 11: transcript_end
    # 12: strand (of chromosome)
    # 13: chromosome_name
    columns=[
        'ensembl_peptide_id',
        'ensembl_gene_id',
        'ensembl_transcript_id',
        'ensembl_exon_id',
        'exon_chrom_start',
        'exon_chrom_end',
        'rank',
        'start_position',
        'end_position',
        'transcript_start',
        'transcript_end',
        'strand',
        'chromosome_name']
    col=dict([(t,i) for i,t in enumerate(columns)])
    g_col=col['ensembl_gene_id']
    p_col=col['ensembl_peptide_id']
    genes_no_proteins=set()
    genes_with_proteins=set()
    for line in database:
        if line[p_col]:
            genes_with_proteins.add(line[g_col])
        else:
            genes_no_proteins.add(line[g_col])
    genes=genes_no_proteins.difference(genes_with_proteins)
    genes=[el+'\n' for el in sorted(list(genes))]

    print "Writing the genes...",options.output_filename
    file(options.output_filename,'w').writelines(genes)

    if options.output_header_filename:
        file(options.output_header_filename,'w').write('ensembl_gene_id\n')
