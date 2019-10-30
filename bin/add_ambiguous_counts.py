#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It adds an extra column with the counts of reads which are ambiguous.



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

import os
import sys
import optparse

def myorder(a,b):
    return (a,b) if a <= b else (b,a)

#####################################
#####################################
#####################################
if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It adds an extra column with the counts of reads which are ambiguous."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """File containing the candidate fusion genes.""")


    parser.add_option("--input_ambiguous",
                      action = "store",
                      type = "string",
                      dest = "input_ambiguous_filename",
                      help = """File containing the pairs of genes and their corresponding number of reads which map ambiguously on each other.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""File containing the candidate fusion genes from the input and an extra column with their associatied counts of common mapping reads.""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    print "Reading...",options.input_filename
    # 0  - Fusion_gene_1
    # 1  - Fusion_gene_2
    # 2  - Count_paired-end_reads
    # 3  - Fusion_gene_symbol_1
    # 4  - Fusion_gene_symbol_2
    # 5  - Information_fusion_genes
    # 6  - Analysis_status -> further or skipped
    draft_fusions = [line.rstrip('\r\n').split('\t') for line in file(options.input_filename,'r').readlines()]
    header = draft_fusions.pop(0) # remove the header
    candidate_fusions = set([ myorder(line[0], line[1]) for line in draft_fusions])

    print "Reading...",options.input_ambiguous_filename
    # add an extra column to the report with the counts of ambiguous counts of reads
    ambiguous = dict()
    if options.input_ambiguous_filename and (os.path.isfile(options.input_ambiguous_filename) or os.path.islink(options.input_ambiguous_filename)):
        ambiguous = [line.rstrip('\r\n').split('\t') for line in file(options.input_ambiguous_filename,'r').readlines() if line.rstrip('\r\n')]
        ambiguous = dict([(myorder(line[0],line[1]),line[2]) for line in ambiguous if myorder(line[0],line[1]) in candidate_fusions ])
    new_draft_fusions = [ header + ['Counts_of_common_mapping_reads']]
    for line in draft_fusions:
        gene_1 = line[0]
        gene_2 = line[1]
        new_draft_fusions.append(line + [ambiguous.get(myorder(gene_1,gene_2),'0')])

    # FINAL REPORT
    print "Writing the merged report...",options.output_filename
    file(options.output_filename,'w').writelines(['\t'.join(line)+'\n' for line in new_draft_fusions])
    #
