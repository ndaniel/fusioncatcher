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


def shake(aline,nc1,nc2):
    # re-aligns the read using the new cuts (nc1,nc2)
    r = aline[:]
    fs = r[col_fs].split("*")
    star = "*"
    fs1 = fs[0]
    fs2 = fs[1]
    if fs2.startswith("N"):
        p = max([fs2.find(ek) for ek in ("NA","NC","NG","NT")])
        star = star + fs2[:p+1]
        fs2 = fs2[p+1:]

    if r[col_g5e] != nc1:
        w = abs(int(r[col_g5e]) - int(nc1))
        r[col_g5e] = nc1
        fs1 = fs1[:-w]
    if r[col_g3s] != nc2:
        w = abs(int(r[col_g3s]) - int(nc2))
        r[col_g3s] = nc2
        fs2 = fs1[w:]
    r[col_fs] = "%s%s%s" % (fs1,star,fs2)
    # debug
    #r[col_g5] = "->"+r[col_g5]
    return r

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It re-arranges the split reads which support the fusion junctions such that they more favourable for finding the fusions."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file containing the reads alignments supporting potential fusion genes.""")


    parser.add_option("--wiggle","-w",
                      action = "store",
                      type = "int",
                      dest = "wiggle",
                      default = 3,
                      help = "The maximum size of the clipping done in the alignment gaps of the reads. "+
                             "Default is '%default'. ")

    parser.add_option("--output",'-o',
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      default = None,
                      help = "The output file containing reads alignments suporting potential fusion genes, such that the split reads are overlapping more.")




    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
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

    data = [line.rstrip('\r\n').split('\t') for line in file(options.input_filename,'r') if line.rstrip('\r\n')]
    header = data.pop(0)
    col = {}
    if header:
        col = dict([(e.lower(),i) for i,e in enumerate(header)])
        col_read = col['short_read']
        col_g5 = col['gene-5end']
        col_sg5 = col['strand_gene-5end']
        col_cg5 = col['chromosome_gene-5end']
        col_g3 = col['gene-3end']
        col_sg3 = col['strand_gene-3end']
        col_cg3 = col['chromosome_gene-3end']
        col_g5s = col['start_chromosome_part-1-of-read-mapped-gene-5end']
        col_g5e = col['end_chromosome_part-1-read-mapped-gene-5end']
        col_g3s = col['start_chromosome_part-2-of-read-mapped-gene-3end']
        col_g3e = col['end_chromosome_part-2-read-mapped-gene-3end']
        col_fs = col['fusion_sequence']


    wiggle = options.wiggle


    result = []
    
    #
    if data:

        genes53 = set([(line[col_g5],line[col_sg5],line[col_cg5],line[col_g3],line[col_sg3],line[col_cg3]) for line in data])

        reads = dict()
        for line in data:
            c = (line[col_g5e],line[col_g3s],line[col_read])
            reads[c] = line


        for (g5,s5,c5,g3,s3,c3) in genes53:
        
            d = [line for line in data if line[col_g5] == g5 and line[col_sg5] == s5 and line[col_cg5] == c5 and line[col_g3] == g3 and line[col_sg3] == s3 and line[col_cg3] == c3]
            
            # find the cuts
            cuts = dict()

            for line in d:
                c = (line[col_g5e],line[col_g3s])
                if not cuts.has_key(c):
                    cuts[c] = dict()
                cuts[c][line[col_read]] = line


            
            
                
            n = len(cuts)
            k = cuts.keys()
            for i in xrange(n-1):
                for j in xrange(i+1,n):
                    cut_a = cuts[k[i]]
                    cut_b = cuts[k[j]]
                    cut_a1 = int(k[i][0])
                    cut_a2 = int(k[i][1])
                    cut_b1 = int(k[j][0])
                    cut_b2 = int(k[j][1])
                    if abs(cut_a1 - cut_b1) <= wiggle and abs(cut_a2 - cut_b2) <= wiggle:
                        cut_new1 = min(cut_a1,cut_b1)
                        if s5 == "-1":
                            cut_new1 = max(cut_a1,cut_b1)
                        cut_new2 = max(cut_a2,cut_b2)
                        if s3 == "-1":
                            cut_new2 = min(cut_a2,cut_b2)
                        cut_new1 = str(cut_new1)
                        cut_new2 = str(cut_new2)
                        #print cut_a1,cut_a2
                        #print cut_b1,cut_b2
                        #print cut_new1,cut_new2
                        #print "--------"
                        # proceed
                        for cc in cut_a.keys():
                            cu = (cut_new1,cut_new2,cc)
                            if cu not in reads:
                                reads[cu] = shake(cut_a[cc],cut_new1,cut_new2)
                        for cc in cut_b.keys():
                            cu = (cut_new1,cut_new2,cc)
                            if cu not in reads:
                                reads[cu] = shake(cut_b[cc],cut_new1,cut_new2)
        # add the results from reads
        result = reads.values()
    result.insert(0,header)
    file(options.output_filename,'w').writelines(['\t'.join(line)+'\n' for line in result])
