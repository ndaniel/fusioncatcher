#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It builds a list of overlapping genes.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2016 Daniel Nicorici

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

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It builds a list of overlapping genes."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

    parser.add_option("--input_genes",
                      action="store",
                      type="string",
                      dest="input_genes_positions",
                      help="""Input file with genes positions.""")

    parser.add_option("--head",
                      action="store",
                      type="string",
                      dest="head",
                      default = "",
                      help="""This is added in the beginning of the output filenames. Default is %default.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="The output directory where the genes sequences "+
                           "are written. Default is '%default'.")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_genes_positions and
            options.output_directory
            ):
        parser.print_help()
        parser.error("Missing argument(s)!")
        sys.exit(1)


    #
    #
    #

    genes_database_filename = options.input_genes_positions
    fully_filename = os.path.join(options.output_directory,'%sfully_overlapping_genes.txt' % (options.head,))
    partially_filename = os.path.join(options.output_directory,'%spartially_overlapping_genes.txt' % (options.head,))
    same_strand_filename = os.path.join(options.output_directory,'%ssame_strand_overlapping_genes.txt' % (options.head,))

    # of type
    # for example file: genes.txt
    #
    # columns:
    #
    # ensembl_gene_id
    # end_position
    # start_position
    # strand
    # chromosome_name

    col = dict()
    col['gene'] = 0
    col['end'] = 1
    col['start'] = 2
    col['strand'] = 3
    col['chr'] = 4

    database = [line.rstrip('\r\n').split('\t') for line in file(genes_database_filename,'r').readlines() if line.rstrip('\r\n')]

    chrom = set([line[col['chr']] for line in database])

    fully = set()
    partially = set()
    same_strand = set()
    tolerance = 40

    for c in chrom:
        cr = col['chr']
        start = col['start']
        end = col['end']
        name = col['gene']
        strand = col['strand']
        data = [(line[name],int(line[start]),int(line[end]),line[strand]) for line in database if line[cr] == c]
        data = [line if line[1] <= line[2] else (line[0],line[2],line[1],line[3]) for line in data]
        data = sorted( data, key=lambda x: (x[1],x[2]))
        print 'chromsome =',c,'has', len(data),'known genes'
        n = len(data)
        for i in xrange(n-1):
            g_1 = data[i][0]
            g_1s = data[i][1]
            g_1e = data[i][2]
            g_1d = data[i][3]
            for j in xrange(i+1,n):
                g_2 = data[j][0]
                g_2s = data[j][1]
                g_2e = data[j][2]
                g_2d = data[j][3]
                if g_2 == g_1:
                    continue
                d1 = g_1e - g_2s
                d2 = g_2e - g_1s
                if d1 > tolerance and d2 > tolerance:
                    k = ""
                    if g_1 <= g_2:
                        k = "%s\t%s\n" % (g_1,g_2)
                    else:
                        k = "%s\t%s\n" % (g_2,g_1)
                    if ((g_1s - g_2s > -tolerance and g_2e - g_1e > -tolerance) or
                        (g_2s - g_1s > -tolerance and g_1e - g_2e > -tolerance)):
                        # fully overlapping
                        fully.add(k)
                        if g_1d == g_2d:
                            same_strand.add(k)
                    else:
                        partially.add(k)
                        if g_1d == g_2d:
                            same_strand.add(k)
    print "Found",len(fully),"fully overlapping genes (on same strand or on different strands)."
    print "Found",len(partially),"partially overlapping genes (on same strand or on different strands)."
    print "Found",len(same_strand),"(partially and fully) overlapping genes on the same strand."

    fully = sorted(fully)
    partially = sorted(partially)
    same_strand = sorted(same_strand)
    file(fully_filename,'w').writelines(fully)
    file(partially_filename,'w').writelines(partially)
    file(same_strand_filename,'w').writelines(same_strand)
    #
