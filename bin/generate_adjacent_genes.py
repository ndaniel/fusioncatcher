#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It builds a list of adjacent genes. Two genes are considered NOT adjacent if
there is no other gene located between the two genes on the same strand of the
chromosome. The gene in between is not allowed to overlap neither of the two
genes.



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

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It builds a list of adjacent genes. Two genes are considered NOT adjacent if
there is no other gene located between the two genes on the same strand of the
chromosome. The gene in between is not allowed to overlap neither of the two
genes."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_genes",
                      action="store",
                      type="string",
                      dest="input_genes_positions",
                      help="""Input file with genes positions.""")


    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the genes sequences are written. Default is '%default'.""")


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
    adjacent_database_filename = os.path.join(options.output_directory,'adjacent_genes.txt')
    readthrough_database_filename = os.path.join(options.output_directory,'readthroughs.txt')

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

    col=dict()
    col['gene']=0
    col['end']=1
    col['start']=2
    col['strand']=3
    col['chr']=4

    database=[line.rstrip('\r\n').split('\t') for line in file(genes_database_filename,'r').readlines()]

    chrom=set([line[col['chr']] for line in database])

    stra=set([line[col['strand']] for line in database])

    adjacent = set()
    readthrough = set()

    cr=col['chr']
    start=col['start']
    end=col['end']
    name=col['gene']
    strand=col['strand']

    for c in chrom:
        for st in stra:
            forward = False if st.startswith('-') else True
            data=[(line[0],int(line[1]),int(line[2]),int(line[3]),line[4]) for line in database if line[cr]==c and line[strand]==st]
            data=sorted( data, key=lambda x: (x[start],x[end]))
            print 'chromsome=',c,'has', len(data),'known genes on',st,'strand.'
            n=len(data)
            for i in range(n-1):
                g_1=data[i][name]
                g_1s=data[i][start]
                g_1e=data[i][end]
                for j in range(i+1,n):
                    g_2=data[j][name]
                    g_2s=data[j][start]
                    g_2e=data[j][end]
                    if g_2==g_1:
                        continue

                    flag=True
                    for k in range(i+1,j):
                        g_3=data[k][name]
                        g_3s=data[k][start]
                        g_3e=data[k][end]
                        if g_2==g_1 or g_1==g_3:
                            continue
                        if g_1e<g_3s and g_3e<g_2s:
                            flag=False
                            break
                    if flag: # they are adjacent
                        adjacent.add('\t'.join(sorted([g_1,g_2]))+'\n')
                        if forward:
                            readthrough.add("%s\t%s\n"%(g_1,g_2))
                        else:
                            readthrough.add("%s\t%s\n"%(g_2,g_1))
                    else:
                        break
    print "Found",len(adjacent)," adjacent genes."
    print "Found",len(readthrough)," reathroughs."

    file(adjacent_database_filename,'w').writelines(sorted(list(adjacent)))
    file(readthrough_database_filename,'w').writelines(sorted(list(readthrough)))
    #
