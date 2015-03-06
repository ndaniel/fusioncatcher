#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It finds a list of genes that might homologous (there is a short read which maps on both genes).



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

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

#########################
def line_from(a_map_filename):
    # it gives chunks from a_map_filename which is assumed to be ordered by the name of transcripts (i.e. column 3)
    fin = None
    if a_map_filename == '-':
        fin = sys.stdin
    elif a_map_filename.lower().endswith('.gz'):
        fin = gzip.open(a_map_filename,'r')
    else:
        fin = open(a_map_filename,'r')

    while True:
        gc.disable()
        lines=fin.readlines(10**8)
        gc.enable()
        if not lines:
            break
        gc.disable()
        lines=[line.rstrip('\r\n').split('\t',3)[:3] for line in lines if line.rstrip('\r\n')]
        gc.enable()
        for line in lines:
            yield line
    fin.close()

#########################
def read_from(a_map_filename):
    last_r=''
    chunk=set()
    last_g = '' # for speed purposes only
    for line in line_from(a_map_filename):
        if not chunk:
            last_r=line[0]
        if last_r!=line[0]: # line[2] is column no 3 in the BOWTIE MAP file which contains the reference sequence name
            yield (last_r,chunk)
            last_r=line[0]
            chunk=set()
            last_g = ''
        #tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103
        #g=[el.split('ge=')[1] for el in line[2].split(';') if el.startswith('ge=')]
        #chunk.add(g[0])
#        g1 = line[2].find('ge=') + 3
#        g2 = line[2].find(';',g1+1)
#        g = line[2][g1:g2]
        g = line[2].partition(';')[2]
        if g != last_g:
            chunk.add(g)
            last_g = g
    if chunk:
        yield (last_r,chunk)

#########################
def is_overlapping(g1,g2,tolerance = 50):
    r = False
    e1 = database[g1]
    e2 = database[g2]
    if e1[2] == e2[2]: # same chromosomes
        if e1[1] - e2[0] > tolerance:
            if e2[1] - e1[0] > tolerance:
                r = True
    return r


if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It finds a list of genes that might homologous (there is a short read which maps on both genes)."""
    version="%prog 0.10 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_map_filename",
                      help="""The input file in Bowtie MAP format (sorted by read name) containing the short reads mapped on the transcripts (can be also STDIN).""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text tab-separated file containing the candidate homologous genes (the genes are sorted alphabetically on the each line).""")

    parser.add_option("--reads",
                      action="store",
                      type="int",
                      dest="reads",
                      default = 1,
                      help="""The minimum number of reads which map simultaneously on two genes in order to be considered as homolog genes. Default is %default.""")


    parser.add_option("--output_offending_reads",
                      action="store",
                      type="string",
                      dest="output_offending_reads_filename",
                      help="""The output text file containing the reads names which mapp simultaneously on transcripts from at least two genes.""")

    parser.add_option("--output_offending_pair_reads",
                      action="store",
                      type="string",
                      dest="output_offending_pair_reads_filename",
                      help="""The output text file containing the reads names (and its mate) which mapp simultaneously on transcripts from at least two genes.""")

    parser.add_option("--input_exons",
                      action="store",
                      type="string",
                      dest="exons_filename",
                      help="""Database with exons position on chromosomes, e.g. 'more_exons_ensembl.txt'. This is used for filtering the UTRs extensions by removing any extension which overlaps with any exons from the database. This is optional.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_map_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    database = {}
    if options.exons_filename:
        #print "Processing the exons database..."
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
                  line[12] # chromosome           3
                  ) for line in exons]

        for line in exons:
            gn = line[0]
            gs = int(line[1])
            ge = int(line[2])
            ch = line[3].upper()

            (gs,ge) = (ge,gs) if gs>ge else (gs,ge)
            if not database.has_key(gn):
                database[gn] = (gs,ge,ch)

    #print "Finding the homolog genes..."
    homolog=dict()
    offenders=set()
    for (a_read,genes) in read_from(options.input_map_filename):
        n=len(genes)
        if n==1:
            continue
        g=list(genes)
        flag = False
        for a in xrange(0,n-1):
            for b in xrange(a+1,n):
                k = ''
                if g[a] < g[b]:
                    k = '%s\t%s' % (g[a],g[b])
                else:
                    k = '%s\t%s' % (g[b],g[a])
                gc.disable()
                homolog[k] = homolog.get(k,0) + 1
                gc.enable()
                if database and not is_overlapping(g[a],g[b]):
                    flag = True
        if (not database) or flag:
            gc.disable()
            offenders.add(a_read)
            gc.enable()

    #print "Writing...",options.output_filename
    #homolog = sorted([k+'\t'+str(v)+'\n' for (k,v) in homolog.items() if v >= options.reads])
    homolog = [k+'\t'+str(v)+'\n' for (k,v) in homolog.items() if v >= options.reads]
    file(options.output_filename,'w').writelines(homolog)

    if options.output_offending_reads_filename:
        #print "Writing...",options.output_offending_reads_filename
        d = sorted(list(offenders))
        file(options.output_offending_reads_filename,'w').writelines([line+'\n' for line in d])

    if options.output_offending_pair_reads_filename:
        #print "Writing...",options.output_offending_pair_reads_filename
        d = set()
        for e in offenders:
            if e.endswith('/1'):
                d.add(e[:-1]+'2')
            elif e.endswith('/2'):
                d.add(e[:-1]+'1')
        gc.disable()
        offenders.update(d)
        gc.enable()
        d = sorted(list(offenders))
        file(options.output_offending_pair_reads_filename,'w').writelines([line+'\n' for line in d])

    #print "The end."
    #
