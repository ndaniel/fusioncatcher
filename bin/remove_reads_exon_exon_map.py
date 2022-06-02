#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given a input MAP file containing reads mapping on exon-exon junctions it will
remove the reads which have the mate mapping (i.e. second MAP input file
containing reads mapping on transcriptom) on other genes than those
containing the exons which form the exon-exon junction.



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
remove_reads_exon_exon.py \
--input_exon_exon     reads_mapped-exon-exon-fusion-genes_sorted-ref.map \
--input_transcriptome reads_filtered_transcriptome_sorted-read.map \
--output_exon_exon    reads_mapped-exon-exon-fusion-genes_sorted-ref_filtered.map
"""
import sys
import os
import optparse
import gzip

#########################
def line_from(a_map_filename):
    # it gives chunks from a_map_filename which is assumed to be ordered by the name of transcripts (i.e. column 3)
    # col 1 => read name
    # col 3 => name sequence on which read is aligning
    fin = None
    if a_map_filename.lower().endswith('.gz'):
        fin = gzip.open(a_map_filename,'r')
    else:
        fin = open(a_map_filename,'r')
    while True:
        lines=fin.readlines(10**8)
        if not lines:
            break
        lines = [line.rstrip('\r\n').split('\t')[:3] for line in lines if line.rstrip('\r\n')]
        for line in lines:
            yield line
    fin.close()

#########################
def read_from(a_map_filename, reads):
    # process only the reads which are in the set reads!
    last_r = ''
    chunk = set()
    last_g = '' # for speed purposes only
    skip = False # output only the reads which appear in 'reads'
    for line in line_from(a_map_filename):
        #if not chunk:
        #    last_r = line[0]
        #    if last_r in reads:
        #        skip = False
        #    else:
        #        skip = True
        if last_r != line[0]: # line[2] is column no 3 in the BOWTIE MAP file which contains the reference sequence name
            if chunk and not skip:
                yield (last_r,chunk)
            last_r = line[0]
            if last_r in reads:
                skip = False
            else:
                skip = True
            chunk = set()
            last_g = ''
        if not skip:
            #tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103
            #g=[el.split('ge=')[1] for el in line[2].split(';') if el.startswith('ge=')]
            #chunk.add(g[0])
#            g1 = line[2].find('ge=') + 3
#            g2 = line[2].find(';',g1+1)
#            g = line[2][g1:g2]
            g = line[2].partition(';')[2]
            if g != last_g:
                chunk.add(g)
                last_g = g
    if last_r and chunk and not skip:
        yield (last_r,chunk)


#
#
#
if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Given a input MAP file containing reads mapping on
                     exon-exon junctions it will remove the reads which have
                     the mate mapping (i.e. second MAP input file
                     containing reads mapping on transcriptom) on other genes
                     than those containing the exons which form the exon-exon
                     junction."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input_exon_exon",
                      action = "store",
                      type = "string",
                      dest = "input_exon_exon_filename",
                      help = """The input MAP file containing the reads mapping
                           on exon-exon junctions.""")

    parser.add_option("--input_transcriptome",
                      action = "store",
                      type = "string",
                      dest = "input_transcriptome_filename",
                      help = """The input MAP file containing the reads mapping
                             on transcriptome.""")

    parser.add_option("--only_pairs",
                      action = "store_true",
                      default = False,
                      dest = "only_pairs",
                      help = """Only reads which form a pair are kept for
                      further analysis (i.e. one read maps on one of the known
                      transcripts of the genes involved in the candidate
                      fusion and its corresponding mate maps on fusion point
                      which is the exon-exon junction). Default is %default.""")

    parser.add_option("--output_exon_exon",
                      action = "store",
                      type = "string",
                      dest = "output_exon_exon_filename",
                      help = """The output text file containing all reads mapping
                             on the exon-exon junctions from the input MAP file
                             except the ones which have been removed (because
                             their mates map on other genes than those from which
                             the exons which form the exon-exon junction).""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_exon_exon_filename and
            options.input_transcriptome_filename and
            options.output_exon_exon_filename
            ):
        parser.print_help()
        sys.exit(1)

    # read the exon-exon junctions
    # get: read name, gene 5 name, gene 3 name
    # col 1: F1000050402361/1
    # col 3: ENSG00000002586-ENSG00000210082;ENST00000381180-ENST00000387347;ENSE00001669439-ENSE00001544497;rank_exon=2-1;length=80;junction=40;index_seq=8
    print "Reading exon-exon junctions mappings...",options.input_exon_exon_filename
    exjuncs = [line.rstrip('\r\n').split('\t',3)[:3] for line in file(options.input_exon_exon_filename,'r').readlines() if line.rstrip('\r\n')]
    ex = dict()
    for line in exjuncs:
        r = line[0] # in exon-exon junction
        m = r[:-1]+'2' if r.endswith('/1') else r[:-1]+'1' # in transcriptome
        g = line[2].split(';',1)[0]
        g2 = g.split('-',1)
        if not ex.has_key(m):
            ex[m] = {'mate':r, 'tr':False, 'gene':dict()}

        if not ex[m]['gene'].has_key(g2[0]):
            ex[m]['gene'][g2[0]] = dict()
        ex[m]['gene'][g2[0]][g] = {'erase':False, 'pair': g2[1]} # later flag it True if it needs to be removed


        if not ex[m]['gene'].has_key(g2[1]):
            ex[m]['gene'][g2[1]] = dict()
        ex[m]['gene'][g2[1]][g] = {'erase':False, 'pair': g2[0]} # later flag it True if it needs to be removed



    # read the transcriptome mappings
    # get: read name and gene on which it maps
    # col 1: F1000050402361/1
    # col 3: tr=ENST00000449131;ge=ENSG00000167995;pn=ENSP00000399709;chr=11;str=+;len=4267
    print "Reading transcriptome mappings...",options.input_transcriptome_filename
    i = 0
    mates = set(ex.keys())

    on_transcriptome = set()
    for (a_read, genes) in read_from(options.input_transcriptome_filename, reads = mates):
        x = ex[a_read]['gene']
        on_transcriptome.add(ex[a_read]['mate']) # its mate found to map on transcriptome
        for (k1,v1) in x.items():
            for (k2,v2) in v1.items():
                # k1 -> gene1
                # k2 -> gene1-gene2
                # v2['pair'] -> gene2
                # v2['erase'] -> True or False
                if ((not v2['erase']) and (k1 not in genes) and (v2['pair'] not in genes)):
                        # mark it for deletion
                        x[k1][k2]['erase'] = True
                        x[v2['pair']][k2]['erase'] = True
                        i = i + 1
                        #print "The read %s from exon-exon-junc (gene-gene: %s) has the mate %s which maps on transcritome on genes: %s" % (ex[a_read]['mate'],k2,a_read,','.join(genes))
        ex[a_read]['gene'] = x
    print "%d lines marked for deletion..." % (i,)

    # arrange nicely the reads which need to be erased
    erase = dict()
    for (a_mate,v1) in ex.iteritems():
        a_read = v1['mate']
        if not erase.has_key(a_read):
            erase[a_read] = set()
        for (k2,v2) in v1['gene'].iteritems():
            for (k3,v3) in v2.iteritems():
                if v3['erase']:
                    erase[a_read].add(k3)
                    #print a_read,k3


    #
    print 'Writing...',options.output_exon_exon_filename
    fi=open(options.input_exon_exon_filename,'r')
    fo=open(options.output_exon_exon_filename,'w')
    while True:
        lines = fi.readlines(10**8)
        if not lines:
            break
        lines = [line.rstrip('\r\n').split('\t') for line in lines if line.rstrip('\r\n')]
        #
        if options.only_pairs:
            lines = [line for line in lines if  (line[0] in on_transcriptome and # read; its mate is found to map on transcriptome
                                                 not (erase.has_key(line[0]) and # read
                                                 line[2].split(';',1)[0] in erase[line[0]])) # junction
                                                ]
        else:
            lines = [line for line in lines if  (not (erase.has_key(line[0]) and # read
                                                 line[2].split(';',1)[0] in erase[line[0]])) # junction
                                                ]
        #
        lines = ['\t'.join(line)+'\n' for line in lines]
        fo.writelines(lines)
    fo.close()
    fi.close()

    print "The end."
