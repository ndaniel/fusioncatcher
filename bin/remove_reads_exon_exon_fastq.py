#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given a input FASTQ file containing unmapped reads it will
remove the reads which have the mate mapping (i.e. second MAP input file
containing reads mapping on transcriptom) on other genes than the list of
candidate fusion genes.



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

#
"""
remove_reads_exon_exon_psl.py \
--input_psl     reads_best_unique_blat_mapped_on_fusion_genes.psl \
--input_transcriptome     reads_filtered_transcriptome_sorted-read.map \
--output_psl    reads_best_unique_blat_mapped_on_fusion_genes_pairs.psl
"""
import sys
import os
import optparse
import gc
import gzip
import itertools

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
    buffer_size = 10 **8
    while True:
        lines=fin.readlines(buffer_size)
        if not lines:
            break
        gc.disable()
        lines = [line.rstrip('\r\n').split('\t')[:3] for line in lines if line.rstrip('\r\n')]
        gc.enable()
        for line in lines:
            yield line
    fin.close()

#########################
def read_from(a_map_filename, reads, genes):
    # process only the reads which are in the set reads!
    last_r = ''
    chunk = False
    last_g = '' # for speed purposes only
    skip = False # output only the reads which appear in 'reads' and map on genes
    for line in line_from(a_map_filename):
        if last_r != line[0]: # line[2] is column no 3 in the BOWTIE MAP file which contains the reference sequence name
            if chunk and not skip:
                yield last_r
#                print "yield",last_r
            last_r = line[0]
            if last_r in reads:
                skip = False
            else:
                skip = True
            chunk = False
            last_g = ''
        if not skip:
            g = line[2].partition(';')[2]
            if g != last_g:
                last_g = g
                if g in genes:
                    chunk = True
    if last_r and chunk and not skip:
        yield last_r


#
#
#
if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Given a input FASTQ file containing unmapped reads mapping
                     it will remove the reads which have
                     the mate mapping (i.e. second MAP input file
                     containing reads mapping on transcriptom) on other genes
                     than the list of candidate fusion genes."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input_fastq",
                      action = "store",
                      type = "string",
                      dest = "input_fastq_filename",
                      help = """The input FASTQ file containing the reads.""")

    parser.add_option("--input_fusions",
                      action = "store",
                      type = "string",
                      dest = "input_fusions_filename",
                      help = """The input file containing the list of candidate fusion genes.""")

    parser.add_option("--input_transcriptome",
                      action = "store",
                      type = "string",
                      dest = "input_transcriptome_filename",
                      help = """The input MAP file containing the reads mapping
                             on transcriptome.""")

    parser.add_option("--output_fastq",
                      action = "store",
                      type = "string",
                      dest = "output_fastq_filename",
                      help = """The output FASTQ file containing all reads which
                             have their mate mapping on the candidate fusion
                             genes.""")

    parser.add_option("--log",
                      action = "store",
                      type = "string",
                      dest = "output_log_filename",
                      help = """The output log file.""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_fastq_filename and
            options.input_fusions_filename and
            options.input_transcriptome_filename and
            options.output_fastq_filename
            ):
        parser.print_help()
        sys.exit(1)

    #print "Reading the list of candidate fusion genes...",options.input_fusions_filename
    fusions = [line.rstrip('\r\n').split('\t') for line in file(options.input_fusions_filename,'r').readlines() if line.rstrip('\r\n')]
    fusions = set(list(itertools.chain(*fusions)))
#    print "1->",len(fusions),fusions

    #print "Reading FASTQ file...",options.input_fastq_filename
    fid = open(options.input_fastq_filename,'r')
    r1 = []
    r2 = []
    buffer_size = 10**8
    i = 0
    j = 0
    while True:
        lines = fid.readlines(buffer_size)
        if not lines:
            break
        for line in lines:
            i = i + 1
            if i == 1:
                j = j + 1
                # this is read id
                r = line[1:-1]
                if r.endswith('/2'):
                    gc.disable()
                    r1.append(r[:-1])
                    gc.enable()
                elif r.endswith('/1'):
                    gc.disable()
                    r2.append(r[:-1])
                    gc.enable()
            elif i == 4:
                i = 0
    fid.close()
    r1 = set(r1)
    r2 = set(r2)
    x = r1.intersection(r2)
    r1.difference_update(x)
    r2.difference_update(x)
    r1 = set(("%s1" % (el,) for el in r1))
    r2 = set(("%s2" % (el,) for el in r2))
    r1.update(r2)
    mates = r1
#    print "2->",len(mates)

    # read the transcriptome mappings
    # get: read name and gene on which it maps
    # col 1: F1000050402361/1
    # col 3: tr=ENST00000449131;ge=ENSG00000167995;pn=ENSP00000399709;chr=11;str=+;len=4267
    #print "Reading transcriptome mappings...",options.input_transcriptome_filename
    remains = set(["@%s1\n" % (a_read[:-1],) if a_read.endswith('2') else "@%s2\n" % (a_read[:-1],) for a_read in read_from(options.input_transcriptome_filename, reads = mates, genes = fusions )])
    if options.output_log_filename:
        if options.output_log_filename == '-':
            fo = sys.stdin
        else:
            fo = open(options.output_log_filename,'a')
        fo.write("\n\nFiltering reads which do not have their mate-read mapping on any candidate fusion gene:\n--------------------------------------------------------\nFrom a total of %d read(s), %d read(s) are left after filtering.\n\n" % (j,len(remains),))
        fo.close()

    #
    #print 'Writing...',options.output_psl_filename
    fi = open(options.input_fastq_filename,'r')
    fo = open(options.output_fastq_filename,'w')
    i = 0
    save = False
    while True:
        lines = fi.readlines(buffer_size)
        if not lines:
            break
        new_lines = []
        for line in lines:
            i = i + 1
            if i == 1:
                save = False
                if line in remains:
                    save = True
                    remains.discard(line)
            elif i == 4:
                i = 0
            if save:
                new_lines.append(line)
        #
        if new_lines:
            fo.writelines(new_lines)
    fo.close()
    fi.close()

    #print "The end."
