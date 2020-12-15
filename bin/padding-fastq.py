#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes a FASTQ file and padds its reads which are shorter than a given size with a given nucleotide.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2020 Daniel Nicorici

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
import gc
import gzip


if __name__=='__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It takes a FASTQ file and padds its reads which are shorter than a given size with a given nucleotide."""
    version="%prog 0.12 beta              Author: Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file (in FASTQ format) containing the short reads to be processed.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTQ file containing the short reads where the reads, shorted than '--size', are padding using '--nucleotide'.""")

    parser.add_option("--nucleotide","-n",
                      action="store",
                      type="string",
                      dest="nucleotide",
                      default="N",
                      help="""The character to be use for padding. Default is '%default'.""")

    parser.add_option("--size","-s",
                      action="store",
                      type="int",
                      dest="size",
                      default = 0,
                      help="""If if this is larger than zero then all the short reads strictly shorter than this threshold will be padded. Reads longer than the threshold will be trimmed from 3 end. Default is '%default'.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")
        sys.exit(1)

    buffer_size = 10**8

    fin = None
    if options.input_filename == '-':
        fin = sys.stdin
    elif options.input_filename.endswith('.gz'):
        fin = gzip.open(options.input_filename,'r')
    else:
        fin = open(options.input_filename,'r')

    fou = None
    if options.output_filename == '-':
        fou = sys.stdout
    elif options.output_filename.endswith('.gz'):
        fou = gzip.open(options.output_filename,'w')
    else:
        fou = open(options.output_filename,'w')
    i = 0
    first = True
    pad = options.nucleotide * options.size
    qual = 'F' * options.size
    n = options.size
    next = False
    while True:
        gc.disable()
        data = fin.readlines(buffer_size)
        gc.enable()
        if not data:
            fou.writelines(data)
            break
        if first:
            # find the lowest quality score and use it for padding
            try:
                qual = min([min(data[j][:-1]) for j in xrange(len(data)) if (j+i)%4 == 3 and len(data[j])>1])
            except:
                qual = 'F'
            qual = qual * n
            first = False
        if n:
            nn = len(data)
            for j in xrange(nn):
                x = (j+i)%4
                if x == 1: # sequence
                    m = len(data[j])-1
                    if m < n:
                        d = data[j][:-1] + pad
                        data[j] = d[:n]+'\n'
                        next = True
                    elif m > n:
                        data[j] = data[j][:n]+'\n'
                        next = True
                elif x == 2:
                    data[j] = "+\n"
                elif x == 3: # quality scores
                    d = data[j][:-1] + qual
                    data[j] = d[:n]+'\n'
                    next = False
            i = i + nn
            while data:
                if data[-1].rstrip('\r\n'):
                    break
                else:
                    data.pop()
        fou.writelines(data)
    fin.close()
    fou.close()
    #
