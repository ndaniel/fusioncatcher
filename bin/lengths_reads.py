#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It writes into a text files the lengths (sorted in descending order) of the reads found in a FASTQ Solexa file.



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
import sys
import os
import optparse
import gc
import gzip

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It writes into a text files the lengths (sorted in descending order) of the reads found in a FASTQ Solexa file."""
    version="%prog 0.10 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ Solexa file (also given thru stdin or as gzipped file).""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text file containg the unique lengths of the reads found in the input file. The unique lengths are sorted in descending order.""")

    parser.add_option("--counts","-c",
                      action="store",
                      type="string",
                      dest="counts_filename",
                      help="""The output text file containg the counts of reads found in the input file.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("Input and output files should be specified!")
        sys.exit(1)


    fin = None
    if options.input_filename == '-':
        fin = sys.stdin
    elif options.input_filename.lower().endswith('.gz'):
        fin = gzip.open(options.input_filename,'r')
    else:
        fin = open(options.input_filename,'r')
    i = 0
    l = set()
    m = -1
    buffer = 10**8
    while True:
        gc.disable()
        lines = fin.readlines(buffer)
        gc.enable()
        if not lines:
            break
        gc.disable()
        l.update(set([len(lines[j]) for j in xrange(len(lines)) if (j+i)%4==1 and m != len(lines[j]) ]))
        gc.enable()
        m = max(l)
        i=i+len(lines)
    fin.close()

    file(options.output_filename,'w').writelines([str(line-1)+'\n' for line in sorted(list(l),reverse=True)])
    if options.counts_filename:
        file(options.counts_filename,'w').write("%d" % ((i+1)/4,))
    #
