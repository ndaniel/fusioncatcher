#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It un-interleaves reads from an input FASTQ files where the reads are interleaved into two output FASTQ file.


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


import os
import sys
import optparse
import itertools
import gc
import gzip


if __name__=='__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It un-interleaves reads from an input FASTQ files where the reads are interleaved into two output FASTQ file."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The first FASTQ input file containing the short reads which are interleaved.""")


    parser.add_option("--forward","-f",
                      action = "store",
                      type = "string",
                      dest = "output_forward_filename",
                      help = """The output FASTQ file where the forward short reads from the input FASTQ.""")

    parser.add_option("--reverse","-r",
                      action = "store",
                      type = "string",
                      dest = "output_reverse_filename",
                      help = """The output FASTQ file where the reverse short reads from the input FASTQ.""")



    (options, args) = parser.parse_args()

    # validate options
    if not (options.output_forward_filename and
            options.output_reverse_filename and
            options.input_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")

    #print "Starting..."

    fi = None
    if options.input_filename == '-':
        fi = sys.stdin
    elif options.input_filename.lower().endswith('.gz'):
        fi = gzip.open(options.input_filename,'r')
    else:
        fi = open(options.input_filename,'r')
        
    ff = None
    if options.output_forward_filename.lower().endswith('.gz'):
        ff = gzip.open(options.output_forward_filename,'w')
    else:
        ff = open(options.output_forward_filename,'w')

    fr = None
    if options.output_reverse_filename.lower().endswith('.gz'):
        fr = gzip.open(options.output_reverse_filename,'w')
    else:
        fr = open(options.output_reverse_filename,'w')

    i = -1
    size_buffer = 10**8

    while True:
        gc.disable()
        lines = fi.readlines(size_buffer)
        gc.enable()
        if not lines:
            break
        f = []
        r = []
        gc.disable()
        for aline in lines:
            i = i + 1
            if i == 8:
                i = 0
            if i < 4:
                f.append(aline)
            else:
                r.append(aline)
        gc.enable()
        ff.writelines(f)
        fr.writelines(r)
    
    if options.input_filename != '-':
        fi.close()
    ff.close()
    fr.close()

    #print "End."
