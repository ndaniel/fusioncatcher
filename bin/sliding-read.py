#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It does a sliding window for each raed such that several shorter reads are generate from the input read.


Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

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
    description = """It does a sliding window for each raed such that several shorter reads are generate from the input read."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The first FASTQ input file containing the short reads which are interleaved.""")

    parser.add_option("--window","-w",
                      action="store",
                      type="int",
                      dest="window_size",
                      help="""The size of window.""")

    parser.add_option("--step","-s",
                      action="store",
                      type="int",
                      dest="step_size",
                      help="""The size of window.""")


    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output FASTQ file.""")



    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename and
            options.step_size and 
            options.window_size
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
        
    fo = None
    if options.output_filename == '-':
        fo = sys.stdout
    elif options.output_filename.lower().endswith('.gz'):
        fo = gzip.open(options.output_filename,'w')
    else:
        fo = open(options.output_filename,'w')


    size_buffer = 10**8

    step = options.step_size
    window = options.window_size
    bucket = []
    r = []
    while True:
        gc.disable()
        lines = fi.readlines(size_buffer)
        gc.enable()
        if not lines:
            break
        gc.disable()
        for aline in lines:
            bucket.append(aline)
            if len(bucket) == 4:
                bucket[1] = bucket[1].rstrip("\r\n")
                bucket[3] = bucket[3].rstrip("\r\n")
                n = len(bucket[1])
                for j in xrange(0,n-window+1,step):
                    r.append("%s%s\n+\n%s\n" % (bucket[0],bucket[1][j:j+window],bucket[3][j:j+window]))
                bucket = []
        gc.enable()
        if r:
            fo.writelines(r)
            r = []
            
    if r:
        fo.writelines(r)

    
    if options.input_filename != '-':
        fi.close()

    if options.output_filename != '-':
        fo.close()

    #print "End."
