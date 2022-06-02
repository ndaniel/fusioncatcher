#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It shuffles the reads from the two input FASTQ files into one output FASTQ file.



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

def reads_from_fastq_file(f_name, size_buffer = 5 * (10**7)):
    fid = None
    if f_name == '-':
        fid = sys.stdin
    elif f_name.lower().endswith('.gz'):
        fid = gzip.open(f_name,'r')
    else:
        fid = open(f_name,'r')
    piece = [None,None,None,None]
    i = 0
    while True:
        gc.disable()
        lines = fid.readlines(size_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            i = i + 1
            piece[i-1] = line
            if i == 4:
                piece[2] = '+\n'
                yield piece
                piece = [None,None,None,None]
                i = 0
    fid.close()



if __name__=='__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It shuffles the reads into one single file (read-A/1 from file reads1.fq, read-A/2 from file reads2.fq, read-B/1 from file reads1.fq ...). """
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input_1","-f",
                      action = "store",
                      type = "string",
                      dest = "input_1_filename",
                      help = """The first FASTQ input file containing the short reads to be interleaved.""")

    parser.add_option("--input_2","-r",
                      action = "store",
                      type = "string",
                      dest = "input_2_filename",
                      help = """The second FASTQ input file containing the short reads to be interleaved.""")

    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output FASTQ file where the short reads from the two input FASTQ files are interleaved.""")



    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_1_filename and
            options.input_2_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")

    print "Starting..."

    ou = open(options.output_filename,'w')
    data = []
    i = 0
    limit = 10**5
    for (it1,it2) in itertools.izip_longest(reads_from_fastq_file(options.input_1_filename),
                                            reads_from_fastq_file(options.input_2_filename)):
        if it1:
            gc.disable()
            data.extend(it1)
            gc.enable()
            i = i + 1
        if it2:
            gc.disable()
            data.extend(it2)
            gc.enable()
            i = i + 1
        if i > limit:
            ou.writelines(data)
            data = []
            i = 0
    if data:
        ou.writelines(data)
        data = []
    ou.close()

    print "End."
