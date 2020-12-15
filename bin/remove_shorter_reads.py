#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes a FASTQ file and removes all reads which are shorter than a given threshold.



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
import string
import gzip
import gc


if __name__=='__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It takes a FASTQ file and removes all reads which are shorter strictly than a given threshold."""
    version = "%prog 0.10 beta              Author: Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file (in FASTQ format) containing the short reads to be processed.""")

    parser.add_option("--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output FASTQ file containing the short reads where the reads which are shorter then a given threshold are removed.""")

    parser.add_option("--threshold",
                      action = "store",
                      type = "int",
                      dest = "threshold",
                      default = 0,
                      help = """If if this is larger than zero then all the short reads strictly shorter than this threshold are removed. Default is %default.""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")
        sys.exit(1)

    size_read_buffer = 10**8
    threshold = options.threshold
    if threshold > 0:
        # input
        fid = None
        f_name = options.input_filename
        if f_name == '-':
            fid = sys.stdin
        elif f_name.lower().endswith('.gz'):
            fid = gzip.open(f_name,'r')
        else:
            fid = open(f_name,'r')
        # output
        f_name = options.output_filename
        fou = None
        if f_name == '-':
            fou = sys.stdout
        elif f_name.lower().endswith('.gz'):
            fou = gzip.open(f_name,'w')
        else:
            fou = open(f_name,'w')

        piece = [None,None,None,None]
        i = 0
        bucket = []
        while True:
            gc.disable()
            lines = fid.readlines(size_read_buffer)
            gc.enable()
            if not lines:
                break
            for a_line in lines:
                piece[i] = a_line
                i = i + 1
                if i == 4:
                    if len(piece[1]) > threshold:
                        gc.disable()
                        bucket.append("%s%s+\n%s"%(piece[0],piece[1],piece[3]))
                        gc.enable()
                    piece = [None,None,None,None]
                    i = 0
            if bucket:
                fou.writelines(bucket)
                bucket = []
        if bucket:
            fou.writelines(bucket)
            bucket = []
        fid.close()
        fou.close()
    else:
        # input
        fid = None
        f_name = options.input_filename
        if f_name == '-':
            fid = sys.stdin
        elif f_name.lower().endswith('.gz'):
            fid = gzip.open(f_name,'r')
        else:
            fid = open(f_name,'r')
        # output
        f_name = options.output_filename
        fou = None
        if f_name == '-':
            fou = sys.stdout
        elif f_name.lower().endswith('.gz'):
            fou = gzip.open(f_name,'w')
        else:
            fou = open(f_name,'w')

        piece = [None,None,None,None]
        i = 0
        bucket = []
        while True:
            gc.disable()
            lines = fid.readlines(size_read_buffer)
            gc.enable()
            if not lines:
                break
            for a_line in lines:
                piece[i] = a_line
                i = i + 1
                if i == 4:
                    gc.disable()
                    bucket.append("%s%s+\n%s"%(piece[0],piece[1],piece[3]))
                    gc.enable()
                    piece = [None,None,None,None]
                    i = 0
            if bucket:
                fou.writelines(bucket)
                bucket = []
        if bucket:
            fou.writelines(bucket)
            bucket = []
        fid.close()
        fou.close()
