#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given input FASTQ file it extracts all the reads ids.



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


import sys
import os
import optparse
import gc



if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""Given input FASTQ file it extracts all the reads ids."""
    version="%prog 0.10 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ format.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text file containinf the reads ids.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    fin = open(options.input_filename,'r')
    fout = open(options.output_filename,'w')
    chunk = []
    buffer_size = 10**8
    i = 0
    while True:
        lines = fin.readlines(buffer_size)
        if not lines:
            break
        n = len(lines)
        lines = [lines[j][1:] for j in xrange(n) if (j+i) % 4 == 0]
        i = i + n
        fout.writelines(lines)
    fin.close()
    fout.close()

    print "The end."
