#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes a sorted bed file (by coordinates) and merge the intervals which are overlapping.



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


import os
import sys
import optparse
import gzip


def give(fi):
    fin = None
    if fi == '-':
        fin = sys.stdin
    elif fi.lower().endswith('.gz'):
        fin = gzip.open(fi,'r')
    else:
        fin = open(fi,'r')

    buffer_size = 10**8
    while True:
        lines = fin.readlines(buffer_size)
        if not lines:
            break
        
        for line in lines:
            if line.rstrip("\r\n"):
                x = line.rstrip("\r\n").split("\t")
                yield (x[0],int(x[1]),int(x[2]))
    yield (0,0,0) # mark for last record
    fin.close()

############################################################################
if __name__=='__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It takes a sorted bed file (by coordinates) and merge the intervals which are overlapping."""
    version="%prog 0.11 beta              Author: Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input (sorted) BED file.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output (sorted) BED file.""")



    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")
        sys.exit(1)

    #

    fou = None
    if options.output_filename == '-':
        fou = sys.stdout
    elif options.output_filename.lower().endswith('.gz'):
        fou = gzip.open(options.output_filename,'w')
    else:
        fou = open(options.output_filename,'w')

    buf = []
    last = (0,0,0)
    for line in give(options.input_filename):

        if line[0] == last[0] and line[1] < last[2]:
            last = (last[0],last[1],max(line[2],last[2]))
            continue
            
        if last[0]:
            buf.append("%s\t%d\t%d\n" % (last[0],last[1],last[2]))
            if len(buf) > 10000:
                fou.writelines(buf)
                buf = []
        elif not line[0]:
            break

        last = line[:]

    if buf:
        fou.writelines(buf)


    fou.close()
    #
