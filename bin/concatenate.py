#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It concatenates several files.
The input files are concatenated and the result is in the output file, which is the last filename in the command line arguments.
It makes sure that the last line ends with newline character (i.e. it adds one if it is missing)!



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

import sys
import os
import gc
import gzip

def concatenate(output_file,input_files):
    """
    All the files from the list of 'input_files' are concatenated into 'output_file'.
    If the last line of a file from the 'input_files' does not end with a new line
    character than a new line character is added.
    """
    fod = None
    if output_file == '-':
        fod = sys.stdout
    elif output_file.lower().endswith('.gz'):
        fod = gzip.open(output_file,'wt')
    else:
        fod = open(output_file,'wt')

    for element in input_files:
        #print element
        fid = None
        if element.lower().endswith('.gz'):
            fid = gzip.open(element,'rt')
        else:
            fid = open(element,'rt')
        while True:
            gc.disable()
            lines = fid.readlines(10**8)
            gc.enable()
            if not lines:
                break
            lines[-1] = lines[-1].rstrip('\r\n') + '\n'
            fod.writelines(lines)
        fid.close()
    fod.close()

if __name__ == '__main__':

    #command line parsing
    # initializing
    # reading command line arguments
    filenames = sys.argv

    if len(filenames) > 1:
        if filenames[1] == '-f':
            infile = [line.rstrip('\r\n') for line in file(filenames[2],'r').readlines() if line.rstrip('\r\n')]
            outfile = filenames[-1]
        else:
            infile = filenames[1:-1]
            outfile = filenames[-1]
    else:
        print >>sys.stderr,"""It concatenates several files.
The input files are concatenated and the result is in the output file, which is the last filename in the command line arguments.
It makes sure that the last line ends with newline character (i.e. it adds one if it is missing)!

If the command line option '-f' is specified then only on input file is given which will contain on each line the input files which should be concatenated!
"""
        print >>sys.stderr,"ERROR: Not enough arguments!"
        sys.exit(1)
    #print "Concatenating the following files:"
    #for el in infile:
    #    print ' - ',el
    #print "into file:", outfile

    concatenate(outfile,infile)
    #
