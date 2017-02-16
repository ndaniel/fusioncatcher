#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given a list of reads names it removes them from a MAP BOWTIE file.



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
import optparse
import gzip
import gc

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""Given a list of reads names it removes/extracts them from a MAP BOWTIE file."""
    version="%prog 0.11 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_reads",
                      action="store",
                      type="string",
                      dest="input_reads_filename",
                      help="""The input file containing the list of reads names to be removed from the MAP BOWTIE file.""")

    parser.add_option("--input_map",
                      action="store",
                      type="string",
                      dest="input_map_filename",
                      help="""The input file in Bowtie MAP format from where the reads will be removed.""")

    parser.add_option("--output_map",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text file containing all reads from the input MAP BOWTIE file except the ones which have been removed.""")

    choices = ('remove','extract')
    parser.add_option("--operation",
                      action = "store",
                      type = "choice",
                      choices = choices,
                      dest = "kind",
                      default = "remove",
                      help = "Type of operation to be performed. The choices are "+
                             "['"+"','".join(choices)+"']. "+
                             "Default is '%default'.")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_reads_filename and
            options.input_map_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    print "Reading...",options.input_reads_filename
    gc.disable()
    offenders = frozenset([line.rstrip('\r\n') for line in file(options.input_reads_filename,'r').readlines() if line.rstrip('\r\n')])
    gc.enable()

    print "Reading...",options.input_map_filename
    print 'Writing...',options.output_filename

    fi = None
    if options.input_map_filename.lower().endswith('.gz'):
        fi = gzip.open(options.input_map_filename,'r')
    elif options.input_map_filename == "-":
        fi = sys.stdin
    else:
        fi = open(options.input_map_filename,'r')

    fo = None
    if options.output_filename.lower().endswith('.gz'):
        fo = gzip.open(options.output_filename,'w')
    elif options.output_filename == "-":
        fo = sys.stdin
    else:
        fo = open(options.output_filename,'w')

    if options.kind == 'remove':
        last = ''
        while True:
            lines = fi.readlines(10**8)
            if not lines:
                break
            bucket = []
            for line in lines:
                r = line.partition('\t')[0]
                if r == last:
                    if do:
                        gc.disable()
                        bucket.append(line)
                        gc.enable()
                elif r in offenders:
                    do = False
                    last = r
                else:
                    do = True
                    last = r
                    gc.disable()
                    bucket.append(line)
                    gc.enable()
            if bucket:
                fo.writelines(bucket)
    elif options.kind == 'extract':
        while True:
            lines = fi.readlines(10**8)
            if not lines:
                break
            bucket = []
            for line in lines:
                r = line.partition('\t')[0]
                if r == last:
                    if do:
                        gc.disable()
                        bucket.append(line)
                        gc.enable()
                elif r in offenders:
                    do = True
                    last = r
                    gc.disable()
                    bucket.append(line)
                    gc.enable()
                else:
                    do = False
                    last = r
            if bucket:
                fo.writelines(bucket)
    fo.close()
    fi.close()

    print "The end."
