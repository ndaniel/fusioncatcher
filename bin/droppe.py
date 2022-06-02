#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It removes the pairs of reads except those that are given in a separate input list.


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


import sys
import optparse
import os


######################################################################
######################################################################
######################################################################
######################################################################

def giveme(f,size_buffer=10**8):
    while True:
        lines = f.readlines(size_buffer)
        if not lines:
            break
        for line in lines:
            yield line.split("\t",1)

if __name__ == '__main__':

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options] --input <tab separated file>  --output <tab separated file>"

    description = """It removes the pairs of reads except those that are given in a separate input list."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2022 Daniel Nicorici

"""

    version = "%prog 0.98 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
)



    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """Tab separeted text file such that one line has two reads which are paired (i.e. generated via 'paste - - - -').""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """Tab separeted text file such that one line has two reads which are paired (i.e. generated via 'paste - - - -').""")


    parser.add_option("-k","--keep",
                      action = "store",
                      type = "string",
                      dest = "keep_filename",
                      help = """Reads ids that should not be removed. Default is %default.""")

    parser.add_option("-f","--info",
                      action = "store",
                      type = "string",
                      dest = "info_filename",
                      help = """Counts of Reads that where removed. Default is %default.""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    #
    fi = sys.stdin
    if options.input_filename != '-':
        fi = open(options.input_filename,'r')

    fo = sys.stdout
    if options.output_filename != '-':
        fo = open(options.output_filename,'w')

    fk = None
    keep = set()
    if options.keep_filename:
        if os.path.isfile(options.keep_filename):
            keep = [e.rstrip("\r\n") for e in file(options.keep_filename,"r") if e.rstrip("\r\n")]
            keep = set([e[:-2] if e[:-1].endswith("/1") or e.endswith("/2") else e for e in keep])

    size_buffer = 10**8
    
    r = ''
    s = ''
    bucket = []
    i = 0
    sc = 0
    if keep:
        for e in giveme(fi):
            i = i + 1
            if r:
                if r[:-1] == e[0][:-1] and r != e[0] and r[-2]=="/":
                    # we have a pair
                    if r[1:-2] not in keep:
                        r = ''
                        s = ''
                        continue
                bucket.append("%s\t%s" % (r,s))
                if len(bucket) > 100000:
                    fo.writelines(bucket)
                    sc = sc + len(bucket)
                    bucket = []
                r = e[0]
                s = e[1]
            else:
                r = e[0]
                s = e[1]
    else:
        for e in giveme(fi):
            i = i + 1
            if r:
                if r[:-1] == e[0][:-1] and r != e[0] and r[-2]=="/":
                    # we have a pair
                    r = ''
                    s = ''
                    continue
                bucket.append("%s\t%s" % (r,s))
                if len(bucket) > 100000:
                    fo.writelines(bucket)
                    sc = sc + len(bucket)
                    bucket = []
                r = e[0]
                s = e[1]
            else:
                r = e[0]
                s = e[1]

    if r:
        bucket.append("%s\t%s" % (r,s))
        sc = sc + len(bucket)

    if bucket:
        fo.writelines(bucket)

    if options.input_filename != '-':
        fi.close()
    if options.output_filename != '-':
        fo.close()

    if options.info_filename:
        file(options.info_filename,"a").write("\n\n\nFrom a total of %d unmapped reads, %d (%.4f%%) were removed because are (unmapped) pairs.\n\n\n" % (i,i-sc,float(i-sc)/float(i)*100))
    
    
#
