#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It removes the pairs of reads where one or both of the reads are shorter than a 
given threshold.



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
import optparse


######################################################################
######################################################################
######################################################################
######################################################################





if __name__ == '__main__':

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options] --input <tab separated file>  --output <tab separated file>"

    description = """It removes the pairs of reads where one or both of the reads are shorter than a given threshold."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2020 Daniel Nicorici

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
                      help = """Tab separeted text file such that one line has two reads which are paired (i.e. generated via 'paste - - - - - - - -').""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """Tab separeted text file such that one line has two reads which are paired (i.e. generated via 'paste - - - - - - - -').""")


    parser.add_option("-l","--length",
                      action = "store",
                      type = "int",
                      dest = "length",
                      default = 25,
                      help = """The minimum length of a read. All reads shorter than this will be removed. Default is %default.""")

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

    l = options.length - 1
    size_buffer = 10**8
    
    while True:
        lines = fi.readlines(size_buffer)
        if not lines:
            break
        lines = [line.split("\t") for line in lines if line]
        lines = ["\t".join(line) for line in lines if len(line[1]) > l and len(line[5]) > l ]
        fo.writelines(lines)
        
        
    if options.input_filename != '-':
        fi.close()
    if options.output_filename != '-':
        fo.close()
    
    
    
#
