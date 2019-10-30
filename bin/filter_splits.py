#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It pre-filtering of splitted reads.

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2019 Daniel Nicorici

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


"""
Input looks like this:


"""

import os
import sys
import optparse
import gzip

def prepare(x):
    y = x.rstrip("\r\n").split("\t")
    return y[0][:-2]+"\t"+y[1]


################################################################################
################################################################################
################################################################################
def main():

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options]"

    description = """It pre-filtering of splitted reads."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2019 Daniel Nicorici

"""



    version = "%prog 0.04 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
                            )



    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file containing the fusion (chromosomal) coordinates for each fusion genes.""")
# reads_filtered_not-mapped-genome_transcriptome_trim2_3end.map.all.ex

    parser.add_option("-x","--extra",
                      action = "store",
                      type = "string",
                      dest = "extra_filename",
                      help = """The input file containing the fusion (chromosomal) coordinates for each fusion genes.""")
# reads_filtered_transcriptome_trim2.txt


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file where the frame predictions are written. """)

    parser.add_option("-m","--map",
                      action = "store",
                      type = "string",
                      dest = "map_filename",
                      help = """The output file where the frame predictions are written. """)

    parser.add_option("-q", "--quiet",
                      action = "store_false",
                      dest = "verbose",
                      default = True,
                      help = "Do not print status messages to stdout.")

    parser.add_option("-a","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (
            options.input_filename and
            options.extra_filename and 
            options.output_filename 
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    fin = open(options.input_filename,'r')

    fou = open(options.output_filename,'w')

    keep = set([prepare(e) for e in file(options.extra_filename,"r") if e.rstrip("\r\n")])

    size_buffer = 10**8
    o = list()
    while True:
        lines = fin.readlines(size_buffer)
        if not lines:
            break
        lines = [e.rstrip("\r\n").split("\t") for e in lines if e.rstrip("\r\n")]
        for e in lines:
            if e[1] != e[2]:
                g1 = e[0][:-2]+"\t"+e[1]
                g2 = e[0][:-2]+"\t"+e[2]
#                if g1 in keep:
#                    o.add((e[0],e[2]))
#                elif g2 in keep:
#                    o.add((e[0],e[1]))
                if (g1 in keep) or (g2 in keep):
                    o.append(e[0]+"\n")
                if len(o) > 10000:
                    fou.writelines(sorted(set(o)))
                    o = []


    if o:
        fou.writelines(sorted(set(o)))

    fin.close()
    fou.close()


#    if options.extra_filename:
#        fex.close()

if __name__ == '__main__':
    main()
    
    
#    
