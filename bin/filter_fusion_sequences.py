#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It filters out the fusion sequences that are given as a input.

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2021 Daniel Nicorici

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
import math


def give(s,le=30):
    r = s.upper().rstrip("\r\n").split("*")
    x = ""
    if len(r) == 2:
        x = r[0][-le:]+"*"+r[1][0:le]
    return x


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

    description = """It filters out the fusion sequences that are given as a input."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2021 Daniel Nicorici

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

    parser.add_option("-b","--banned",
                      action = "store",
                      type = "string",
                      dest = "banned_filename",
                      help = """The input file containing the banned fusion sequences.""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file where the frame predictions are written. """)

    parser.add_option("-r","--removed",
                      action = "store",
                      type = "string",
                      dest = "removed_filename",
                      help = """The output file where the removed fusions are written. """)

    parser.add_option("-w","--window",
                      action = "store",
                      type = "int",
                      dest = "window",
                      default = 30,
                      help = """The length of the window. Default is %default.""")

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
            options.banned_filename and
            options.output_filename and
            options.removed_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")


    #
    ban = set([give(e.rstrip("\r\n").upper(),options.window) for e in file(options.banned_filename,"r") if e.rstrip("\r\n") and not e.startswith(">")])
    
    f = [e.rstrip("\r\n").split("\t") for e in file(options.input_filename,"r") if e.rstrip("\r\n")]
    rem = []
    res = []
    if f:
        h = f.pop(0)
        rem.append(h)
        res.append(h)

    for e in f:
        found = False
        if give(e[14],options.window) in ban:
            rem.append(e)
        else:
            res.append(e)
    
    file(options.output_filename,"w").writelines(["\t".join(e)+"\n" for e in res])
    file(options.removed_filename,"w").writelines(["\t".join(e)+"\n" for e in rem])



if __name__ == '__main__':
    main()
    
    
#    
