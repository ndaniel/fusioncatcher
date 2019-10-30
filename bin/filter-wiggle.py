#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a list of chromosomal coordinates of fusion genes and filters out the ones that might be an artefact introduced by up-stream wiggle method.

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

import os
import sys
import optparse
import gzip







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

    description = """It takes as input a list of chromosomal coordinates of fusion genes and filters out the ones that might be an artefact introduced by up-stream wiggle method."""

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


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
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
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    d = [l.rstrip("\r\n").split("\t") for l in file(options.input_filename,"r") if l.rstrip("\r\n")]
    #
    r = []
    wiggle = 20
    if d:
        header = d.pop(0)
        r.append(header)
        if d:
            # process
            bucket = dict()
            for e in d:
                g1 = e[0] # gene 1
                g2 = e[1] # gene 2
                e1 = e[10] # ensembl 1
                e2 = e[11] # ensembl 2
                p1 = tuple(e[8].split(":")) # position 1
                p2 = tuple(e[9].split(":")) # position 2
                cp = int(e[4]) # count pairs
                cr = int(e[5]) # count reads
                lab = e[2].find("exon-exon")
                k = e1+e2
                if bucket.has_key(k):
                    flag = True
                    for x in bucket[k]:
                        if lab == -1 and x[0] >= cp and x[1] >= cr and x[2][0] == p1[0] and x[3][0] == p2[0] and x[2][2] == p1[2] and x[3][2] == p2[2] and abs(int(x[2][1])-int(p1[1])) <= wiggle and abs(int(x[3][1])-int(p2[1])) <= wiggle:
                            flag = False
                            break
                    if flag:
                        r.append(e)
                        bucket[k].append((cp,cr,p1,p2))
                else:
                    bucket[k] = [(cp,cr,p1,p2)]
                    r.append(e)
            
    else:
        r = d
    file(options.output_filename,"w").writelines(['\t'.join(e)+'\n' for e in r])



if __name__ == '__main__':
    main()
    
    
#    
