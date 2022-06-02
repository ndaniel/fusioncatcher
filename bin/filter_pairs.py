#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It pre-filtering of splitted reads.

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


"""
Input looks like this:
0/2	ENSG00000157554	ENSG00000184012
0/2	ENSG00000184012	ENSG00000157554
1/2	ENSG00000184012	ENSG00000157554
2/2	ENSG00000157554	ENSG00000184012
2/2	ENSG00000184012	ENSG00000157554

"""

import os
import sys
import optparse
import gzip

def giveme(fi,buffer_size=10**8):
    bucket = []
    k = ''
    while True:
        lines = fi.readlines(buffer_size)
        if not lines:
            break
        for e in lines:
            x = e.rstrip("\r\n").split("\t")
            r = x[0].partition("/")[0]
            if k != r:
                if bucket:
                    yield bucket
                    bucket = []
                k = r
            bucket.append(x)
    if bucket:
        yield bucket




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
Copyright (c) 2009-2022 Daniel Nicorici

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


    parser.add_option("-g","--output_good",
                      action = "store",
                      type = "string",
                      dest = "output_good_filename",
                      help = """The output file where the frame predictions are written. """)

    parser.add_option("-b","--output_bad",
                      action = "store",
                      type = "string",
                      dest = "output_bad_filename",
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
            options.output_good_filename and 
            options.output_bad_filename 
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    fin = None
    if options.input_filename == '-':
        fin = sys.stdin
    else:
        fin = open(options.input_filename,'r')

    fgou = open(options.output_good_filename,'w')
    fbou = open(options.output_bad_filename,'w')

#    fex = None
#    if options.extra_filename:
#        fex = open(options.extra_filename,"w") 

    r_bad = []
    r_good = []
    for e in giveme(fin):
        le = len(e)
        if le == 1:
            continue
        pair1 = [u for u in e if u[0].endswith("/1")]
        pair2 = [u for u in e if u[0].endswith("/2")]
        if not (pair1 and pair2):
            continue
        bad = True
        if le < 7:
            r1 = pair1[0][0]
            r2 = pair2[0][0]
            r1gs = set()
            for u in pair1:
                r1gs.add(u[1])
                r1gs.add(u[2])
            r2gs = set()
            for u in pair2:
                r2gs.add(u[1])
                r2gs.add(u[2])
            x = [v for v in sorted(r1gs.intersection(r2gs)) if not v.startswith(">")]
            #y = r1gs.union(r2gs)
            if x and len(x) == 2: #and len(y) != 1:
                bad = False
                r_good.append("%s\n%s\n" % (r1,r2))
                if len(r_good) > 100000:
                    fgou.writelines(r_good)
                    r_good = []
        if bad:
            r_bad.append("%s\n%s\n" % (r1,r2))
            if len(r_bad) > 100000:
                fbou.writelines(r_bad)
                r_bad = []

    if r_good:
        fgou.writelines(r_good)
        
    if r_bad:
        fbou.writelines(r_bad)


    if options.input_filename != '-':
        fin.close()

    fgou.close()

    fbou.close()


#    if options.extra_filename:
#        fex.close()

if __name__ == '__main__':
    main()
    
    
#    
