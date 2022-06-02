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
            if k != x[0]:
                if bucket:
                    yield (k,tuple(bucket))
                    bucket = []
                k = x[0]
            bucket.append((x[1],x[2]))
    if bucket:
        yield (k,tuple(bucket))




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


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file where the frame predictions are written. """)

    parser.add_option("-x","--extra",
                      action = "store",
                      type = "string",
                      dest = "extra_filename",
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

    fin = None
    if options.input_filename == '-':
        fin = sys.stdin
    else:
        fin = open(options.input_filename,'r')

    fou = None
    if options.output_filename == '-':
        fou = sys.stdout
    else:
        fou = open(options.output_filename,'w')

    fex = None
    if options.extra_filename:
        fex = open(options.extra_filename,"w") 

    r = []
    ex = []
    for e in giveme(fin):
        k = e[0]
        t = e[1]
        dic = dict()
        for x in t:
            if not dic.has_key(x[1]):
                dic[x[1]] = set()
            dic[x[1]].add(x[0])
        ks = sorted(dic.keys())
        if len(ks) > 2:
            print >>sys.stder,"ERROR: not expected so many keys!"
            sys.exit(1)
        if len(ks) == 2:
            ka = ks[0]
            kb = ks[1]
            da = dic[ka]
            db = dic[kb]
            intersect = da.intersection(db)
            lda = len(da)
            ldb = len(db)
            if (not intersect) and lda < 4 and ldb < 4:
                r.append(k+"\n")
                if len(r) > 100000:
                    fou.writelines(r)
                    r = []
                u1 = '>1'
                if lda > 2:
                    u1 = '>2'
                if lda == 1:
                    u1 = da.pop()
                u2 = ">1"
                if ldb > 2:
                    u2 = '>2'
                if ldb == 1:
                    u2 = db.pop()
                ex.append("%s\t%s\t%s\n" % (k,u1,u2))
                if len(ex) > 100000:
                    fex.writelines(ex)
                    ex = []


    if r:
        fou.writelines(r)
        
    if ex:
        fex.writelines(ex)


    if options.input_filename != '-':
        fin.close()

    if options.output_filename != '-':
        fou.close()

    if options.extra_filename:
        fex.close()

if __name__ == '__main__':
    main()
    
    
#    
