#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It prepares the fines needed to analyze one candidate fusion gene at a time.



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
import os
import gc
import optparse


######################################################################
######################################################################
######################################################################
######################################################################


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It prepares the fines needed to analyze one candidate fusion gene at a time."""
    version = "%prog 0.14 beta"

    parser = optparse.OptionParser(usage       = usage,
                                   description = description,
                                   version     = version)

    parser.add_option("-i","--fusions",
                      action = "store",
                      type = "string",
                      dest = "input_fusions",
                      help = """The input list of fusions for which the supporting reads should be extracted.""")

    parser.add_option("-r","--reads",
                      action = "store",
                      type = "string",
                      dest = "input_supporting_reads",
                      help = """The file containing the supporting paired-reads.""")

    parser.add_option("-1","--output-list-files-1",
                      action = "store",
                      type = "string",
                      dest = "output_list_filename1",
                      help = """It is a list of files containing the extracted information, regarding preliminary candidate fusion genes.""")

    parser.add_option("-2","--output-list-files-2",
                      action = "store",
                      type = "string",
                      dest = "output_list_filename2",
                      help = """It is a list of files containing the extracted information, regarding reads ids which support the preliminary candidate fusion genes.""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.input_fusions and
            options.input_supporting_reads and
            options.output_list_filename1 and 
            options.output_list_filename2
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    # read the fusions of interest
    f = set([tuple(line.rstrip("\r\n").split("\t")[:2]) for line in file(options.input_fusions,"r") if line.rstrip("\r\n")])
    
    # reads the file with the reads
    r = [line.rstrip("\r\n").split("\t") for line in file(options.input_supporting_reads,"r") if line.rstrip("\r\n")]
    r.pop(0)
    
    o1 = []
    o2 = []
    i = -1
    for line in r:
        if (line[2],line[3]) in f:
            i = i + 1
            r1 = options.output_list_filename1 + "."+str(i)
            r2 = options.output_list_filename2 + "."+str(i)
            file(r1,"w").write("%s\t%s" % (line[2],line[3]))
            d = line[5].split(",")
            d2 = []
            for e in d:
                d2.append(e+"/1\n")
                d2.append(e+"/2\n")
            file(r2,"w").writelines([e for e in d2])
            o1.append(r1+"\n")
            o2.append(r2+"\n")
    file(options.output_list_filename1,"w").writelines(o1)
    file(options.output_list_filename2,"w").writelines(o2)
    #
