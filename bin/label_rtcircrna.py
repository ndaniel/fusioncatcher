#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It labels the fusions with label such as reciprocal for reciprocal fusion genes and rt_circ_rna for the reciprocal fusion genes which are also a readthrough.

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





#
#
#
def rtcircrna(input_file,
             output_file,
             verbose = True):


    if verbose:
        print >>sys.stderr,"Reading the list of fusion genes..."
    # read the list of fusion genes and their chromosomal positions
    data = [line.rstrip('\r\n').split('\t') for line in file(input_file,'r').readlines() if line.rstrip('')]
    header = data.pop(0)
    if data:
        g = [(line[0],line[1],True if line[2].find("readthrough")!=-1 else False) for line in data]
        r = set([(e[0],e[1]) for e in g])
        reciprocal = set([e for e in r if (e[1],e[0]) in r])
        rtcircrna = [(e[0],e[1]) for e in g if ((e[0],e[1]) in reciprocal) and e[2]]
        rtcircrna = set(rtcircrna+[(e[1],e[0]) for e in rtcircrna])
        
        r = [header]
        for e in data:
            t = e
            if (e[0],e[1]) in reciprocal:
                x = ',reciprocal' if e[2] else 'reciprocal'
                e[2] = e[2] + x
            if (e[0],e[1]) in rtcircrna:
                x = ',rt_circ_rna' if e[2] else 'rt_circ_rna'
                e[2] = e[2] + x
            r.append(t)
        file(output_file,'w').writelines(['\t'.join(line)+'\n' for line in r])
    else:
        file(output_file,'w').write('\t'.join(header)+'\n')



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

    description = """It labels the fusions with label such as reciprocal for reciprocal fusion genes and rt_circ_rna for the reciprocal fusion genes which are also a readthrough."""

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


    #
    rtcircrna(
             options.input_filename,
             options.output_filename,
             options.verbose)


if __name__ == '__main__':
    main()
    
    
#    
