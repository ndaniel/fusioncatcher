#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a list of pairs of genes and it orders it and removes the duplicates.

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

def fix(fin,fou,tab=" "*8):
    d = [line for line in file(fin,'r').readlines() if line.strip()]
    dhash = [line for line in d if line.strip().startswith("#")]
    d = [line.upper().rstrip('\r\n').replace("'","").replace('"','').replace("[","").replace("]","").replace(" ","").partition(",") for line in d if not line.strip().startswith("#")]
    d = [(line[0],line[2].partition(",")[0].strip()) for line in d]
    d = sorted(set([(line[0],line[1]) if line[0]<line[1] else (line[1],line[0]) for line in d]))
    d = ["%s['%s','%s'],\n" % (tab,line[0],line[1]) for line in d if line[0] != line[1]]
    d[-1] = d[-1][:-2]+'\n'
    d.extend(dhash)
    file(fou,'w').writelines(d)


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

    description = """It takes as input a list of pairs of genes and it orders it and removes the duplicates."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2022 Daniel Nicorici

"""



    version = "%prog 0.02 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
                            )



    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input list of gene pairs.""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The sorted output of gene pairs. """)


    parser.add_option("-a","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")


    #
    fix(options.input_filename,
        options.output_filename)


if __name__ == '__main__':
    main()

#
