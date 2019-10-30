#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Add the custom human genes which are known to have random inserts at fusion junction.



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
import sys
import os
import optparse
import symbols

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Add the custom human genes which are known to have random inserts at fusion junction."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the list of banned candidate fusion genes is generated, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the list of banned candidate fusion genes is generated. Default is '%default'.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    data = sorted([line.rstrip("\r\n") for line in file(os.path.join(options.output_directory,'ig_loci.txt'),"r")])
    data.append("ENSG00000071564") # TCF3
    data.append("ENSG00000126746") # ZNF384
    data.append("ENSG00000066468") # FGFR2
    data.append("ENSG00000122025") # FLT3

    data = sorted(set([line+"\n" for line in data]))
    file(os.path.join(options.output_directory,'gap_fusions.txt'),'w').writelines(data)
    #
    #
    #


