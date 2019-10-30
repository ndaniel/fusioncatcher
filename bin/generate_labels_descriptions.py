#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the list of banned candidate fusion genes. This list is hard coded inhere.



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

def expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It generates the description of labels which are used to described the found fusion genes."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the descriptions of labels is generated. Default is '%default'.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    pipeline_path = os.path.dirname(expand(sys.argv[0]))
    manual = os.path.abspath(os.path.join(pipeline_path,"..","doc","manual.md"))

    #
    #
    #

    print "Reading the manual..."
    d = [line for line in file(manual,'r').readlines()]
    
    tables = []
    flag = False
    for line in d:
        if line.strip().lower().startswith('table 1'):
            flag = True
        if flag and line.strip().lower().startswith('#'):
            flag = False
        if flag:
            tables.append(line)

    print "writing the labels description..."
    file(os.path.join(options.output_directory,'final-list_candidate-fusion-genes.caption.md.txt'),'w').writelines(tables)
    #
