#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It counts the fragments for the preliminary fusions genes (only in case that they were fragmented using fragment_fastq.py).


Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

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




if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It counts the fragments for the preliminary fusions genes (only in case that they were fragmented using fragment_fastq.py). The preliminary fusions which have strictly less than the minimum value will be written in the output file."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

    parser.add_option("--fusion-reads","-f",
                      action = "store",
                      type = "string",
                      dest = "fusion_reads",
                      help = "Input file containing preliminary fusion genes and the supporting paired-reads. ")

    parser.add_option("--minimum","-m",
                      action = "store",
                      type = "int",
                      dest = "minimum",
                      default = 0,
                      help = "Preliminary fusion genes which have the number of paired-reads strictly less than "+
                             "this value will be written in the output. "+
                             "Default is '%default'.")

    parser.add_option("--fragments","-g",
                      action = "store",
                      type = "string",
                      dest = "fragments",
                      help = "The output file containing the preliminary fusion genes which have supporting paired-reads stricly less the given minimum value.")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.fusion_reads and
            options.fragments
            ):
        parser.print_help()
        parser.error("Missing argument(s)!")
        sys.exit(1)


    #
    #
    #
    fusions = [line.rstrip('\r\n').split('\t') for line in file(options.fusion_reads,'r').readlines() if line.rstrip('\r\n')]
    fusions.pop(0) # remove header

    fragments = []
    for line in fusions:
        o = set(line[5].split(','))
        on = len(o)
        f = set([el.partition('_')[2] if el.find('_') != -1 else el for el in line[5].split(',')])
        fn = len(f)
        if on != fn and fn < options.minimum:
            fragments.append([line[2],line[3]])


    file(options.fragments,'w').writelines(['\t'.join(line)+'\n' for line in fragments])

    #
