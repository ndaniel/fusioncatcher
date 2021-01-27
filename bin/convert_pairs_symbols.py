#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It converts a input list of gene symbols which form a pair (there are two gene
symbols per line separated by tab) into their Ensembl gene ids.



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
import sys
import os
import optparse
import symbols

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It converts a input list of gene symbols which form a pair (there are two gene
symbols per line separated by tab) into their Ensembl gene ids."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage, description = description, version = version)

    parser.add_option("--input",
                      action = "store",
                      type = "string",
                      dest = "input",
                      help="""The input file containing on each line two gene symbols separated by tab.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output",
                      help="""The output file which will contain the Ensembl gene ids, where are two genes ids per line separated by tab.""")

    parser.add_option("--filter",
                      action="store",
                      type="string",
                      dest="filter",
                      help="""Input file containing the Ensembl Ids of genes pairs which should be removed, where are two genes ids per line separated by tab.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input and options.output):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    print "Reading the input file...", options.input
    mygenes = [line.rstrip('\r\n').split('\t') for line in file(options.input,'r').readlines() if line.rstrip('\r\n')]

    data = []
    if mygenes:

        file_symbols1 = os.path.join(os.path.dirname(options.output),'genes_symbols.txt')
        file_symbols2 = os.path.join(os.path.dirname(options.output),'synonyms.txt')

        loci1 = symbols.generate_loci(file_symbols1)
        loci2 = symbols.generate_loci(file_symbols2)

        genes1 = symbols.read_genes_symbols(file_symbols1)
        genes2 = symbols.read_genes_symbols(file_symbols2)

        d = []
        for (g1,g2) in mygenes:
            if g1 and g2 and g1.upper() != g2.upper():
                ens1 = symbols.ensembl(g1.upper(),genes1,loci1)
                ens2 = symbols.ensembl(g2.upper(),genes1,loci1)
                if not ens1:
                    ens1 = symbols.ensembl(g1.upper(),genes2,loci2)
                if not ens2:
                    ens2 = symbols.ensembl(g2.upper(),genes2,loci2)
                    
                if ens1 and ens2:
                    for e1 in ens1:
                        for e2 in ens2:
                            if e1 and e2 and e1 != e2:
                                d.append([e1,e2])

        data = ['\t'.join(sorted(line)) + '\n' for line in d]
        data = list(set(data))
        data = sorted(data)
        if options.filter:
            removal = set([line for line in file(options.filter,'r').readlines()])
            data = [line for line in data if line not in removal]

        print "%d known genes found in the input file" % (len(data),)

    file(options.output,'w').writelines(data)
    #
