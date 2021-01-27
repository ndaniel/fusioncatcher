#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the symbols of all genes from the Ensembl database.



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
#using biomart in Python
import sys
import os
import optparse



if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the symbols of all genes from the Ensembl database."""
    version = "%prog 0.20 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the genes positions are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the genes positions are stored. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #
    symbols_filename = os.path.join(options.output_directory,'genes_symbols.txt')
    gtf_filename = os.path.join(options.output_directory,'organism.gtf')

    # print keeping only the chromosomes (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT)
    #chromosomes = [str(i) for i in range(1,1000)] + ['X','Y','MT','UN']
    #chromosomes = set(el.upper() for el in chromosomes)


    print "Starting..."

    g = [line.rstrip('\r\n').split("\t") for line in file(gtf_filename,'r').readlines() if (not line.startswith("#")) and line.rstrip('\r\n')]
    data = set()
    for line in g:
        if line[2] != 'gene':
            continue
        d = line[8].split(';')
        gene_id = ''
        gene_name = ''
        for e in d:
            e = e.strip()
            if e.startswith('gene_id '):
                gene_id = e[:-1].partition(' "')[2].upper()
            elif e.startswith('gene_name '):
                gene_name = e[:-1].partition(' "')[2].upper()
        if options.organism.lower() == "saccharomyces_cerevisiae":
            x = options.organism.upper().split('_')
            templ = "ENS"+x[0][0]+x[1][0:2]+"G"
            data.add((templ+gene_id.upper(),gene_name))
        elif gene_id.startswith('ENS'):
            data.add((gene_id,gene_name))

    # keep only the genes which start with ENS
    data = ['\t'.join(line)+'\n' for line in data]
    data = sorted(list(set(data)))
    file(symbols_filename,'w').writelines(data)


    print "End."
    #
