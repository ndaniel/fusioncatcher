#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It joins a list of homologous gene (found using find_homolog_genes.py).



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
import sys
import os
import optparse
import gc


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It joins a list of homologous gene (found using find_homolog_genes.py) by adding also their corresponding counts."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_1",
                      action="store",
                      type="string",
                      dest="input_1_filename",
                      help="""The first input file containing homolog genes (generated using find_homolog_genes.py).""")

    parser.add_option("--input_2",
                      action="store",
                      type="string",
                      dest="input_2_filename",
                      help="""The second input file containing homolog genes (generated using find_homolog_genes.py).""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text tab-separated file containing the candidate homologous genes (the genes are sorted alphabetically on the each line).""")

    parser.add_option("--rest",
                      action="store",
                      type="string",
                      dest="output_rest_filename",
                      help="""The output text tab-separated file containing the candidate homologous genes (the genes are sorted alphabetically on the each line) which were found to have less reads than the threshold specified using --reads.""")


    parser.add_option("--all",
                      action="store",
                      type="string",
                      dest="output_all_filename",
                      help="""The output text tab-separated file containing all candidate homologous genes (the genes are sorted alphabetically on the each line) and the corresponding counts of ambiguous reads.""")

    parser.add_option("--reads",
                      action="store",
                      type="int",
                      dest="reads",
                      default = 1,
                      help="""The minimum number of reads which map simultaneously on two genes in order to be considered as homolog genes (after the joining is done). Default is %default.""")




    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_1_filename and
            options.input_2_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    gc.disable()
    homolog = []
    if os.path.isfile(options.input_1_filename) or os.path.islink(options.input_1_filename):
        homolog = [line.rstrip('\r\n').split('\t') for line in file(options.input_1_filename,'r') if line.rstrip('\r\n')]
    gc.enable()
    gc.disable()
    homolog = dict([("%s\t%s"%(line[0],line[1]),int(line[2])) for line in homolog])
    gc.enable()
    data2 = []
    if os.path.isfile(options.input_2_filename) or os.path.islink(options.input_2_filename):
        data2 = file(options.input_2_filename,'r').readlines()
    for line in data2:
        li = line.rstrip('\r\n').split('\t')
        if li:
            k = "%s\t%s"% (li[0],li[1])
            v = int(li[2])
            homolog[k] = homolog.get(k,0) + v

    reads = options.reads
    if options.output_rest_filename:
        print "Writing...",options.output_rest_filename
        h = [k+'\t'+str(v)+'\n' for (k,v) in homolog.items() if v < reads]
        file(options.output_rest_filename,'w').writelines(h)

    if options.output_all_filename:
        print "Writing...",options.output_all_filename
        h = [k+'\t'+str(v)+'\n' for (k,v) in homolog.items()]
        file(options.output_all_filename,'w').writelines(h)

    print "Writing...",options.output_filename
    homolog = [k+'\t'+str(v)+'\n' for (k,v) in homolog.items() if v >= reads]
    file(options.output_filename,'w').writelines(homolog)

    print "The end."
