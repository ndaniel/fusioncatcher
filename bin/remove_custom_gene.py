#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Remove human genes which are wrongly annotated.



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
import sys
import os
import optparse
import symbols

def remove(outdir,
           transcript_id = ''):

    headex = dict([(line.rstrip('\r\n'),i) for i,line in enumerate(file(os.path.join(outdir,'exons_header.txt'),'r').readlines()) if line.rstrip('\r\n')])
    exons = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(outdir,'exons.txt'),'r').readlines() if line.rstrip('\r\n')]
    n = len(exons)
    egid = headex['ensembl_gene_id']
    exons = [line for line in exons if line[egid] != gene_id]
    m = len(exons)
    if n == m:
        print >>sys.stderr, "WARNING: %s gene not found in exons.txt!" % (gene_id,)
    else:
        file(os.path.join(outdir,'exons.txt'),'w').writelines(['\t'.join(line)+'\n' for line in exons])

    headge = dict([(line.rstrip('\r\n'),i) for i,line in enumerate(file(os.path.join(outdir,'genes_header.txt'),'r').readlines()) if line.rstrip('\r\n')])
    genes = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(outdir,'genes.txt'),'r').readlines() if line.rstrip('\r\n')]
    n = len(genes)
    ggid = headge['ensembl_gene_id']
    genes = [line for line in genes if line[ggid] != gene_id]
    m = len(genes)
    if n == m:
        print >>sys.stderr, "WARNING: %s gene not found in genes.txt!" % (gene_id,)
    else:
        file(os.path.join(outdir,'genes.txt'),'w').writelines(['\t'.join(line)+'\n' for line in genes])

    # remove it from the GTF file also
    gtf = [line for line in file(os.path.join(outdir,'organism.gtf'),'r').readlines()]
    n = len(gtf)
    gtf = [line for line in gtf if line.find(gene_id) == -1]
    m = len(gtf)
    if n == m:
        print >>sys.stderr, "WARNING: %s gene not found in organism.gtf!" % (gene_id,)
    else:
        file(os.path.join(outdir,'organism.gtf'),'w').writelines(gtf)




if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Remove human genes which are wrongly annotated."""
    version = "%prog 0.13 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the list of allowed candidate fusion genes is generated, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the list of allowed candidate fusion genes is generated. Default is '%default'.""")

    parser.add_option("--genes_ensembl",
                      action="store",
                      type="string",
                      dest="genes_ensembl",
                      help="""List of genes which should be removed from Ensembl.""")

    parser.add_option("--genes_refseq",
                      action="store_true",
                      dest="genes_refseq",
                      help="""Built-in list of fusion genes which should be removed from RefSeq.""")

    parser.add_option("--genes_ucsc",
                      action="store_true",
                      dest="genes_ucsc",
                      help="""Build-in of fusion genes which should be removed from UCSC.""")

    parser.add_option("--genes_gencode",
                      action="store_true",
                      dest="genes_gencode",
                      help="""Build-in of fusion genes which should be removed from Gencode.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    print "Remove genes which are wrongly annotated..."

    if options.genes_ensembl:
        gl = [line.rstrip("\r\n") for line in file(options.genes_ensembl,'r').readlines() if line.rstrip("\r\n")]
        for g in gl:
            remove(options.output_directory,g)

    if options.organism.lower() == 'homo_sapiens':
        #remove(options.output_directory,"ENST00000467125") # ENST00000467125 is listed as GOPC, but is actually the GOPC/ROS1
        #remove(options.output_directory,"ENST00000507166") # ENST00000507166 is listed as FIP1L1 but is really the FIP1L1/PDGFRA
        pass

    #
