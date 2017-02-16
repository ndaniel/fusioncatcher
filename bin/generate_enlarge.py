#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the list of genes which need to be enlarge and covered fully.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2017 Daniel Nicorici

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
    description = """It generates the list of genes which need to be enlarge and covered fully."""
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


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    print "Generating the list of allowed/known fusion genes..."
    enlarge = dict()

    # manual curation from papers

    enlarge['rattus_norvegicus'] = []

    enlarge['mus_musculus'] = []

    enlarge['canis_familiaris'] = []

    enlarge['homo_sapiens'] = [
        'ADRBK1',
        'APRIL',
        'BACH2',
        'BANK1',
        'BCL1',
        'BCL10',
        'BCL11A',
        'BCL2',
        'BCL3',
        'BCL5',
        'BCL6',
        'BCL8',
        'BCL9',
        'BMI1',
        'BRD4',
        'C14ORF49',
        'CBFA2T3',
        'CCND1',
        'CCND2',
        'CCND3',
        'CCNE1',
        'CD44',
        'CDK6',
        'CEBPA',
        'CEBPB',
        'CEBPD',
        'CEBPE',
        'CEBPG',
        'CHST11',
        'CNN3',
        'CRLF2',
        'DDX6',
        'DEGS2',
        'DUSP22',
        'EBF1',
        'EGFR',
        'EPOR',
        'ERBB2',
        'ERVW-1',
        'ERVWE1',
        'ETV6',
        'FCGR2B',
        'FCRL4',
        'FOXP1',
        'GPR34',
        'ID4',
        'IGF2BP1',
        'IL31RA',
        'IRF4',
        'IRF8',
        'IRTA1',
        'LAPTM5',
        'LHX4',
        'MAF',
        'MAFB',
        'MIR125B1',
        'MMSET',
        'MRPL21',
        'MUC1',
        'NFKB2',
        'PAFAH1B2',
        'PAX5',
        'PCSK7',
        'RHOH',
        'SPIB',
        'TNFSF13',
        'WWOX',
        'ZFP36L1',
#        'NPM1', # already enlarged
        'CIC',
        'DUX4'

]



    data = enlarge.get(options.organism.lower(),[])
    if data:

        #file_symbols = os.path.join(options.output_directory,'genes_symbols.txt')
        file_symbols = os.path.join(options.output_directory,'synonyms.txt')
        loci = symbols.generate_loci(file_symbols)

        genes = symbols.read_genes_symbols(file_symbols)

        d = []
        for g in data:
            ens = symbols.ensembl(g.upper(),genes,loci)
            if ens:
                for e in ens:
                    d.append(e)
            else:
                print "   - Original:",g

        data = [line + '\n' for line in d]
        data = sorted(set(data))

        print "%d genes to be enlarged and covered" % (len(data),)


    file(os.path.join(options.output_directory,'enlarge.txt'),'w').writelines(data)
    #
