#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the list of candidate fusion genes. This list is hard coded
in here and it is manually curated from:


Greger et al. Tandem RNA Chimeras Contribute to Transcriptome Diversity in 
Human Population and Are Associated with Intronic Genetic Variants, 
Plos One, Aug 2014, http://dx.doi.org/10.1371/journal.pone.0104567 



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

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
    description = """It generates the list of pre-candidate fusion genes from 1000 genomes project."""
    version = "%prog 0.12 beta"

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

    parser.add_option("--skip-filter-overlap",
                      action="store_true",
                      dest="skip_filter_overlap",
                      default = False,
                      help="""If set then it filters out the known fusion genes where the (i) genes are fully overlapping, or (ii) the genes are partially overlapping and are on the same strand. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    print "Generating the list of 1000 genomes fusion genes..."
    fusions = dict()

    # manual curation from papers

    fusions['rattus_norvegicus'] = []

    fusions['mus_musculus'] = []

    fusions['canis_familiaris'] = []

    fusions['homo_sapiens'] = [
        ['ACTB','POTEE'],
        ['ACTB','POTEM'],
        ['AIG1','PARL'],
        ['AP5S1','MAVS'],
        ['ARHGAP19','SLIT1'],
        ['ARL4A','MTHFD1L'],
        ['BPTF','LRRC37A3'],
        ['C11ORF48','INTS5'],
        ['C1ORF189','TOX'],
        ['C2ORF27A','NBEA'],
        ['C6ORF72','PPIL4'],
        ['C7ORF55','LUC7L2'],
        ['CCL22','CX3CL1'],
        ['CENPE','BDH2'],
        ['CHURC1','FNTB'],
        ['CLN6','CALML4'],
        ['CNPY2','CS'],
        ['COPE','CERS1'],
        ['COPE','LASS1'],
        ['CORO7','PAM16'],
        ['COX5A','EDC3'],
        ['CTBS','GNG5'],
        ['CTSC','RAB38'],
        ['DDX5','POLG2'],
        ['EDARADD','ENO1'],
        ['EEF1A1','XPOT'],
        ['ELAVL1','TIMM44'],
        ['ENTPD1','CC2D2B'],
        ['FAM18B2','CDRT4'],
        ['FARSB','TRIM61'],
        ['FKBP1A','SDCBP2'],
        ['GNG5','CTBS'],
        ['GPI','PDCD2L'],
        ['HACL1','COLQ'],
        ['HAUS4','PRMT5'],
        ['HILPDA','EFCAB3'],
        ['HMSD','SERPINB8'],
        ['HSPE1','MOBKL3'],
        ['IFNAR2','IL10RB'],
        ['IFRD1','C7ORF53'],
        ['ISY1','RAB43'],
        ['JAK3','INSL3'],
        ['KIAA0101','CSNK1G1'],
        ['KIAA0494','ATPAF1'],
        ['LMAN2','MXD3'],
        ['LRRC33','PIGX'],
        ['LSP1','TNNT3'],
        ['MAPKAPK5','ACAD10'],
        ['MED8','ELOVL1'],
        ['METTL10','FAM53B'],
        ['METTL21B','TSFM'],
        ['NAIP','OCLN'],
        ['NDUFA13','YJEFN3'],
        ['NDUFB8','SEC31B'],
        ['NHP2L1','LLPH'],
        ['NRXN1','EIF2AK2'],
        ['NSUN4','FAAH'],
        ['PEX26','TUBA8'],
        ['PFKFB4','SHISA5'],
        ['PKHD1L1','EBAG9'],
        ['PLEKHO2','ANKDD1A'],
        ['POLA2','CDC42EP2'],
        ['POLR1A','REEP1'],
        ['PPIP5K1','CATSPER2'],
        ['PPRC1','NOLC1'],
        ['PRH1','PRR4'],
        ['PRIM1','NACA'],
        ['PRKAA1','TTC33'],
        ['PRKCB','YBX1'],
        ['PRR11','C17ORF71'],
        ['PRR13','PCBP2'],
        ['PXMP2','PGAM5'],
        ['RBM14','RBM4'],
        ['RHOQ','LRR1'],
        ['RNASET2','RPS6KA2'],
        ['RRM2','C2ORF48'],
        ['S1PR2','DNMT1'],
        ['SAV1','GYPE'],
        ['SDHAF2','C11ORF66'],
        ['SDHD','TEX12'],
        ['SLC35A3','HIAT1'],
        ['SLC39A1','CRTC2'],
        ['SLC43A3','PRG2'],
        ['SMC4','BCL6'],
        ['SMG1','ARL6IP1'],
        ['SNTB2','VPS4A'],
        ['SP100','HMGB1'],
        ['SUMO2','HN1'],
        ['SYNJ2BP','COX16'],
        ['TAGLN2','CCDC19'],
        ['TAP2','HLA-DOB'],
        ['TFG','GPR128'],
        ['TMBIM4','LLPH'],
        ['TNFAIP8L2','SCNM1'],
        ['TOMM5','FBXO10'],
        ['TOPORS','DDX58'],
        ['TPD52L2','DNAJC5'],
        ['TRIP12','SLC16A14'],
        ['TSC22D4','C7ORF61'],
        ['TSTD1','F11R'],
        ['TYK2','CDC37'],
        ['UBA2','WTIP'],
        ['UBE2J1','GABRR2'],
        ['UBE2J2','FAM132A'],
        ['UCHL3','LMO7'],
        ['UQCRQ','LEAP2'],
        ['VBP1','BRCC3'],
        ['VKORC1','PRSS53'],
        ['YARS2','NAP1L1'],
        ['ZNF175','CTU1'],
        ['ZNF343','SNRPB'],
        ['ZNF562','RBAK']

]



    data = fusions.get(options.organism.lower(),[])
    if data:

        #file_symbols = os.path.join(options.output_directory,'genes_symbols.txt')
        file_symbols = os.path.join(options.output_directory,'synonyms.txt')
        loci = symbols.generate_loci(file_symbols)

        genes = symbols.read_genes_symbols(file_symbols)

        d = []
        for (g1,g2) in data:
            if g1.upper() != g2.upper():
                ens1 = symbols.ensembl(g1.upper(),genes,loci)
                ens2 = symbols.ensembl(g2.upper(),genes,loci)
                if ens1 and ens2:
                    for e1 in ens1:
                        for e2 in ens2:
                            if e1 != e2:
                                d.append([e1,e2])

        data = ['\t'.join(sorted(line)) + '\n' for line in d]
        data = list(set(data))

        print "%d known fusion genes found in manually currated database" % (len(data),)

        if not options.skip_filter_overlap:
            d1 = []
            overlappings = ['ensembl_fully_overlapping_genes.txt',
                            'ensembl_same_strand_overlapping_genes.txt',
#                            'refseq_fully_overlapping_genes.txt',
#                            'refseq_same_strand_overlapping_genes.txt',
#                            'ucsc_fully_overlapping_genes.txt',
#                            'ucsc_same_strand_overlapping_genes.txt',
#                            'pairs_pseudogenes.txt',
#                            'paralogs.txt'
                            ]
            for ov in overlappings:
                p = os.path.join(options.output_directory,ov)
                print "Parsing file:",p
                if os.path.isfile(p):
                    d2 = sorted(set([tuple(sorted(line.rstrip('\r\n').split('\t'))) for line in file(p,'r').readlines() if line.rstrip('\r\n')]))
                d1.extend(d2)
            d = set()
            for line in d1:
                (a,b) = (line[0],line[1])
                if a > b:
                    (a,b) = (b,a)
                d.add("%s\t%s\n" % (a,b))
            skipped = [line for line in data if line in d]
            data = [line for line in data if not line in d]
            file(os.path.join(options.output_directory,'1000genomes_known_but_overlapping.txt'),'w').writelines(sorted(skipped))

            print "%d known fusion genes left after removing the overlappings" % (len(data),)

    file(os.path.join(options.output_directory,'1000genomes.txt'),'w').writelines(sorted(data))
    #
