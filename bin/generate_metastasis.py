#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the list of candidate fusion genes. This list is hard coded
in here and it is manually curated from:


Priestley et al., Pan-cancer whole-genome analyses of metastatic solid tumours, Nature, 2019, 
https://doi.org/10.1038/s41586-019-1689-y




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
    description = """It generates the list of pre-candidate fusion genes from metastasis."""
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

    print "Generating the list of metastasis fusion genes..."
    fusions = dict()

    # manual curation from papers

    fusions['rattus_norvegicus'] = []

    fusions['mus_musculus'] = []

    fusions['canis_familiaris'] = []

    fusions['homo_sapiens'] = [
        ['ACTB','BRAF'],
        ['ADAM32','FGFR1'],
        ['ADAM9','NRG1'],
        ['AFF3','ERG'],
        ['AGAP3','BRAF'],
        ['AGK','BRAF'],
        ['AGMO','ETV1'],
        ['AGR2','ETV1'],
        ['ANKDD1A','ERG'],
        ['ARMC8','ERBB4'],
        ['ASF1B','MAST1'],
        ['ASH2L','FGFR1'],
        ['ATAD2','NTRK3'],
        ['ATRNL1','ERG'],
        ['B3GALT5','ERG'],
        ['BACE2','ERG'],
        ['BACH2','ERG'],
        ['BCOR','KMT2A'],
        ['BCOR','SH3KBP1'],
        ['BCOR','STMN2'],
        ['BCR','GNAZ'],
        ['BCR','OBSCN'],
        ['BRD4','HYDIN'],
        ['BRD4','SEMA6A'],
        ['C16ORF45','ALK'],
        ['CCDC6','RET'],
        ['CCDC94','MAST1'],
        ['CD74','NRG1'],
        ['CHCHD3','BRAF'],
        ['CIC','MEGF8'],
        ['CRTC1','MAML2'],
        ['CSMD1','NRG1'],
        ['CSMD3','ABL1'],
        ['CTAGE5','BRAF'],
        ['CTNNB1','PLAG1'],
        ['CUL1','BRAF'],
        ['DCUN1D2','AKT3'],
        ['DENND2A','BRAF'],
        ['DIAPH2','ERBB4'],
        ['EGFR','C1ORF101'],
        ['EGFR','SEPT14'],
        ['EHF','ETV1'],
        ['EIF2AK2','ALK'],
        ['EIF3E','RSPO2'],
        ['EML4','ALK'],
        ['EML4','NTRK3'],
        ['EP300','EP300'],
        ['EPB41','BRAF'],
        ['ERBB4','ERBB4'],
        ['ESR1','ABCB1'],
        ['ESR1','AKAP12'],
        ['ESR1','ARNT2'],
        ['ESR1','CLINT1'],
        ['ESR1','ESRRG'],
        ['ESR1','GRIP1'],
        ['ESR1','LPP'],
        ['ESR1','NCOA1'],
        ['ESR1','PLEKHG1'],
        ['ESR1','SCNN1B'],
        ['ESR1','TCF12'],
        ['ESR1','TNRC6B'],
        ['ETV6','CRLS1'],
        ['EWSR1','ATF1'],
        ['EWSR1','FLI1'],
        ['EWSR1','SSX1'],
        ['EYA1','FGFR1'],
        ['EZR','ROS1'],
        ['FGFR1','VPS13B'],
        ['FGFR2','ATE1'],
        ['FGFR2','BICC1'],
        ['FGFR2','CREB5'],
        ['FGFR2','MCU'],
        ['FGFR2','ROCK1'],
        ['FGFR2','SPEG'],
        ['FGFR2','TACC2'],
        ['FGFR2','TBC1D4'],
        ['FGFR2','TENC1'],
        ['FGFR3','TACC3'],
        ['FIP1L1','PDGFRA'],
        ['FSTL5','NRG1'],
        ['FUS','DDIT3'],
        ['FUS','ERG'],
        ['FUS','TFCP2'],
        ['GATA4','FGFR1'],
        ['GOPC','ROS1'],
        ['HIP1','BRAF'],
        ['HMGA2','CLIC5'],
        ['HMGA2','TTC28'],
        ['HMGN2P46','MYC'],
        ['HNRNPA2B1','ETV1'],
        ['HNRNPH1','ETV4'],
        ['IFT88','ERBB4'],
        ['IKZF2','ERBB4'],
        ['JAK2','GNAQ'],
        ['JAK2','TTC28'],
        ['KCNJ15','ERG'],
        ['KMT2A','ATP5L'],
        ['KMT2A','BCOR'],
        ['KMT2A','CCDC84'],
        ['LAMB1','MET'],
        ['LAMB4','MET'],
        ['LHFPL2','ALK'],
        ['LIPC','ERBB4'],
        ['LSAMP','NRG1'],
        ['MKRN1','BRAF'],
        ['MPP4','EP300'],
        ['MXD1','ALK'],
        ['MYB','NFIB'],
        ['NAB2','STAT6'],
        ['NDRG1','ERG'],
        ['NRG1','FGFR1'],
        ['NRG1','NRG1'],
        ['NSD1','ZNF346'],
        ['NTRK2','NTRK2'],
        ['NTRK3','NTRK3'],
        ['PARP12','BRAF'],
        ['PAX3','FOXO1'],
        ['PAX3','SGPP2'],
        ['PLCH1','PIK3CB'],
        ['PLXND1','TMCC1'],
        ['PRTG','NRG1'],
        ['PSD3','NRG1'],
        ['PTPRK','RSPO3'],
        ['RABEP1','USP6'],
        ['RANBP2','ALK'],
        ['RAP1GAP2','ALK'],
        ['RUNX1','ERG'],
        ['SEPT8','AFF4'],
        ['SEZ6','ERG'],
        ['SHANK2','PLAG1'],
        ['SKAP2','MET'],
        ['SLC44A1','BRAF'],
        ['SLC45A3','CHST10'],
        ['SLC45A3','ERG'],
        ['SLC45A3','ETV4'],
        ['SLC45A3','FLI1'],
        ['SLC45A3','NUCKS1'],
        ['SND1','BRAF'],
        ['SPAG17','ALK'],
        ['SS18','KCTD1'],
        ['SS18','SSX1'],
        ['STAT3','ETV4'],
        ['SV2B','NTRK3'],
        ['TMEM40','RAF1'],
        ['TMPRSS2','BACE2'],
        ['TMPRSS2','ERG'],
        ['TMPRSS2','ETV4'],
        ['TMPRSS2','NDUFAF2'],
        ['TMPRSS2','PNPLA7'],
        ['TMPRSS2','RIPK4'],
        ['TPM3','ASH1L'],
        ['TPM3','KIAA1217'],
        ['TPR','NTRK1'],
        ['UBTF','ETV4'],
        ['UNC5D','NRG1'],
        ['VTI1A','TCF7L2'],
        ['WASF2','BRAF'],
        ['WBP2','ETV4'],
        ['YWHAE','CRK'],
        ['YWHAE','IQCK'],
        ['YWHAE','PITPNA'],
        ['YWHAE','SMG6'],
        ['ZCCHC24','MYC']

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
            file(os.path.join(options.output_directory,'metastasis_known_but_overlapping.txt'),'w').writelines(sorted(skipped))

            print "%d known fusion genes left after removing the overlappings" % (len(data),)

    file(os.path.join(options.output_directory,'metastasis.txt'),'w').writelines(sorted(data))
    #
