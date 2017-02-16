#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the list of candidate fusion genes. This list is hard coded
in here and it is manually curated from:


Bao Z.S. et al., RNA-seq of 272 gliomas revealed a novel, recurrent PTPRZ1-MET 
fusion transcript in secondary glioblastomas, Genome Research, 2014, 
http://dx.doi.org/10.1101/gr.165126.113




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
    description = """It generates the list of pre-candidate fusion genes from glioblastomas."""
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

    print "Generating the list of glioblastomas fusion genes..."
    fusions = dict()

    # manual curation from papers

    fusions['rattus_norvegicus'] = []

    fusions['mus_musculus'] = []

    fusions['canis_familiaris'] = []

    fusions['homo_sapiens'] = [
        ['A2M','NOL10'],
        ['ACER3','PAK1'],
        ['ADNP2','ITGB3BP'],
        ['AGAP2','TSFM'],
        ['AHCYL2','TMEM178B'],
        ['ANO4','GALNT9'],
        ['AP2A2','SBF2'],
        ['APEH','SHISA5'],
        ['AVIL','VSNL1'],
        ['B4GALNT3','DENND5B'],
        ['BAZ2A','LRP1'],
        ['BAZ2A','METTL1'],
        ['BCAT1','FIP1L1'],
        ['BCR','LZTR1'],
        ['BIN2','CELA1'],
        ['BRPF3','CLPSL1'],
        ['C12ORF5','KCNC2'],
        ['C15ORF57','CBX3'],
        ['C1ORF21','MIXL1'],
        ['C2CD3','SMAD9'],
        ['C2ORF48','PRMT8'],
        ['CBL','FBXO2'],
        ['CCDC136','FKBP9L'],
        ['CCM2','OGDH'],
        ['CD81','SPAG6'],
        ['CDK17','KCNC2'],
        ['CDK4','TSFM'],
        ['CHFR','SMARCD1'],
        ['CHST3','PTCH2'],
        ['CLASRP','SYMPK'],
        ['CLSTN1','KAZN'],
        ['COG3','GPC6'],
        ['COL22A1','BMP7'],
        ['COL4A2','FAM155A'],
        ['COQ10A','NINJ2'],
        ['CPNE5','ANXA6'],
        ['CS','STAT6'],
        ['CTC1','TMEM117'],
        ['CTDSP2','SLC6A15'],
        ['CTDSP2','XRCC6BP1'],
        ['CTNND1','TMEM135'],
        ['CUL1','CNTNAP2'],
        ['CXADR','C21ORF37'],
        ['DENND4C','ADAMTSL1'],
        ['DEPDC5','ALDH18A1'],
        ['DGKE','CACNA1G'],
        ['DHX35','LPAR1'],
        ['DNAJC3','RNASEH2B'],
        ['EDC3','COX5A'],
        ['EGFR','VSTM2A'],
        ['ELAVL2','PLGRKT'],
        ['ELK1','GPRC5A'],
        ['EML1','EVL'],
        ['ESAM','MSANTD2'],
        ['ETV6','NELL2'],
        ['FAM117A','LGALS9'],
        ['FAM155A','COL4A1'],
        ['FAM222A','MVK'],
        ['FAM222A','SIRT4'],
        ['FARP2','DTNB'],
        ['FGFR3','TACC3'],
        ['FOXP2','CAPZA2'],
        ['FOXP2','MET'],
        ['FRMD4A','PFKP'],
        ['GLTP','VPS29'],
        ['GOLPH3','APBB1IP'],
        ['GPR162','CCDC39'],
        ['GRIP1','RNFT2'],
        ['GRIP1','SOX5'],
        ['HUWE1','RIBC1'],
        ['IFT80','MLH1'],
        ['IFT81','ETNK1'],
        ['IL1RAP','FGF12'],
        ['IP6K2','USP19'],
        ['IQSEC1','GP1BA'],
        ['JAK1','HIVEP3'],
        ['KIAA0430','SPR'],
        ['KSR1','PKN3'],
        ['LAMA1','PTPRM'],
        ['LANCL2','DNAH11'],
        ['LHFP','SERP2'],
        ['LINC00662','AP1M2'],
        ['LMBR1L','GALNT9'],
        ['LNX2','LHFP'],
        ['LPCAT1','DENND3'],
        ['LRP1','BAZ2A'],
        ['LRRFIP2','GOLGA4'],
        ['MACF1','AHDC1'],
        ['MARCH6','MYO10'],
        ['MARCH9','SLCO1B7'],
        ['MAX','TTC7B'],
        ['MECOM','GAP43'],
        ['MED13L','GRIP1'],
        ['MED13L','YEATS4'],
        ['MGAT5','KIAA0825'],
        ['MID1','ARHGAP6'],
        ['MLL3','CHGB'],
        ['MRPS28','ASPH'],
        ['MSANTD2','ESAM'],
        ['MTAP','C9ORF92'],
        ['MYO6','SENP6'],
        ['NAB2','TAC3'],
        ['NEK6','RXRA'],
        ['NFATC3','CPNE2'],
        ['NMU','LNX1'],
        ['NPAT','RPS11'],
        ['NUS1','SLC25A13'],
        ['PAK4','AKT2'],
        ['PAK4','NFKBIB'],
        ['PER3','CAMTA1'],
        ['PITPNA','RAD51C'],
        ['PITPNM2','ACSS3'],
        ['PLA2G6','CRYBB1'],
        ['PLAGL2','HCK'],
        ['PLGRKT','HSPBP1'],
        ['PODXL','BGN'],
        ['POLN','SH3BP2'],
        ['PPP1R12A','C12ORF42'],
        ['PPP1R12A','OS9'],
        ['PPP6R1','MLLT3'],
        ['PRMT8','LAPTM4A'],
        ['PTEN','COL17A1'],
        ['PTN','DGKI'],
        ['PTPN1','CHMP4B'],
        ['PTPRZ1','MET'],
        ['PVRL2','SNTG1'],
        ['PXDC1','NFS1'],
        ['RAB3C','CNTN5'],
        ['RASSF7','SCUBE2'],
        ['RBMS2','LRP1'],
        ['RDH16','R3HDM2'],
        ['RIMKLB','MICALL2'],
        ['RNF213','SLC26A11'],
        ['RNFT2','MYO1A'],
        ['RPAP3','CPNE8'],
        ['RRN3P2','ANK2'],
        ['RTN3','TRPC6'],
        ['SCHIP1','ACOT9'],
        ['SCMH1','LIMA1'],
        ['SCUBE2','ZNF195'],
        ['SEC61G','LANCL2'],
        ['SFMBT1','CCRL2'],
        ['SFSWAP','TCP11L2'],
        ['SH3GL1','CHAF1A'],
        ['SIGLEC8','MBOAT7'],
        ['SLC2A8','KPNA3'],
        ['SLC6A8','GABRA3'],
        ['SMARCC2','TSFM'],
        ['SORL1','CNTNAP2'],
        ['SP1','TMEM132D'],
        ['SREBF2','EFCAB6'],
        ['SRGAP3','TATDN2'],
        ['SSH2','RTN1'],
        ['SSPN','GRIN2B'],
        ['ST7','CTTNBP2'],
        ['STAT2','CD9'],
        ['STAU1','CABP5'],
        ['STK3','VPS13B'],
        ['TCF7L1','KIF1B'],
        ['TFDP2','LMBRD1'],
        ['TFRC','IL1RAP'],
        ['TGFB1','SAE1'],
        ['TMEM117','CTC1'],
        ['TOM1L2','IGSF9B'],
        ['TPM3','ADAR'],
        ['TRMT6','C11ORF34'],
        ['TSFM','R3HDM2'],
        ['TTC19','ADORA2B'],
        ['TTLL11','FIBCD1'],
        ['UBE3C','CNTNAP2'],
        ['URI1','SLC6A20'],
        ['VEZT','CORO1C'],
        ['VHL','BRK1'],
        ['VSNL1','AVIL'],
        ['ZBTB17','NEB'],
        ['ZC3H11A','REN'],
        ['ZMIZ1','MAT1A'],
        ['ZNF140','CUX2'],
        ['ZNF160','ZNF415'],
        ['ZNF195','SCUBE2'],
        ['ZNF564','FARSA'],
        ['ZNF564','TUT1'],
        ['ZNF787','GPR4'],
        ['ZNF804A','PHLPP1'],
        ['ZNF814','ZNF417']


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
            file(os.path.join(options.output_directory,'gliomas_known_but_overlapping.txt'),'w').writelines(sorted(skipped))

            print "%d known fusion genes left after removing the overlappings" % (len(data),)

    file(os.path.join(options.output_directory,'gliomas.txt'),'w').writelines(sorted(data))
    #
