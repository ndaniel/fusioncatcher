#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the list of known tumor genes (that are labeled as proto-oncogene or tumor 
suppresor genes) from UniProt database. This list is hard coded
in here and it is manually curated from: http://www.uniprot.org


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
    description = """It generates the list of known tumor genes from UniProt database."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the list of genes is generated, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the list of genes is generated. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    print "Generating the list of UNIPROT tumor genes..."
    fusions = dict()

    # manual curation from papers

    fusions['rattus_norvegicus'] = []

    fusions['mus_musculus'] = []

    fusions['canis_familiaris'] = []

    fusions['homo_sapiens'] = [
'ABL1',
 'AFDN',
 'AFF1',
 'AFF4',
 'AGAP2',
 'AGR2',
 'AIM2',
 'AKAP13',
 'AKIP1',
 'AKT1',
 'AKT2',
 'ALK',
 'APC',
 'ARAF',
 'ARHGAP20',
 'ARHGAP26',
 'ARHGEF12',
 'ARHGEF5',
 'ARID3B',
 'ASPSCR1',
 'ATM',
 'AURKA',
 'AXIN1',
 'AXL',
 'BANP',
 'BAX',
 'BCAS3',
 'BCAS4',
 'BCL10',
 'BCL2',
 'BCL3',
 'BCL6',
 'BCL9',
 'BCR',
 'BIN1',
 'BMI1',
 'BRAF',
 'BRCA1',
 'BRCA2',
 'BRCC3',
 'BRD7',
 'BRI3BP',
 'BRMS1',
 'BTG1',
 'BUB1B',
 'C10ORF90',
 'C10ORF99',
 'CADM1',
 'CADM4',
 'CARS',
 'CBFB',
 'CBL',
 'CCAR2',
 'CCDC6',
 'CCND1',
 'CCNL1',
 'CDC73',
 'CDK2AP1',
 'CDKN1B',
 'CDKN1C',
 'CDKN2A',
 'CDKN2B',
 'CDKN2D',
 'CDT1',
 'CHD5',
 'CHEK2',
 'CMC4',
 'CREB3L2',
 'CREBL2',
 'CRK',
 'CSF1R',
 'CSNK2A3',
 'CTCF',
 'CYLD',
 'DAB2',
 'DAB2IP',
 'DAPK3',
 'DCC',
 'DCUN1D1',
 'DDIT3',
 'DDX6',
 'DEC1',
 'DEK',
 'DFNA5',
 'DIS3L2',
 'DLC1',
 'DLEC1',
 'DLEU1',
 'DMBT1',
 'DMTF1',
 'DMTN',
 'DOCK4',
 'DPH1',
 'EFNA1',
 'EGFR',
 'ELAC2',
 'ELF4',
 'ELL',
 'ENTPD5',
 'EPB41L3',
 'EPHB2',
 'EPS15',
 'ERG',
 'ERRFI1',
 'ETS1',
 'ETS2',
 'ETV1',
 'ETV6',
 'EVI2A',
 'EVI2B',
 'EWSR1',
 'EXT1',
 'EXT2',
 'FAM120A',
 'FAM83A',
 'FAM83B',
 'FAM83D',
 'FCGR2B',
 'FER',
 'FES',
 'FGF3',
 'FGF4',
 'FGF5',
 'FGF6',
 'FGFR2',
 'FGR',
 'FH',
 'FHIT',
 'FLCN',
 'FLI1',
 'FLT3',
 'FOS',
 'FOXO1',
 'FOXO3',
 'FOXO4',
 'FRAT1',
 'FRK',
 'FSTL3',
 'FUS',
 'FYN',
 'GAS7',
 'GFI1B',
 'GLI1',
 'GMPS',
 'GNAS',
 'GPR68',
 'GPRC5A',
 'HCK',
 'HIC1',
 'HIF3A',
 'HLF',
 'HMGA2',
 'HOPX',
 'HOXA9',
 'HPGD',
 'HRAS',
 'HTATIP2',
 'IL2',
 'ING1',
 'ING4',
 'IRF1',
 'JAK2',
 'JAZF1',
 'JUN',
 'KANK1',
 'KAT6A',
 'KCTD11',
 'KDSR',
 'KIT',
 'KLK10',
 'KMT2A',
 'KRAS',
 'LATS1',
 'LATS2',
 'LCK',
 'LETMD1',
 'LGR6',
 'LHX4',
 'LIMD1',
 'LIN9',
 'LITAF',
 'LMO1',
 'LMO2',
 'LYL1',
 'LYN',
 'LZTS1',
 'MAF',
 'MAFA',
 'MAFB',
 'MAP3K8',
 'MAPKAPK5',
 'MAS1',
 'MCC',
 'MCF2',
 'MCF2L',
 'MCTS1',
 'MDM2',
 'MDS2',
 'MECOM',
 'MERTK',
 'MET',
 'MFHAS1',
 'MKL1',
 'MLH1',
 'MLLT1',
 'MLLT10',
 'MLLT11',
 'MLLT3',
 'MLLT6',
 'MN1',
 'MOS',
 'MSH2',
 'MTCP1',
 'MTSS1',
 'MTUS1',
 'MUC1',
 'MUTYH',
 'MXI1',
 'MYB',
 'MYC',
 'MYCN',
 'MYEOV',
 'MYH11',
 'NAT6',
 'NBL1',
 'NCKIPSD',
 'NCOA1',
 'NCOA4',
 'NDRG2',
 'NET1',
 'NEURL1',
 'NF1',
 'NF2',
 'NFKB2',
 'NKX3-1',
 'NPM1',
 'NPRL2',
 'NR4A3',
 'NRAS',
 'NSD1',
 'NSD2',
 'NSD3',
 'NTRK1',
 'NUP214',
 'OLIG2',
 'PAF1',
 'PALB2',
 'PANO1',
 'PARK7',
 'PATZ1',
 'PAX3',
 'PAX5',
 'PAX7',
 'PBRM1',
 'PBX1',
 'PCM1',
 'PDCD4',
 'PDGFB',
 'PDGFD',
 'PDGFRA',
 'PDGFRB',
 'PHB',
 'PHLDA3',
 'PHLPP1',
 'PHLPP2',
 'PICALM',
 'PIK3CA',
 'PIM1',
 'PIM2',
 'PIM3',
 'PINX1',
 'PLAG1',
 'PLEKHG2',
 'PLEKHO1',
 'PLK2',
 'PLPP5',
 'PML',
 'PMS1',
 'PMS2',
 'PNN',
 'POU2AF1',
 'PRCC',
 'PRKCA',
 'PRKCD',
 'PRKCDBP',
 'PRKCI',
 'PRR5',
 'PTCH1',
 'PTEN',
 'PTTG1',
 'PTTG2',
 'PTTG3P',
 'PYCARD',
 'PYHIN1',
 'RAB8A',
 'RAF1',
 'RAP1A',
 'RARA',
 'RARB',
 'RASA1',
 'RASL10A',
 'RASSF1',
 'RASSF2',
 'RASSF4',
 'RASSF5',
 'RB1',
 'RB1CC1',
 'RBL1',
 'RBL2',
 'RBM15',
 'RBMX',
 'RECK',
 'REL',
 'RET',
 'RHOA',
 'RHOB',
 'RNF213',
 'ROS1',
 'RPS6KA2',
 'RRAS2',
 'RUNX1',
 'RUNX1T1',
 'SASH1',
 'SDHA',
 'SEC31A',
 'SET',
 'SH3GL1',
 'SIK1',
 'SIRT4',
 'SKI',
 'SLC5A8',
 'SMARCB1',
 'SPECC1',
 'SPI1',
 'SRC',
 'SS18',
 'SSX1',
 'SSX2',
 'ST20',
 'STARD13',
 'STIL',
 'STK11',
 'STYK1',
 'SUFU',
 'SUSD2',
 'SUSD6',
 'SUZ12',
 'SYNPO2',
 'TAF15',
 'TAL1',
 'TAL2',
 'TBC1D3',
 'TBRG1',
 'TCF3',
 'TCHP',
 'TCL1A',
 'TCL1B',
 'TCP10L',
 'TCTA',
 'TET2',
 'TFE3',
 'TFG',
 'TFPT',
 'TLX1',
 'TMEM127',
 'TNFRSF17',
 'TOP1',
 'TP53',
 'TP53INP1',
 'TP73',
 'TPM3',
 'TPR',
 'TRIM24',
 'TRIM27',
 'TRIM37',
 'TSC1',
 'TSC2',
 'TUSC2',
 'TXNIP',
 'UFL1',
 'USP4',
 'USP6',
 'VAV1',
 'VHL',
 'VWA5A',
 'WDR11',
 'WISP1',
 'WNT1',
 'WNT3',
 'WT1',
 'WWOX',
 'WWTR1',
 'XAF1',
 'XRN1',
 'YAP1',
 'YES1',
 'ZBTB16',
 'ZBTB7C',
 'ZDHHC17',
 'ZFYVE19',
 'ZMYND11',
 'ZNF320',
 'ZNF521'


]

    file(os.path.join(options.output_directory,'tumor_genes.txt'),'w').write("")

    data = fusions.get(options.organism.lower(),[])
    if data:

        #file_symbols = os.path.join(options.output_directory,'genes_symbols.txt')
        file_symbols = os.path.join(options.output_directory,'synonyms.txt')
        loci = symbols.generate_loci(file_symbols)

        genes = symbols.read_genes_symbols(file_symbols)

        d = []
        for g in data:
            ens = symbols.ensembl(g.upper(),genes,loci)
            if ens:
                d.extend(ens)

        data = [line + '\n' for line in d]
        data = sorted(set(data))

        print "%d genes found in manually curated database" % (len(data),)

    file(os.path.join(options.output_directory,'tumor_genes.txt'),'w').writelines(sorted(data))
    #
