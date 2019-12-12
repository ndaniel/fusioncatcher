#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It converts the FusionCatcher output into VCF format.


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


import os
import sys
import optparse
import itertools
import gc
import gzip
import datetime


if __name__=='__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It converts the FusionCatcher output into VCF format."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """Input file containing the candidate fusion genes in FusionCatcher's format, e.g. 'final-list_candidate-fusion-genes.vcf'.""")

    parser.add_option("--info","-n",
                      action = "store",
                      type = "string",
                      dest = "input_info_filename",
                      help = """Input file information about the FusionCatcher, e.g. 'info.txt'.""")


    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output VCF file.""")



    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")

    #print "Starting..."

    fi = None
    if options.input_filename == '-':
        fi = sys.stdin
    else:
        fi = open(options.input_filename,'r')
        
    fo = None
    if options.output_filename == '-':
        fo = sys.stdout
    else:
        fo = open(options.output_filename,'w')

    fc = [e.rstrip("\r\n").split("\t") for e in fi.readlines() if e.rstrip("\r\n")]
    header = dict([(e,i) for i,e in enumerate(fc.pop(0))])
    h_gs1 = header["Gene_1_symbol(5end_fusion_partner)"]
    h_gs2 = header["Gene_2_symbol(3end_fusion_partner)"]
    h_desc = header["Fusion_description"]
    #"Counts_of_common_mapping_reads"
    h_pairs = header["Spanning_pairs"]
    h_reads = header["Spanning_unique_reads"]
    #"Longest_anchor_found"
    #h_method = header["Fusion_finding_method"]
    h_pos1 = header["Fusion_point_for_gene_1(5end_fusion_partner)"]
    h_pos2 = header["Fusion_point_for_gene_2(3end_fusion_partner)"]
    h_ge1 = header["Gene_1_id(5end_fusion_partner)"]
    h_ge2 = header["Gene_2_id(3end_fusion_partner)"]
    #"Exon_1_id(5end_fusion_partner)"
    #"Exon_2_id(3end_fusion_partner)"
    h_seq = header["Fusion_sequence"]
    h_eff = header["Predicted_effect"]

    #
    #
    #
    s_fc_version = ""
    s_genome_fasta = ""
    s_data_directory = ""
    s_ensembl_version = ""
    s_files = ""
    s_command = ""
    # read the info file
    if options.input_info_filename:
        info = [e.rstrip("\r\n") for e in file(options.input_info_filename,"r")]
        for e in info:
            if e.startswith("Software version:"):
                s_fc_version = e.partition(": ")[2].replace(".py ","_v")
                break
        # read genome FASTA
        i = -1
        for e in info:
            i = i + 1
            if e.startswith("Genome FASTA files:"):
                s_genome_fasta = "ftp://"+info[i+2].replace("//","/")
        # read data directory
        for e in info:
            if e.startswith("data_directory ="):
                s_data_directory = e.partition(" = ")[2]
                break
        # read data directory
        for e in info:
            if e.startswith("Ensembl database version:"):
                s_ensembl_version = "ensembl_v"+e.partition(": ")[2]
                break
        # read input files
        flag = False
        rx = []
        for e in info:
            if e.startswith("Input files"):
                flag = True
                continue
            if flag:
                if e.startswith("---"):
                    continue
                if e:
                    rx.append(e)
                else:
                    break
        s_files = ",".join(rx)
        # read command line
        flag = False
        rx = []
        for e in info:
            if e.startswith("Command line used for launching FusionCatcher:"):
                flag = True
                continue
            if flag:
                if e.startswith("---"):
                    continue
                if e:
                    rx.append(e)
                else:
                    break
        s_command = " ".join([e.rstrip(" \\") for e in rx])


    #
    #
    #
    s_date = str(datetime.date.today())
    s_date = s_date.replace("-","")
    vcf = []
    vcf.append("##fileformat=VCFv4.3\n")
    vcf.append("##fileDate=%s\n" % (s_date,)) # ##fileDate=20090805
    vcf.append("##source=%s\n" % (s_fc_version)) # ##source=myImputationProgramV3.1
    vcf.append("##reference=%s\n" % (s_genome_fasta,)) # reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    vcf.append("##assembly=%s\n" % (s_genome_fasta,))
    vcf.append("##phasing=none\n")
    vcf.append("##geneAnnotation=%s\n" % (s_ensembl_version,))
    vcf.append("##inputFiles=%s\n" % (s_files,))
    vcf.append("##commandLine=%s\n" % (s_command,))
    vcf.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    vcf.append('##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of the mate breakend">\n')
    vcf.append('##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">\n')
    vcf.append('##INFO=<ID=SP,Number=1,Type=Integer,Description="Spanning pairs supporting the fusion gene">\n')
    vcf.append('##INFO=<ID=SR,Number=1,Type=Integer,Description="Spanning reads mapping on fusion junction">\n')
    vcf.append('##INFO=<ID=FUSION,Number=1,Type=String,Description="Fusion gene">\n')
    vcf.append('##INFO=<ID=END,Number=1,Type=String,Description="End of the fusion to which it belongs, that is 5 or 3 end">\n')
    vcf.append('##INFO=<ID=GENE5,Number=1,Type=String,Description="Gene symbol that is at the 5 prime end of the fusion">\n')
    vcf.append('##INFO=<ID=GENE3,Number=1,Type=String,Description="Gene symbol that is at the 3 prime end of the fusion">\n')
    vcf.append('##INFO=<ID=GENS5,Number=1,Type=String,Description="Ensembl Gene Id  that is at the 5 prime end of the fusion">\n')
    vcf.append('##INFO=<ID=GENS3,Number=1,Type=String,Description="Ensembl Gene Id that is at the 3 prime end of the fusion">\n')
    vcf.append('##INFO=<ID=EFFECT,Number=1,Type=String,Description="Predicted effect of the fusion">\n')
    vcf.append('##INFO=<ID=DESC,Number=1,Type=String,Description="Description of the fusion">\n')
    vcf.append('##INFO=<ID=PUBLISHED,Number=0,Type=Flag,Description="The fusion gene has been mentioned already in published scientific articles or publicly available scientific databases">\n')
    vcf.append('##INFO=<ID=STRAND,Number=1,Type=String,Description="The strand on which the gene, taking part in the fusion, is.">\n')
    vcf.append("#CHROM POS ID REF ALT QUAL FILTER INFO\n".replace(" ","\t"))

    # find unique fusion genes
    ggu = set()
    dup = set()
    for e in fc:
        gs1 = e[h_gs1]
        gs2 = e[h_gs2]
        xid = gs1+"--"+gs2
        if xid in ggu:
            ggu.remove(xid)
            dup.add(xid)
        elif xid not in dup:
            ggu.add(xid)
    print 'unique:',ggu
    print 'duplicates:',dup

    i = -1 
    for e in fc:
        i = i + 1
        pos1 = e[h_pos1].split(":")
        pos2 = e[h_pos2].split(":")
        gs1 = e[h_gs1]
        gs2 = e[h_gs2]
        ge1 = e[h_ge1]
        ge2 = e[h_ge2]
        sp = e[h_pairs]
        sr = e[h_reads]
        eff= e[h_eff]
        desc = e[h_desc]
        igh = True if e[h_seq].find("*NNNNNN") != -1 else False

        xid = gs1+"--"+gs2
        unique = True if xid in ggu else False
        print xid,unique
        
        xid5 = xid+"__5"
        xid3 = xid+"__3"
        xidu = xid
        if not unique:
            xid5 = xid5+"__"+str(i)
            xid3 = xid3+"__"+str(i)
            xidu = xidu+"__"+str(i)
        
        alt = "N["+pos2[0]+":"+pos2[1]+"["
        if pos2[2] == "-":
            alt = "N]"+pos2[0]+":"+pos2[1]+"]"
        if igh:
            alt = 'N'*20 + alt
        info = ['SVTYPE=BND',
                'MATEID='+xid3,
                'EVENT='+xidu,
                'PUBLISHED' if desc.find('known') !=-1 else '',
                'FUSION='+xid,
                'END=5',
                'STRAND='+pos1[2],
                'SR='+sr,
                'SP='+sp,
                'GENE5='+gs1,
                'GENS5='+ge1,
                'GENE3='+gs2,
                'GENS3='+ge2,
                'EFFECT='+eff,
                'DESC='+desc
                ]
        info = ';'.join([u for u in info if u])
        first = [pos1[0],pos1[1],xid5,'N',alt,'.','PASS',info]

        alt = "["+pos1[0]+":"+pos1[1]+"[N"
        if pos1[2] == "-":
            alt = "]"+pos1[0]+":"+pos1[1]+"]N"
        if igh:
            alt = alt + 'N'*20
        info = ['SVTYPE=BND',
                'MATEID='+xid5,
                'EVENT='+xidu,
                'PUBLISHED' if desc.find('known') !=-1 else '',
                'FUSION='+xid,
                'END=3',
                'STRAND='+pos2[2],
                'SR='+sr,
                'SP='+sp,
                'GENE5='+gs1,
                'GENS5='+ge1,
                'GENE3='+gs2,
                'GENS3='+ge2,
                'EFFECT='+eff,
                'DESC='+desc
                ]
        info = ';'.join([u for u in info if u])
        second = [pos2[0],pos2[1],xid3,'N',alt,'.','PASS',info]
        
        vcf.append('\t'.join(first)+"\n")
        vcf.append('\t'.join(second)+"\n")
    #
    #
    #
    fo.writelines(vcf)
    if options.output_filename != '-':
        fo.close()
    if options.input_filename == '-':
        fi.close()









