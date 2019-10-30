#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It finds a candidate list of fusion genes.



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

#
"""
# Columns' description of the MAP output file:
#    1.      Name of read that aligned
#    2.      Reference strand aligned to, + for forward strand, - for reverse
#    3.      Name of reference sequence where alignment occurs, or numeric ID if no name was provided
#    4.      0-based offset into the forward reference strand where leftmost character of the alignment occurs
#    5.      Read sequence (reverse-complemented if orientation is -). If the read was in colorspace, then the sequence shown in this column is the sequence of decoded nucleotides, not the original colors. See the Colorspace alignment section for details about decoding. To display colors instead, use the --col-cseq option.
#    6.      ASCII-encoded read qualities (reversed if orientation is -). The encoded quality values are on the Phred scale and the encoding is ASCII-offset by 33 (ASCII char !). If the read was in colorspace, then the qualities shown in this column are the decoded qualities, not the original qualities. See the Colorspace alignment section for details about decoding. To display colors instead, use the --col-cqual option.
#    7.      If -M was specified and the prescribed ceiling was exceeded for this read, this column contains the value of the ceiling, indicating that at least that many valid alignments were found in addition to the one reported. Otherwise, this column contains the number of other instances where the same sequence aligned against the same reference characters as were aligned against in the reported alignment. This is not the number of other places the read aligns with the same number of mismatches. The number in this column is generally not a good proxy for that number (e.g., the number in this column may be '0' while the number of other alignments with the same number of mismatches might be large).
#    8.      Comma-separated list of mismatch descriptors. If there are no mismatches in the alignment, this field is empty. A single descriptor has the format offset:reference-base>read-base. The offset is expressed as a 0-based offset from the high-quality (5') end of the read.


"""
import sys
import os
import optparse
import gc
import gzip

#########################
def line_from(a_map_filename):
    # it gives chunks from a_map_filename which is assumed to be ordered by the
    # name of reads (i.e. column 1 = read name)
    fin = None
    if a_map_filename.lower().endswith('.gz'):
        fin = gzip.open(a_map_filename,'r')
    else:
        fin = open(a_map_filename,'r')
    buffer_size = 10**8
    while True:
        gc.disable()
        lines=fin.readlines(buffer_size)
        gc.enable()
        if not lines:
            break
        gc.disable()
        lines=[line.rstrip('\r\n').split('\t',4)[0:4] for line in lines if line.rstrip('\r\n')]
        gc.enable()
        for line in lines:
            yield line
    fin.close()

#########################
def read_from(a_map_filename):
    last_r=''
    chunk=dict()
#    sequ = dict()
    for line in line_from(a_map_filename):
        rr=line[0].partition('/') # read name with /1 or /2
        r=rr[0] # read name without /1 or /2
        m=rr[2] # /1 or /2
        p = line[3] # position
        if not chunk:
            last_r=r
        if last_r!=r:
            # line[2] is column no 3 in the BOWTIE MAP file which contains
            # the "reference sequence name" = "transcript name"
            yield (last_r,chunk)
            last_r = r
            chunk = dict()
            seq = dict()
        #tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103
#        w = line[2].split(';')
#        g = ''
#        t = ''
#        for el in w:
#            if el.startswith('ge='):
#                g = el[3:]
#                if t:
#                    break
#            elif el.startswith('tr='):
#                t = el[3:]
#                if g:
#                    break
#        g=[el[3:] for el in w if el.startswith('ge=')]
#        t=[el[3:] for el in w if el.startswith('tr=')]
        w = line[2].partition(';')
        t = w[0]
        g = w[2]
        s=line[1] # strand aligned to
        #sequ[m] = line[4] # read sequence
        if chunk.has_key(m):
            chunk[m].append((g,t,s,p))
        else:
            chunk[m]=[(g,t,s,p)]
    if chunk:
        yield (last_r,chunk)

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It finds a candidate list of fusion genes."""
    version="%prog 0.11 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_map_filename",
                      help="""The input file in Bowtie MAP format (sorted by read name, i.e. column number 1) containing the short reads mapped on the transcripts. This should contain paired-end data. """)

    parser.add_option("--filter_gene_pairs",
                      action="store",
                      type="string",
                      dest="input_filter_gene_pairs_filename",
                      help="""The input text file tab separated format containg gene pairs which should be removed from the output fusion-genes list.""")

    parser.add_option("--filter_genes",
                      action="store",
                      type="string",
                      dest="input_filter_genes_filename",
                      help="""The input text file containing genes (one gene on each line) which should be removed from the output fusion-genes list.""")

    parser.add_option("--input_hugo",
                      action="store",
                      type="string",
                      dest="input_hugo_filename",
                      help="""The input database used for linking ENSEMBL GENE ID to HUGO gene names.""")

    parser.add_option("--output_fusion_genes",
                      action="store",
                      type="string",
                      dest="output_fusion_genes_filename",
                      help="""The output text tab-separated file containing the candidate fusion genes (the genes are sorted alphabetically on the each line).""")

    parser.add_option("--output_fusion_transcripts",
                      action="store",
                      type="string",
                      dest="output_fusion_transcripts_filename",
                      help="""The output text tab-separated file containing the candidate fusion transcripts (the genes are sorted alphabetically on the each line).""")

    parser.add_option("--output_fusion_reads",
                      action="store",
                      type="string",
                      dest="output_fusion_reads_filename",
                      help="""The output text tab-separated file containing the candidate fusion genes and transcripts together with the ids/names of supporting reads.""")

    parser.add_option("--output_fusion_reads_split",
                      action="store",
                      type="string",
                      dest="output_fusion_reads_split_filename",
                      help="""A file containing paths to candidate fusion genes and transcripts together with the ids/names of supporting reads.""")


    parser.add_option("--output_fusion_reads_simple",
                      action="store",
                      type="string",
                      dest="output_fusion_reads_simple_filename",
                      help="""The output text file containing one each line reads ids of the supporting reads.""")


    parser.add_option("--output_missing_mate_reads",
                      action="store",
                      type="string",
                      dest="output_missing_mate_reads_filename",
                      help="""The output text tab-separated file containing the reads which have their mate read not mapped together to the gene name on which they map.""")



    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_map_filename and
            options.input_hugo_filename and
            options.output_fusion_genes_filename
            ):
        parser.print_help()
        sys.exit(1)

    homologs=set()
    if options.input_filter_gene_pairs_filename:
        print "Reading...",options.input_filter_gene_pairs_filename
        homologs=set([line.rstrip('\r\n') for line in file(options.input_filter_gene_pairs_filename,'r') if line.rstrip('\r\n')])

    no_proteins=set()
    if options.input_filter_genes_filename:
        print "Reading...",options.input_filter_genes_filename
        no_proteins=set([line.rstrip('\r\n') for line in file(options.input_filter_genes_filename,'r') if line.rstrip('\r\n')])


    print "Reading...",options.input_hugo_filename
    # it is assumed that it is in this format
#    ensembl_gene_id
#    hgnc_symbol
#    description
#    description
    hugo=dict()
    for el in [line.rstrip('\r\n').split('\t')[0:2] for line in file(options.input_hugo_filename,'r').readlines() if line.rstrip('\r\n')]:
        g=el[0]
        h=el[1].replace(' ','_')
        if hugo.has_key(g):
            hugo[g]=hugo[g]+','+h
        else:
            hugo[g]=h

    flag_transcripts = False
    if options.output_fusion_transcripts_filename:
        flag_transcripts = True

    print "Processing and reading...",options.input_map_filename
    fusion_transcripts=dict()
    fusion_genes=dict()
    fusion_reads = dict()
    i_paired=0
    #
    missing_mates = ['missing_mate_read\tfound_mate_read\tfound_mate_read_maps_on_following_genes\n']
    #
    g_a = None
    t_a = None
    s_a = None
    g_b = None
    t_b = None
    s_b = None
    unique_positions = set()
    for (a_read,piece) in read_from(options.input_map_filename):

        #
        i_paired = i_paired+1
        if i_paired % 5000000 == 0:
            print '...',i_paired,'paired-end reads processed'
        flag = False
        if piece.has_key('1'):
            a=zip(*piece['1'])
            g_a=a[0]
            t_a=a[1]
            s_a=a[2]
            p_a = a[3]
        else:
            flag = True
        if piece.has_key('2'):
            b=zip(*piece['2'])
            g_b=b[0]
            t_b=b[1]
            s_b=b[2]
            p_b = b[3]
        else:
            flag = True

        if flag:
            if piece.has_key('1'):
                missing_mates.append("%s/2\t%s/1\t%s\n" % (a_read,a_read,','.join(set(g_a))))
            elif piece.has_key('2'):
                missing_mates.append("%s/1\t%s/2\t%s\n" % (a_read,a_read,','.join(set(g_b))))
            continue
        if set(g_a).intersection(set(g_b)): # if there are common genes skip the paired-end read
            continue

        unique_genes=set()
        for i in xrange(len(g_a)):
            for j in xrange(len(g_b)):
                g = '%s\t%s' % (g_a[i],g_b[j]) if g_a[i] < g_b[j] else '%s\t%s' % (g_b[j],g_a[i])
                z1 = "%s+++%s" % (t_a[i],p_a[i])
                z2 = "%s+++t%s" % (t_b[j],p_b[j])
                z = "%s-+-%s" % (z1,z2) if z1 < z2 else "%s-+-%s" % (z2,z1)
                #g = '\t'.join(sorted([g_a[i],g_b[j]]))
                if ( s_a[i] == s_b[j] or
                    (g in homologs) or
                    (g_a[i] in no_proteins) or
                    (g_b[j] in no_proteins)):
                    continue
                else:
                    if g not in unique_genes: # local unique gene fusions
                        if (z not in unique_positions):
                            fusion_genes[g] = fusion_genes.get(g,0)+1

                        if not fusion_reads.has_key(g):
                            fusion_reads[g] = set()
                        fusion_reads[g].add(a_read)

                        unique_genes.add(g)
                    unique_positions.add(z)

                    if flag_transcripts:

                        t1 = t_a[i]+'\t'+g_a[i]
                        t2 = t_b[j]+'\t'+g_b[j]
                        t = '%s\t%s' % (t1,t2) if t1 < t2 else '%s\t%s' % (t2,t1)
                        #t='\t'.join(sorted([t_a[i]+'\t'+g_a[i],t_b[j]+'\t'+g_b[j]]))
                        fusion_transcripts[t] = fusion_transcripts.get(t,0)+1
#                        if not fusion_reads.has_key(t):
#                            fusion_reads[g] = list()
#                       fusion_reads[t].append(a_read)


    print "Writing...",options.output_fusion_genes_filename
    fo=open(options.output_fusion_genes_filename,'w')
    data=['\t'.join(map(str,line)) for line in fusion_genes.items()]
    data=[line.split('\t') for line in data]
    data = sorted(data, key=lambda x:(int(x[2]),x[0],x[1]), reverse = True)
    data=[line+[hugo[line[0]],hugo[line[1]]] for line in data]
    data=['\t'.join(line)+'\n' for line in data]
    data.insert(0,'Fusion_gene_1\tFusion_gene_2\tCount_paired-end_reads\tFusion_gene_symbol_1\tFusion_gene_symbol_2\n')
    fo.writelines(data)
    fo.close()

    if flag_transcripts:
        print "Writing...",options.output_fusion_transcripts_filename
        fo=open(options.output_fusion_transcripts_filename,'w')
        data=['\t'.join(map(str,line))+'\n' for line in sorted(fusion_transcripts.items(), key=lambda x:(int(x[1]),x[0]), reverse = True)]
        data.insert(0,'Fusion_transcript_1\tFusion_gene_1\tFusion_transcript_2\tFusion_gene_2\tCount_paired-end_reads\n')
        fo.writelines(data)
        fo.close()


    print "Writing...",options.output_fusion_reads_filename
    fo=open(options.output_fusion_reads_filename,'w')
    data=[line[0].split('\t')+[str(len(line[1])), ','.join(line[1])] for line in fusion_reads.items()]
    data = sorted(data,key=lambda x:(int(x[2]),x[0],x[1]), reverse = True)
    data=[[hugo[line[0]],hugo[line[1]]]+line for line in data]
    data=['\t'.join(line)+'\n' for line in data]
    data.insert(0,'Fusion_gene_symbol_1\tFusion_gene_symbol_2\tFusion_gene_1\tFusion_gene_2\tCount_paired-end_reads\tSupporting_paired-read_ids\n')
    fo.writelines(data)
    fo.close()

    if options.output_fusion_reads_split_filename:
        xlist=[]
        xdata = data[:]
        xdata.pop(0)
        for xline in xdata:
            xl = xline.rstrip("\r\n").split("\t")
            kx = "%s.%s--%s__%s--%s" % (options.output_fusion_reads_split_filename,xl[0],xl[1],xl[2],xl[3])
            rx = []
            for ry in sorted(set(xl[5].split(","))):
                rx.append("%s/1\n" %(ry,))
                rx.append("%s/2\n" %(ry,))
            file(kx,"w").writelines(rx)
            xlist.append(kx+"\n")
        file(options.output_fusion_reads_split_filename,'w').writelines(xlist)

    if options.output_fusion_reads_simple_filename:
        fo = open(options.output_fusion_reads_simple_filename,'w')
        data = set()
        for line in fusion_reads.values():
            for l in line:
                data.add("%s/1\n"% (l,))
                data.add("%s/2\n"% (l,))
        data = sorted(data)
        fo.writelines(data)
        fo.close()

    if options.output_missing_mate_reads_filename:
        print "Writing...",options.output_missing_mate_reads_filename
        file(options.output_missing_mate_reads_filename,"w").writelines(missing_mates)


    print "The end."
    #
