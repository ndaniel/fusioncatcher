#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It analyzes the mappings of reads on exon-exon junctions.



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


# bowtie outputs one alignment per line. Each line is a collection of 8 fields separated by tabs; from left to right, the fields are:
#
# 1.   Name of read that aligned
#
# 2.   Reference strand aligned to, + for forward strand, - for reverse
#
# 3.   Name of reference sequence where alignment occurs, or numeric ID if no name was provided
#
# 4.   0-based offset into the forward reference strand where leftmost character of the alignment occurs
#
# 5.   Read sequence (reverse-complemented if orientation is -).
#
#      If the read was in colorspace, then the sequence shown in this column is the sequence of decoded nucleotides, not the original colors. See the Colorspace alignment section for details about decoding. To display colors instead, use the --col-cseq option.
#
# 6.   ASCII-encoded read qualities (reversed if orientation is -). The encoded quality values are on the Phred scale and the encoding is ASCII-offset by 33 (ASCII char !).
#
#      If the read was in colorspace, then the qualities shown in this column are the decoded qualities, not the original qualities. See the Colorspace alignment section for details about decoding. To display colors instead, use the --col-cqual option.
#
# 7.   If -M was specified and the prescribed ceiling was exceeded for this read, this column contains the value of the ceiling, indicating that at least that many valid alignments were found in addition to the one reported.
#
#      Otherwise, this column contains the number of other instances where the same sequence aligned against the same reference characters as were aligned against in the reported alignment. This is not the number of other places the read aligns with the same number of mismatches. The number in this column is generally not a good proxy for that number (e.g., the number in this column may be '0' while the number of other alignments with the same number of mismatches might be large).
#
# 8.  Comma-separated list of mismatch descriptors. If there are no mismatches in the alignment, this field is empty. A single descriptor has the format offset:reference-base>read-base. The offset is expressed as a 0-based offset from the high-quality (5') end of the read.

#########################
def line_from(a_map_filename):
    # it gives chunks from a_map_filename
    fin=open(a_map_filename,'r')
    while True:
        lines=fin.readlines(10**8)
        if not lines:
            break
        lines=[line.rstrip('\r\n').split('\t') for line in lines if line.rstrip('\r\n')]
        for line in lines:
            yield line
    fin.close()

#########################
def exon_exon_from(a_map_filename):
# it assumes that the MAP BOWTIE file is sorted by reference sequence column
    last_ee = None
    chunk = []
    chunk2 = []
    for line in line_from(a_map_filename):
        if not chunk:
            last_ee = line[2]
        if last_ee != line[2]:
        # line[2] is column no 3 in the BOWTIE MAP file which contains the
        # 'reference sequence' (i.e. exon-exon junctions) name
            yield (last_ee, chunk, chunk2)
            last_ee = line[2]
            chunk = []
            chunk2 = []
        r = len(line[4]) # read length
        rid = line[0] # read id
        mi = line[7].strip() # mismatches
        if mi:
            mi = len(mi.split(','))
        else:
            mi = 0
        start = int(line[3])+1  # start postion of mapped read on ref seq
        end = start+r-1 # end postion of mapped read on ref seq
        chunk.append((start,end))
        chunk2.append((rid,mi))
    if chunk:
        yield (last_ee,chunk, chunk2)

#########################
def todict(a_list):
    #parsing the transcript info file into a dictionary,
    # e.g. ex=123;ge=34feq;....
    d = dict()
    for el in a_list:
        a = el.partition('=')
        d[a[0]]=a[2]
    return d
#########################


if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It analyzes the mappings of reads on exon-exon junctions."""
    version="%prog 0.10 beta"

    parser = optparse.OptionParser(
                usage = usage,
                description = description,
                version = version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in BOWTIE MAP format where the rows
                         are sorted by the column containing the reference
                       sequence, which is gene-gene,transcript-transcript,
                       exon-exon. """)

    parser.add_option("--input_hugo",
                      action="store",
                      type="string",
                      dest="input_hugo_filename",
                      help="""The input database used for linking ENSEMBL GENE
                      ID to HUGO gene names.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output file containing the summary of the
                      reads mapping on the exon-exon junctions.""")

    parser.add_option("--output_henrik",
                      action="store",
                      type="string",
                      dest="output_henrik_filename",
                      help="""The output file containing extra information
                      regarding the reads' mappings on the exon-exon
                      junctions.""")

    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.input_hugo_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the options has not been specified.")
        sys.exit(1)

    print "Reading...",options.input_hugo_filename
    # it is assumed that it is in this format
    #    ensembl_gene_id
    #    hgnc_symbol
    #    description
    #    description
    hugo=dict()
    for el in [line.strip('\r\n').split('\t')[0:2]
                for line in file(options.input_hugo_filename,'r').readlines()
                    if line.strip('\r\n')]:
        g=el[0]
        h=el[1].replace(' ','_')
        if hugo.has_key(g):
            hugo[g]=hugo[g]+','+h
        else:
            hugo[g]=h

    print "Processing...",options.input_filename
    header=['5end_gene',
            '3end_gene',
            '5end_transcript',
            '3end_transcript',
            '5end_exon',
            '3end_exon',
            '5end_rank_exon_transcript',
            '3end_rank_exon_transcript',
            'counts_reads',
            'counts_unique_reads',
            'mean_lengths_of_shorter_overlaping_parts_of_reads_over_junction',
            'longest_of_shorter_overlaping_parts_of_reads_over_junction',
            'lengths_of_shorter_overlaping_parts_of_reads_over_junction',
            'average_number_mismatches_per_read',
            '5end_gene_symbol',
            '3end_gene_symbol',
            '.',
            '.',
            '.',
            'length_exon-exon_junction',
            'position_junction',
            'index_sequence_exon-exon_junction',
            '.',
            '.',
            '.',
            'id_sequence_in_fasta_file'
            ]
    data=[]
    flag=False
    if options.output_henrik_filename:
        flag=True
        henrik=[]

    exon_exon_unique = set()
    for (ee,positions,info) in exon_exon_from(options.input_filename):
        #example of ee
        #ENSG00000233001-ENSG00000133121;ENST00000439831-ENST00000336934;ENSE00001662779-ENSE00001798260;rank_exon=1-2;length=80;junction=40;index_seq=1
        temp=todict([el for el in ee.split(';') if el.find('=')!=-1])
        len_junct=int(temp['length'])
        junct=int(temp['junction'])
        idx=temp['index_seq']
        ranks=temp['rank_exon'].split('-')
        ux = ee.partition(';rank_exon=')[0].replace(';','-').split('-') # take everything up to rank_exon
        temp = ux[:]
        temp.extend(ranks) # add ranks

#        positions=[(i,j) for (i,j) in positions if i<=junct and j>=junct ] # test for overlap over junction
        posix = []
        infox = []
        for ix in xrange(len(positions)):
            if positions[ix][0] <= junct - 10 and positions[ix][1] >= junct + 10:
                posix.append(positions[ix])
                infox.append(info[ix])
        if not posix:
            continue

        cor=len(posix)
        uniq_posix=list(set(posix))
        uniq_cor=len(uniq_posix)
        uniq_infox = []
        for un in uniq_posix:
            uefa = sorted([infox[wk] for wk in xrange(len(posix)) if posix[wk]==un], key = lambda xyz:(xyz[1],xyz[0]))
            uniq_infox.append(uefa[0])


        d=[min(abs(junct-i+1),abs(j-junct)) for (i,j) in uniq_posix]
        avg=float(sum(d))/uniq_cor
        m = max(d)
        mm = ','.join(map(str,sorted(d,reverse=True)))
        mismatches = float(sum([mix for (rix,mix) in uniq_infox]))/len(uniq_infox)

        mykey = ee.split(';rank_exon=')[0].split(';')
        #mykey = mykey[0]+'-'+mykey[2]+'-'+str(cor)+'-'+str(uniq_cor)+'-'+str(avg)+'-'+str(m)
        mykey = mykey[0]+'-'+mm
        if mykey in exon_exon_unique:
            continue
        else:
            exon_exon_unique.add(mykey)
        if uniq_cor == 1 and ((mismatches > 0 and float(m) / mismatches <= 5) or (m > 14)):
            continue

        temp.extend([str(cor),str(uniq_cor),str(avg),str(m),mm,str(mismatches),hugo[temp[0]],hugo[temp[1]]]) # add counts_reads, counts_unique_reads, etc.
        temp.extend(['.','.','.',str(len_junct),str(junct),idx,'.','.','.',ee]) # add empty columns and the original id of the sequence
        data.append(temp)

        if flag:
            henrik.append(ee.split(';rank_exon=')[0].replace(';','-').replace('-','\t')+'\t'+','.join([str(x)+'-'+str(y) for (x,y) in posix])+'\t'+','.join([str(rix) for (rix,mix) in infox]))


    data=sorted(data,key=lambda x: (-float(x[9]),-float(x[10]),-float(x[11]),-float(x[8]),x[0],x[1],x[2],x[3],x[4],x[5]))
    data.insert(0,header)

    print "Writing...",options.output_filename
    file(options.output_filename,'w').writelines(['\t'.join(line)+'\n' for line in data])

    if flag:
        print "Writing...",options.output_henrik_filename
        file(options.output_henrik_filename,'w').writelines([line+'\n' for line in henrik])

    print "The end."
