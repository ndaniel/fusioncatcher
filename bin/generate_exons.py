#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a fasta file containing a genome and exonic positions on the 
genome and it gives the exon sequences as a fasta file.



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

#
"""
The format of the file contaning exonic positions is as follows, e.g. more_exons_ensembl.txt

ensembl_peptide_id
ensembl_gene_id
ensembl_transcript_id
ensembl_exon_id
exon_chrom_start
exon_chrom_end
rank
start_position
end_position
transcript_start
transcript_end
strand
chromosome_name
cds_start
cds_end
5_utr_start
5_utr_end
3_utr_start
3_utr_end

The output is in file: "exons_ensembl.fa" as follows:
>ENSE00001494058;length=97
CTTTCTTTAGCTTGCCCATGGTGATGTGAAGATGAGAAGAAATAGCAAGGCCCAACCAGT
TCTTCATCTGGAGACAGTTCAACGTTCTGCAAACCAG 

"""

import os
import sys
import Bio
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet
import optparse


def give_seq(seq,start,end,strand):
    if start > end:
        (start,end) = (end,start)
    s = ""
    if strand == '1' or strand == '+1' or strand == '+':
        s = seq[start-1:end]
    else:
        s = seq[start-1:end].reverse_complement()
    return s

def get_cols(t):
    return (t[12],t[3],int(t[4]),int(t[5]),t[11]) # chromosome, exon_id, exon_start, exons_end, exon_strand


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It takes as input a fasta file containing a genome and exonic positions on the genome and it gives the exon sequences as a fasta file."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_exons",
                      action="store",
                      type="string",
                      dest="input_exons_positions",
                      help="""Input file with exons positions.""")

    parser.add_option("--input_genome",
                      action="store",
                      type="string",
                      dest="input_genome_fasta",
                      help="""Input FASTA genome containing the genome sequences.""")


    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the exons sequences are written. Default is '%default'.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_exons_positions and 
            options.input_genome_fasta and 
            options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #



    data = set([get_cols(line.rstrip('\r\n').split('\t')) for line in file(options.input_exons_positions,'r').readlines() if line.rstrip('\r\n')])
    dc = dict()
    for line in data:
        if line[0] not in dc:
            dc[line[0]] = []
        dc[line[0]].append(line[1:])


    exons = []
    output_handle = open(os.path.join(options.output_directory,'exons.fa'), "w")
    for record in Bio.SeqIO.parse(open(options.input_genome_fasta, "rU"), "fasta") :
        print 'chromosome %s has length %d' % (record.id, len(record.seq))
        if record.id in dc:
            s = record.seq
            d = dc[record.id]
            for line in d:
                seq = give_seq(s,line[1],line[2],line[3])
                if seq:
                    exon = line[0]+';length=' + str(len(seq))
                    exons.append(Bio.SeqRecord.SeqRecord(seq,id=exon,name="",description=""))
            Bio.SeqIO.write(exons, output_handle, "fasta")
            exons = []
    output_handle.close()

