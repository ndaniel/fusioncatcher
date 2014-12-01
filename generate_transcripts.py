#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It ansembles the transcripts sequences from exons. Extra information about slice sites is included.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2014 Daniel Nicorici

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
import datetime
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet
import gc


##########################################################
def chunk_from(a_list,column):
    # it assumes that a_list is ordered asceding on 'column'
    deposit=[]
    if a_list:
        a=a_list[0][column]
        for elem in a_list:
            if a==elem[column]:
                deposit.append(elem)
            else:
                yield deposit
                deposit=[elem]
                a=elem[column]
        if deposit:
            yield deposit


#####################################################################
#####################################################################
#####################################################################
# main

if __name__=='__main__':
    print "Starting..."

    #command line parsing

    usage="%prog [options]"
    description="""It ansembles the transcripts sequences from exons. Extra information about slice sites is included."""
    version="%prog 0.11 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_fasta_exons",
                      action="store",
                      type="string",
                      dest="exons_fasta_filename",
                      help="""A FASTA file containing all the exon sequences (e.g. 'exons.fa').""")

    parser.add_option("--input_database",
                      action="store",
                      type="string",
                      dest="database_filename",
                      help="""A text file containg all the information regarding exons, genes, proteins and their positions (e.g. 'exons.txt')""")

    parser.add_option("--output_fasta",
                      action="store",
                      type="string",
                      dest="transcripts_fasta_filename",
                      help="""The output file containing all the assembled transcripts in FASTA format.""")

    parser.add_option("--output_extra",
                      action="store",
                      type="string",
                      dest="transcripts_extra_filename",
                      help="""Extra information about transcripts.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the assembled transcripts are written. Default is '%default'.""")




    (options,args)=parser.parse_args()

    # validate options
    if not (options.exons_fasta_filename and
            options.database_filename and
            options.output_directory
            ):
        parser.print_help()
        parser.error("Missing argument(s)!")
        sys.exit(1)

    if not options.transcripts_fasta_filename:
        options.transcripts_fasta_filename = os.path.join(options.output_directory,'transcripts.fa')

    if not options.transcripts_extra_filename:
        options.transcripts_extra_filename = os.path.join(options.output_directory,'transcripts.txt')


    print "Reading the exon sequences...", options.exons_fasta_filename
    # build a dictionay out of exons
    handle=open(options.exons_fasta_filename, "rU")
    exon_dict=dict()
    gc.disable()
    for record in Bio.SeqIO.parse(handle, "fasta"):
        exon_dict[record.id.split(';')[0].upper()]=record.seq.tostring().upper()
    gc.enable()
    handle.close()

    print "Reading the exon information...", options.database_filename
    gc.disable()
    database=[line.rstrip('\r\n').split('\t') for line in file(options.database_filename,'r').readlines() if line.rstrip('\r\n')]
    gc.enable()
    # the database is expected to have the following columns
    # 1: ensembl_peptide_id
    # 2: ensembl_gene_id
    # 3: ensembl_transcript_id
    # 4: ensembl_exon_id
    # 5: exon_chrom_start
    # 6: exon_chrom_end
    # 7: rank (of exon in the transcript)
    # 8: start_position (of gene)
    # 9: end_position (of gene)
    # 10: transcript_start
    # 11: transcript_end
    # 12: strand (of chromosome)
    # 13: chromosome_name
    # 14: cds_start
    # 15: cds_end
    columns=[
        'ensembl_peptide_id',
        'ensembl_gene_id',
        'ensembl_transcript_id',
        'ensembl_exon_id',
        'exon_chrom_start',
        'exon_chrom_end',
        'rank',
        'start_position',
        'end_position',
        'transcript_start',
        'transcript_end',
        'strand',
        'chromosome_name'
#        'cds_start',
#        'cds_end',
#        '5_utr_start',
#        '5_utr_end',
#        '3_utr_start',
#        '3_utr_end'
        ]
    col=dict([(t,i) for i,t in enumerate(columns)])
    # sort ascendigly the database on 'transcript' and 'exon rank'
    tr=col['ensembl_transcript_id']
    ra=col['rank']
    database=sorted(database, key=lambda x:(x[tr], int(x[ra])))
    sequences=[]
    extra=[]

    for chunk in chunk_from(database,tr):
        transcript=chunk[0][col['ensembl_transcript_id']]
        gene=chunk[0][col['ensembl_gene_id']]
        protein=chunk[0][col['ensembl_peptide_id']]
        chromosome=chunk[0][col['chromosome_name']]
        strand=chunk[0][col['strand']]
        if strand=='1' or strand=='+1':
            strand='+'
        elif strand=='-1' or strand=='-':
            strand='-';
        info_exon=[]
        s=[]
        length_transcript=1

        for exon in chunk:
            exon_id=exon[col['ensembl_exon_id']]
            exon_seq=exon_dict[exon_id.upper()]
            exon_length=len(exon_seq)
            s.append(exon_seq)
            sc=exon[col['exon_chrom_start']]
            ec=exon[col['exon_chrom_end']]
            st=length_transcript
            et=length_transcript+exon_length-1
            rank = exon[col['rank']]
            info=['ex='+exon_id,
                  'r='+rank,
                  'sc='+sc,
                  'ec='+ec,
                  'st='+str(st),
                  'et='+str(et)
                    ]
            length_transcript=length_transcript+exon_length
            info=','.join(info)
            info_exon.append(info)

        s=''.join(s)
        id = "%s;%s" % (transcript,gene)
#            ['tr='+transcript,
#             'ge='+gene
#             'pn='+protein,
#             'chr='+chromosome,
#             'str='+strand,
#             'len='+str(len(s))
#             ]
        idplus = ['tr='+transcript,
                  'ge='+gene,
                  'pn='+protein,
                  'chr='+chromosome,
                  'str='+strand,
                  'len='+str(len(s))
                 ]


        #id1=';'.join(id)
        id1 = id
        #id2=';'.join([id]+idplus+info_exon)
        id2=';'.join(idplus+info_exon)
        if len(s) > 100:
            sequences.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s,Bio.Alphabet.IUPAC.ambiguous_dna),id=id1,name="",description=""))
            extra.append(id1+'\t'+id2+'\n')


    print "Writing the transcript sequences...",options.transcripts_fasta_filename
    handle=open(options.transcripts_fasta_filename,"w")
    Bio.SeqIO.write(sequences,handle,"fasta")
    handle.close()


    print "Writing the extra information...",options.transcripts_extra_filename
    file(options.transcripts_extra_filename,'w').writelines(extra)
    #
