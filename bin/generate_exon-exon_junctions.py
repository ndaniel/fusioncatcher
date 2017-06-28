#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It generates the exon-exon junctions from transcripts for a list of genes.



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
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet

def ispoly(sec):
    f = False
    if sec:
        if sec[0] * len(sec) == sec:
            f = True
    return f

def todict(a_list):
    #parsing the transcript info file into a dictionary, e.g. ex=123;ge=34feq;....
    d = dict()
    for el in a_list:
        a = el.partition('=')
        d[a[0]] = a[2]
    return d


if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It generates the exon-exon junctions from transcripts for a list of genes."""
    version="%prog 0.11 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_fusion_genes",
                      action="store",
                      type="string",
                      dest="input_fusion_genes_filename",
                      help="""The input file in text tab delimited format containing the fusion genes candidates produced by 'extract_fusion_genes.py'. This is optional and if it is not specified all exon-exon junctions will be generated for all genes from the database.""")

    parser.add_option("--input_fasta_transcripts",
                      action="store",
                      type="string",
                      dest="input_fasta_filename",
                      help="""The input FASTA files containing the transcripts, e.g. data/ensembl/transcripts.fa.""")

    parser.add_option("--input_database_transcripts",
                      action="store",
                      type="string",
                      dest="input_database_filename",
                      help="""The input text file containg information regarding the transcripts, e.g. data/ensembl/transcripts.txt.""")

    parser.add_option("--overlap_read",
                      action="store",
                      type="int",
                      dest="overlap_read",
                      default=15,
                      help="""The minimum length of the overlapping region between the read the exon-exon juntion.
Length_of_the_exon-exon_juntion = 2 * (length_reads - overlap_read). The joint point is at the middle point of the exon-exon junction. Default value is %default.""")

    parser.add_option("--length_reads_filename",
                      action="store",
                      type="string",
                      dest="length_reads_filename",
                      help="""A text file containing on the first line the length of the reads.""")

    parser.add_option("--length_reads",
                      action="store",
                      type='int',
                      dest="length_reads",
                      help="""The length of the reads.""")

    parser.add_option("--unique_cut_sequences",
                      action="store_true",
                      dest="unique_cut_seq",
                      default=False,
                      help="""It outputs only the unique cut sequences (of exon-exon junctions). Default value is %default. This is a dangerous option and it is highly recommended to be set on false always.!""")

    parser.add_option("--unique_cut_sequences_same_pair",
                      action="store_true",
                      dest="unique_cut_seq_same_pair",
                      default=False,
                      help="""It outputs only the unique cut sequences (of exon-exon junctions) within the given pair of genes. Default value is %default.""")

    parser.add_option("--output_cut_junction",
                      action="store",
                      type="string",
                      dest="output_cut_filename",
                      help="""A FASTA file containing the exon-exon junctions where the cutting is done according to the options '--overlap_read' and '--length_reads'.""")

    parser.add_option("--output_full_junction",
                      action="store",
                      type="string",
                      dest="output_full_filename",
                      help="""A FASTA file containing the exon-exon junctions where the cutting is not done.""")

    parser.add_option("--output_unique_cut_sequences_same_pair",
                      action="store",
                      type='string',
                      dest="output_unique_cut_seq_same_pair_filename",
                      help="""In case the option '--unique_cut_sequences_same_pair' is used it outputs in a file all the names of the sequences for which the sequences are the same.""")

    parser.add_option("--output_count_seq",
                      action="store",
                      type='string',
                      dest="output_count_seq_filename",
                      help="""If used then the number of sequences from the output FASTA file (i.e. --output_cut_junction) will be reported.""")

    parser.add_option("--output_count_nuc",
                      action="store",
                      type='string',
                      dest="output_count_nuc_filename",
                      help="""If used then the number of nucleotides of all sequences from the output FASTA file (i.e. --output_cut_junction) will be reported.""")


    (options,args)=parser.parse_args()

    # validate options
    if not (
            options.input_fasta_filename and
            options.input_database_filename
            ):
        parser.print_help()
        parser.error("One of the options '--input_fasta_transcripts' or '--input_database_transcripts' has not been specified.")

    if (options.length_reads and
        options.length_reads_filename
       ):
        parser.error("Only one of the options '--length_reads_filename' and '--length_reads' should be used and not both!")

    if (options.unique_cut_seq and
        options.unique_cut_seq_same_pair
        ):
        parser.error("Only one of the options '--unique_cut_sequences' and '--unique_cut_sequences_same_pair' should be set to true and not both!")

    if (options.output_unique_cut_seq_same_pair_filename and
        options.unique_cut_seq_same_pair
        ):
        parser.error("The option '--output_unique_cut_sequences_same_pair' can be used only if the option '--unique_cut_sequences_same_pair' is set to 'true'!")

    cx=0
    if options.length_reads or options.length_reads_filename:
        cx=cx+1
    if options.output_cut_filename:
        cx=cx+1
    if cx>0 and cx<2:
        parser.error("All or none of the options '--overlap_read' and '--length_reads'/'--length_reads_filename' and '--output_cut_filename' should be used!")

    if not (options.output_cut_filename or options.output_full_filename):
        parser.error("No output specified!")

    print "Starting..."

    overlap=options.overlap_read

    length_reads=None
    if options.length_reads_filename:
        print "Reading...",options.length_reads_filename
        length_reads=int(file(options.length_reads_filename,'r').readline().strip())
    elif options.length_reads:
        length_reads=options.length_reads


    print "Reading...",options.input_fasta_filename
# Asume:
# --input_fasta_transcripts 'more_transcript.fa'
#
#>tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103
#ACGGGGGCGGGCCCGCGGTGACGTCGGGAGGGCAGCGACGCGCGGAGGCGGCGGCGGAGC
#CTCCTCCTGCTGCTGCTGCGCCCCATCCCCCCGCGGCCGGCCAGTTCCAGCCCGCACCCC
#GCGTCGGTGCCCGCGCCCCTCCCCGGGCCCCGCCATGGGCCTCACCGTGTCCGCGCTCTT
#TTCGCGGATCTTCGGGAAGAAGCAGATGCGGATTCTCATGGTTGGCTTGGATGCGGCTGG
#...
    handle=open(options.input_fasta_filename,"rU")
    transcript=dict()
    gene=dict()
    for record in Bio.SeqIO.parse(handle, "fasta"):
        # temp=todict(record.id.split(';'))
        #t = temp['tr']
        #g = temp['ge']
        temp = record.id.partition(';')
        t = temp[0]
        g = temp[2]
        transcript[t]=str(record.seq).upper()
        if not gene.has_key(g):
            gene[g]=[]
        gene[g].append(t)
    handle.close()

# --input_database_transcripts 'more_transcript.txt'
#tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103	tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103;ex=ENSE00001872691,r=1,sc=127228399,ec=127228619,st=1,et=221;ex=ENSE00000720374,r=2,sc=127229137,ec=127229217,st=222,et=302;ex=ENSE00001737440,r=3,sc=127229539,ec=127229648,st=303,et=412;ex=ENSE00000720381,r=4,sc=127230120,ec=127230191,st=413,et=484;ex=ENSE00000720384,r=5,sc=127231017,ec=127231142,st=485,et=610;ex=ENSE00000882271,r=6,sc=127231267,ec=127231759,st=611,et=1103;cds_st=;cds_et=;cds_sc=;cds_ec=;
#...
    print "Reading...",options.input_database_filename
    info_tr=[line.rstrip('\r\n').split('\t') for line in file(options.input_database_filename,'r') if line.rstrip('\r\n')]
    #info_tr=dict([tuple([todict(line[0].split(';'))['tr'],line[1]]) for line in info_tr])
    info_tr = dict([ ([line[0].partition(';')[0], line[1]]) for line in info_tr])
    for k in info_tr.keys():
        v=info_tr[k]
        t=todict([x for x in v.split(';') if x.find(',')==-1])
#        t['len']=int(t['len'])
        exons=[]
        for el in [x for x in v.split(';') if x.find(',')!=-1]:
            tt=todict(el.split(','))
            tt['st']=int(tt['st'])
            tt['et']=int(tt['et'])
#            tt['sc']=int(tt['sc'])
#            tt['ec']=int(tt['ec'])
            tt['r']=int(tt['r'])
            exons.append(tt)
        exons=sorted(exons,key=lambda x: (x['r']))
        t['exons']=exons
        info_tr[k]=t

    if options.input_fusion_genes_filename:
        print "Reading...",options.input_fusion_genes_filename
        pairs=set([line.rstrip('\r\n') for line in file(options.input_fusion_genes_filename,'r') if line.rstrip('\r\n')])
        pairs=sorted([line.split('\t') for line in pairs])
    else: # generate a list of all gene-gene pairs for exon-exon combinations
        print "Generating a list of all gene-gene pairs for exon-exon combinations..."
        pairs=[]
        for g,t in gene.iteritems():
            if ([1 for tt in t if len(info_tr[tt]['exons'])>1]): # There should be more than one exon per transcript
                pairs.append([g,g])
        pairs=sorted(pairs)
        print len(pairs),"pairs of genes generated!"

    print "Generating exon-exon junctions..."
    seq_index=dict()
    sequences_cut=[]
    i=0
    handle_full=None
    if options.output_full_filename:
        handle_full=open(options.output_full_filename,"w")
    if options.output_cut_filename:
        handle_cut=open(options.output_cut_filename,"w")
        sequences_cut=[]
    if options.output_cut_filename:
        cut_size=length_reads-overlap

    options_output_unique_cut_seq_same_pair_filename=False
    if options.output_unique_cut_seq_same_pair_filename:
        options_output_unique_cut_seq_same_pair_filename=True
        handle_id_links=open(options.output_unique_cut_seq_same_pair_filename,'w')

    options_output_full_filename=False
    if options.output_full_filename:
        options_output_full_filename=True
    options_output_cut_filename=False
    if options.output_cut_filename:
        options_output_cut_filename=True
    options_unique_cut_seq_same_pair=False
    if options.unique_cut_seq_same_pair:
        options_unique_cut_seq_same_pair=True
    options_unique_cut_seq=False
    if options.unique_cut_seq:
        options_unique_cut_seq=True

    count_seq = 0
    count_nuc = 0
    uniq_tt_g=set()
    uniq_seq_same_pair=dict()
    for (a,b) in pairs: # for each pair
        i=i+1
        if i%1000==0:
            print "...done",i,"gene-gene pairs out of a total of",len(pairs)
        uniq_seq_same_pair=dict()
        for ta in gene[a]: # transcript A
            sa=transcript[ta] # sequence transcript A
            ia=info_tr[ta]['exons'] # information transcript A
            for tb in gene[b]: # transcript B
                sb=transcript[tb] # sequence transcript B
                ib=info_tr[tb]['exons'] # information transcript B
                #print '(',a,',',b,')','(',ta,',',tb,')','has (',len(ia),',',len(ib),') exons'
                if a==b: # avoid some cases of combinations done for the pair of the same gene
                    ttg=tuple((sorted([ta,tb]))+[a])
                    if ttg in uniq_tt_g:
                        continue
                    else:
                        uniq_tt_g.add(ttg)
                for ea in ia:
                    eas=ea['st'] # exon A start
                    eae=ea['et']
                    ean=ea['ex']
                    ear=ea['r']
                    for eb in ib:
                        ebs=eb['st'] # exon B start
                        ebe=eb['et']
                        ebn=eb['ex']
                        ebr=eb['r']
                        if ean==ebn:
                            continue
                        if a==b and ta==tb and ear>=ebr:
                                continue
                        #print ean,ebn,ear,ebr
                        if options_output_full_filename:
                            sequences_full=[]
                            # first exon-exon junction
                            seq=sa[:eae]+sb[ebs-1:]
                            junction_id=a+'-'+b+';'+ta+'-'+tb+';'+ean+'-'+ebn+';rank_exon='+str(ear)+'-'+str(ebr)+';length='+str(len(seq))+';junction='+str(eae)
                            sequences_full.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq,Bio.Alphabet.IUPAC.ambiguous_dna),id=junction_id,name="",description=""))
                            # second exon-exon junction
                            seq=sb[:ebe]+sa[eas-1:]
                            junction_id=b+'-'+a+';'+tb+'-'+ta+';'+ebn+'-'+ean+';rank_exon='+str(ebr)+'-'+str(ear)+';length='+str(len(seq))+';junction='+str(ebe)
                            sequences_full.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq,Bio.Alphabet.IUPAC.ambiguous_dna),id=junction_id,name="",description=""))
                            # writing
                            Bio.SeqIO.write(sequences_full,handle_full,"fasta")
                        if options_output_cut_filename:
                            #
                            # first exon-exon junction
                            #
                            aa=eae-cut_size
                            if aa<0:
                                aa=0
                            bb=ebs-1+cut_size
                            if bb>len(sb):
                                bb=len(sb)
                            s1=sa[aa:eae].upper()
                            s2=sb[ebs-1:bb].upper()
                            seq=s1+s2

                            if ispoly(s1) or ispoly(s2):
                                doit = False

                            idx=None # the sequence that are the same have the same index

                            doit=True
                            flag=False
                            if seq_index.has_key(seq):
                                idx=seq_index[seq]
                                flag=True
                            else:
                                idx=len(seq_index)
                                seq_index[seq]=idx

                            junction_id=a+'-'+b+';'+ta+'-'+tb+';'+ean+'-'+ebn+';rank_exon='+str(ear)+'-'+str(ebr)+';length='+str(len(seq))+';junction='+str(len(s1))+';index_seq='+str(idx)

                            if options_unique_cut_seq_same_pair:
                                if uniq_seq_same_pair.has_key(seq):
                                    uniq_seq_same_pair[seq].append(junction_id)
                                    doit=False
                                else:
                                    uniq_seq_same_pair[seq]=[junction_id]
                            if options_unique_cut_seq and flag:
                                doit=False

                            if doit:
                                sequences_cut.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq,Bio.Alphabet.IUPAC.ambiguous_dna),id=junction_id,name="",description=""))
                                count_seq = count_seq + 1
                                count_nuc = count_nuc + len(seq)
                            #
                            # second exon-exon junction
                            #
                            aa=ebe-cut_size
                            bb=eas-1+cut_size
                            if aa<0:
                                aa=0
                            if bb>len(sa):
                                bb=len(sa)
                            s1=sb[aa:ebe].upper()
                            s2=sa[eas-1:bb].upper()
                            seq=s1+s2

                            if ispoly(s1) or ispoly(s2):
                                doit = False

                            idx=None # the sequence that are the same have the same index

                            doit=True
                            flag=False
                            if seq_index.has_key(seq):
                                idx=seq_index[seq]
                                flag=True
                            else:
                                idx=len(seq_index)
                                seq_index[seq]=idx

                            junction_id=b+'-'+a+';'+tb+'-'+ta+';'+ebn+'-'+ean+';rank_exon='+str(ebr)+'-'+str(ear)+';length='+str(len(seq))+';junction='+str(len(s1))+';index_seq='+str(idx)

                            if options_unique_cut_seq_same_pair:
                                if uniq_seq_same_pair.has_key(seq):
                                    uniq_seq_same_pair[seq].append(junction_id)
                                    doit=False
                                else:
                                    uniq_seq_same_pair[seq]=[junction_id]
                            if options_unique_cut_seq and flag:
                                doit=False

                            if doit:
                                sequences_cut.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq,Bio.Alphabet.IUPAC.ambiguous_dna),id=junction_id,name="",description=""))
                                count_seq = count_seq + 1
                                count_nuc = count_nuc + len(seq)
                            # write the cut sequences
                            if len(sequences_cut)>100000:
                                #print "...saving fasta cuts..."
                                Bio.SeqIO.write(sequences_cut,handle_cut,"fasta")
                                sequences_cut=[]
                                #print "...finished saving..."
        if options_output_unique_cut_seq_same_pair_filename and uniq_seq_same_pair:
            a_temp=[]
            for (a_seq,a_list) in uniq_seq_same_pair.iteritems():
                first_id=a_list[0]
                my_set=list(sorted(set(a_list).difference(set([first_id]))))
                a_temp.append(first_id+'\t'+'|'.join(my_set)+'\n')
            handle_id_links.writelines(a_temp)
            del a_temp


    if options.output_unique_cut_seq_same_pair_filename:
        handle_id_links.close()
    if options.output_full_filename:
        handle_full.close()

    if options.output_cut_filename:
        if sequences_cut:
            Bio.SeqIO.write(sequences_cut,handle_cut,"fasta")
        handle_cut.flush()
        handle_cut.close()
        if os.path.getsize(options.output_cut_filename) == 0:
            file(options.output_cut_filename,"w").write(">empty\nACGTACGTACGTACGTA")

    if options.output_count_seq_filename:
        file(options.output_count_seq_filename,"w").write("%d" % (count_seq,))

    if options.output_count_nuc_filename:
        file(options.output_count_nuc_filename,"w").write("%d" % (count_nuc,))

    print "The end."
    #
