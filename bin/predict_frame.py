#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a list of chromosomal coordinates of fusion genes and it predicts
if their transcripts are in frame or not.

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2020 Daniel Nicorici

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
import gzip
import itertools
import Bio
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet


#
#
#
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def dna2prot(seq):
    seq = seq.strip().upper().replace('N','A')
    prot = ''

    for i in xrange(0,len(seq),3):
        r = codon_table.get(seq[i:i+3],'*')
        if r != '*':
            prot = prot + r
        else:
            break
    return prot


#
#
#
def predict(aexon,acds,achrom,aposition,astrand):
    txt = 'intergenic'
    # exon
    f = False
    info = []
    start = dict()
    if not aexon:
        txt = "---"
    else:
        for at in aexon:
            l = None
            start[at] = 0
            nuc = 0
            s = sum([q[2] - q[1] + 1 for q in aexon[at]])
            for e in aexon[at]:
                if astrand == e[3] and achrom == e[4]:
                    if l and e[1] > aposition and l < aposition and (not txt.startswith('exon')) and (not txt.startswith('UTR')) :
                        txt = 'intronic'
                    l = e[2]
                    if e[1] <= aposition and e[2] >= aposition:
                        if acds:
                            txt = 'UTR'
                            f  = True
                            if astrand == '+':
                                start[at] = nuc + aposition - e[1] + 1
                            elif astrand == '-':
                                start[at] = s - (nuc + aposition - e[1] + 1) + 1
                            else:
                                print >>sys.stderr,"WARNING: unknown strand found!"
                        elif not txt.startswith('UTR'):
                            txt = 'exonic(no-known-CDS)'
                        break
                    nuc = nuc + e[2] - e[1] + 1

        if acds:
            if f:
                for at in acds:
                    nuc = 0
                    ac = acds[at]
                    n = len(ac)
                    # length CDS
                    s = sum([q[2] - q[1] + 1 for q in ac])
                    # is the CDS reliable?
                    reliable = True if s % 3 == 0 else False

                    for i,e in enumerate(ac):
                        if astrand == e[3] and achrom == e[4]:
                            if aposition < e[1] and i == 0:
                                break
                            elif aposition > e[2]:
                                nuc = nuc + e[2] - e[1] + 1
                            elif e[1] <= aposition and e[2] >= aposition:
                                nuc = nuc + aposition - e[1] + 1

                                if not reliable and txt != 'CDS':
                                    txt = 'CDS(no-known-start-or-end)'
                                    continue

                                islastnucleotide = False
                                if i == n-1 and aposition == e[2]: # if it is then last exon and the last nucleotide
                                    islastnucleotide = True

                                if astrand == '+':
                                    # length of first UTR (from 5 to 3 on genome)
                                    utr = start[at] - nuc
                                    info.append((nuc,at,islastnucleotide,s,utr))
                                elif astrand == '-':
                                    islastnucleotide = True if nuc == 1 else False
                                    #length of first UTR (from 5 to 3 on genome)
                                    utr = start[at] - (s-nuc+1)
                                    info.append((s-nuc+1,at,islastnucleotide,s,utr))
                                else:
                                    print >>sys.stderr,"WARNING: unknown strand found!"
                                txt = 'CDS'
                                break
    # info - 0 : position in CDS
    # info - 1 : transcript id
    # info - 2 : is the last nucleotide? True or False
    # info - 3 : length of CDS
    # info - 4 : length UTR
    return (txt,info)


#
#
#
def add_line(aline,adict):
    chrom = aline[0]
    pstart = int(aline[3])
    pend = int(aline[4])
    strand = aline[6]
    if pstart > pend:
        (pstart,pend) = (pend,pstart)
    ids = [l.replace('"','').replace("'","").strip().split(' ') for l in aline[8].split(";") if l]
    ids = dict([l for l in ids if len(l) == 2])
    g_id = ids["gene_id"]
    t_id = ids["transcript_id"]
    e_no = int(ids["exon_number"])
    if g_id not in adict:
        adict[g_id] = dict()
    if t_id not in adict[g_id]:
        adict[g_id][t_id] = []
    adict[g_id][t_id].append((e_no,pstart,pend,strand,chrom))


#
#
#
def predict_frame(gtf_file,
                  transcripts_file,
                  input_file,
                  output_file,
                  compress_transcripts = True,
                  verbose = True):


    if verbose:
        print >>sys.stderr,"Reading the list of fusion genes..."
    # read the list of fusion genes and their chromosomal positions
    data = [line.rstrip('\r\n').split('\t') for line in file(input_file,'r').readlines() if line.rstrip('')]
    fusion1 = [tuple([line[10]] + line[8].split(':')) for line in data[1:]]
    fusion2 = [tuple([line[11]] + line[9].split(':')) for line in data[1:]]
    myg = set([line[10] for line in data[1:]]+[line[11] for line in data[1:]])

    if myg:

        tr2fa = dict()
        if transcripts_file:
            if verbose:
                print >>sys.stderr,"Reading the transcripts' sequences..."
            handle = open(transcripts_file,"rU")
            for record in Bio.SeqIO.parse(handle, "fasta"):
                temp = record.id.partition(';')
                t = temp[0]
                g = temp[2]
                if g in myg:
                    tr2fa[t] = str(record.seq).upper()
            handle.close()

        # read and pre-process the GTF file
        if verbose:
            print >>sys.stderr,"Parsing the GTF file..."
        g = []
        if gtf_file == '-':
            g = [line.rstrip('\r\n').split("\t") for line in file(sys.stdin,"r").readlines() if (not line.startswith("#")) and line.rstrip()]
        elif gtf_file.endswith('.gz'):
            z = gzip.open(gtf_file,"r")
            g = [line.rstrip('\r\n').split("\t") for line in z.readlines() if (not line.startswith("#")) and line.rstrip()]
            z.close()
        else:
            g = [line.rstrip('\r\n').split("\t") for line in file(gtf_file,"r").readlines() if (not line.startswith("#")) and line.rstrip()]



        # get all exons per gene as a dictionary
        if verbose:
            print >>sys.stderr,"Building the database of exons and CDSes..."
        exon = dict()
        cds = dict()
        q = 0
        for line in g:
            if line[2] == 'gene':
                q = q + 1
            flag = False
            x = line[8].partition(";")
            if x[0] and x[0].startswith("gene_id "):
                x = x[0][9:-1]
            if x not in myg:
                continue
            if line[2] == 'exon':
                add_line(line,exon)
            elif line[2] == 'CDS':
                add_line(line,cds)
        # sorting the database
        for g in exon:
            for t in exon[g]:
                exon[g][t].sort(key=lambda x: x[1]) # order by start position
        for g in cds:
            for t in cds[g]:
                cds[g][t].sort(key=lambda x: x[1])

        if verbose:
            print >>sys.stderr,"Predicting effect..."
        # predict frame
        result = ['Predicted_effect']
        transcript = ['Predicted_fused_transcripts']
        protein = ['Predicted_fused_proteins']
        transcript_full = ['Predicted_fused_transcripts']
        for i in xrange(len(fusion1)):
            f1 = fusion1[i]
            f2 = fusion2[i]
            (t1,info1) = predict(exon.get(f1[0],None),cds.get(f1[0],None),f1[1],int(f1[2]),f1[3])
            (t2,info2) = predict(exon.get(f2[0],None),cds.get(f2[0],None),f2[1],int(f2[2]),f2[3])
            tr = ''
            t3 = ''
            pr = ''
            trf = ''
            if t1 == t2 and t1 == 'CDS' and info1 and info2:
                r = [(elem[0][1], # transcript 1 id
                      elem[0][4]+elem[0][0], # end position in transcript 1
                      elem[1][1], # transcript 2 id
                      elem[1][4]+elem[1][0], # # start position in transcript 2
                      ((elem[0][0]%3)-((elem[1][0]-1)%3)), # just test if they are the same ## is divisible by 3?
                      elem[0][4], # length UTR in transcript 1
                      elem[1][4]  # length UTR in transcript 2
                      ) for elem in itertools.product(info1,info2)]
                rr = ["%s:%d/%s:%d" % (el[0],el[1],el[2],el[3]) for el in r if el[4] == 0]
                if rr:
                    t3 = 'in-frame'
                    tr = ';'.join(rr)
                    if tr2fa:
                        pr = []
                        trf = []
                        for el in r:
                            if el[4] == 0:
                                s1 = tr2fa.get(el[0],None)
                                s2 = tr2fa.get(el[2],None)
                                if s1 and s2:
                                    p1 = s1[el[5]:el[1]]
                                    p2 = s2[el[3]-1:]
                                    t = p1+p2
                                    p = dna2prot(t)
                                    #pr.append("%s/%s*%s"% (p,p1,p2))
                                    pr.append("%s" % (p,))
                                    trf.append("%s" % (t,))
                                else:
                                    pr.append("")
                                    trf.append("")
                        pr = ';'.join(pr)
                        trf = ';'.join(trf)
                else:
                    t3 = 'out-of-frame'
                    rr = ["%s:%d/%s:%d" % (el[0],el[1],el[2],el[3]) for el in r if el[4] != 0]
                    tr = ';'.join(rr)
                    if tr2fa:
                        pr = []
                        trf = []
                        for el in r:
                            if el[4] != 0:
                                s1 = tr2fa.get(el[0],None)
                                s2 = tr2fa.get(el[2],None)
                                if s1 and s2:
                                    p1 = s1[el[5]:el[1]]
                                    p2 = s2[el[3]-1:]
                                    t = p1+p2
                                    p = dna2prot(t)
                                    #pr.append("%s/%s*%s"% (p,p1,p2))
                                    pr.append("%s"% (p,))
                                    trf.append("%s" % (t,))
                                else:
                                    pr.append("")
                                    trf.append("")
                        pr = ';'.join(pr)
                        trf = ';'.join(trf)
            else:
                if t1 == 'CDS' and t2 != 'CDS':
                    t1 = 'CDS(truncated)'
                    for elk in info1:
                        if elk[2]:
                            t1 = 'CDS(complete)'
                            break
                elif t2 == 'CDS' and t1 != 'CDS':
                    t2 = 'CDS(truncated)'
                    for elk in info2:
                        if elk[0] == 1:
                            t2 = 'CDS(complete)'
                            break
            y = ''
            if t3:
                y = t3
            else:
                y = "%s/%s" % (t1,t2)
            result.append(y)
            transcript.append(tr)
            protein.append(pr)
            transcript_full.append(trf)
        data = zip(*data)
        data.append(result)
        if compress_transcripts:
            data.append(transcript)
        else:
            data.append(transcript_full)
        data.append(protein)
        data = zip(*data)
        file(output_file,'w').writelines(['\t'.join(line)+'\n' for line in data])
    else:
        file(output_file,'w').write('\t'.join(data[0])+'\n')



################################################################################
################################################################################
################################################################################
def main():

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options]"

    description = """It takes as input a list of chromosomal coordinates of fusion genes and it predicts
if their transcripts are in frame or not."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2020 Daniel Nicorici

"""



    version = "%prog 0.04 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
                            )



    parser.add_option("-g","--gtf",
                      action = "store",
                      type = "string",
                      dest = "gtf_filename",
                      help = """The input GTF file containing the genome annotation.""")

    parser.add_option("-t","--transcripts",
                      action = "store",
                      type = "string",
                      dest = "transcripts_filename",
                      help = """The input FASTA file contains the DNA sequences of all ENSEMBL transcripts.""")

    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file containing the fusion (chromosomal) coordinates for each fusion genes.""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file where the frame predictions are written. """)

    parser.add_option("-c", "--compress-transcripts",
                      action = "store_true",
                      dest = "compress_transcripts",
                      default = False,
                      help = "Compress the transcript sequences.")

    parser.add_option("-q", "--quiet",
                      action = "store_false",
                      dest = "verbose",
                      default = True,
                      help = "Do not print status messages to stdout.")

    parser.add_option("-a","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.gtf_filename and
            options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")


    #
    predict_frame(options.gtf_filename,
                  options.transcripts_filename,
                  options.input_filename,
                  options.output_filename,
                  options.compress_transcripts,
                  options.verbose)


if __name__ == '__main__':
    main()
    
    
#    
