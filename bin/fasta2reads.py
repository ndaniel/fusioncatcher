#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a file in FASTA format and it covers uniformly with paired-reads which are written in a FASTQ file.



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
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.Alphabet
import string
import gzip
import gc

ttable = string.maketrans("ACGTYRSWKMBDHV-.","TGCARYSWMKVHDB-.") # global

#
#
#
def reversecomplement(seq):
    return seq.rstrip('\r\n').translate(ttable)[::-1]

#
#
#
class tofastq:

    def __init__(self,file_name,size_buffer=10**8):
        self.file_name = file_name
        if file_name:
            if file_name.lower().endswith('.gz'):
                self.file_handle = gzip.open(file_name,'w')
            else:
                self.file_handle = open(file_name,'w')
        self.size_buffer = size_buffer
        self.data = []
        self.size = 0
        
    def add_line(self,line):
        #line = line.rstrip('\r\n')+'\n'
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size = self.size + len(line)
        if self.size>self.size_buffer:
            self.__write_buffer()
            
    def add_lines(self,lines):
        gc.disable()
        self.size = self.size + sum([len(line) for line in lines])
        self.data.extend(lines)
        gc.enable()
        if self.size > self.size_buffer:
            self.__write_buffer()
            
    def __write_buffer(self):
        self.file_handle.writelines(self.data)
        self.size = 0
        self.data = []
        
    def is_filename_valid(self):
        if self.file_name:
            return True
        else:
            return False
            
    def close(self):
        if self.is_filename_valid():
            if self.data:
                self.__write_buffer()
            self.file_handle.close()
            self.file_name = None
            
    def __del__(self):
        self.close()


def fasta2reads(
        input_filename,
        output1_filename,
        output2_filename,
        output3_filename = None,
        step = 1,
        gap = 50,
        length = 100,
        snp = 0
        ):
        

    fq1 = tofastq(output1_filename)
    fq2 = tofastq(output2_filename)
    fasta = None
    if output3_filename:
        fasta = open(output3_filename, "w")
    q = "I"*length
    i = -1
    a = -1
    fa = []
    f = 2 * length + gap
    b = None
    if gap > -1:
        b = 'N'*gap
    k = 0
    for record in Bio.SeqIO.parse(open(input_filename, "r"), "fasta") :
        a = a + 1
        s = str(record.seq)
        s = s.upper()
        n = len(s)
        k = 0
        for j in xrange(0,n-f+1,step):
            k = k + 1
            i = i + 1
            r1 = s[j:j+length]
            r2 = s[j+f-length:j+f]
            r2r = reversecomplement(r2)
            fq1.add_line("@%d\n%s\n+\n%s\n" % (i,r1,q))
            fq2.add_line("@%d\n%s\n+\n%s\n" % (i,r2r,q))
            if fasta:
                if gap > 0:
                    fa.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("%s%s%s" % (r1,b,r2)),id=str(i),name="",description=""))
                else:
                    fa.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(s[j:j+f]),id=str(i),name="",description=""))
                if fa and len(fa) > 10000:
                    Bio.SeqIO.write(fa, fasta, "fasta")
                    fa = []
        print >>sys.stderr, "  * Sequence no.",a+1,"of length", n,"bp has produced",k,"paired-reads."
        
    if fasta:
        if fa:
            Bio.SeqIO.write(fa, fasta, "fasta")
        fasta.close()
    fq1.close()
    fq2.close()


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It takes as input a file in FASTA format and it covers uniformly with paired-reads which are written in a FASTQ file."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage = usage, description = description, version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTA format.""")


    parser.add_option("--output1","-1",
                      action="store",
                      type="string",
                      dest="output1_filename",
                      help="""The output file in FASTQ format.""")

    parser.add_option("--output2","-2",
                      action="store",
                      type="string",
                      dest="output2_filename",
                      help="""The output file in FASTQ format.""")

    parser.add_option("--output3","-3",
                      action="store",
                      type="string",
                      dest="output3_filename",
                      help="""The output file in FASTA format.""")


    parser.add_option("-s","--step",
                      action = "store",
                      type = "int",
                      dest = "step",
                      default = 1,
                      help = """Step for the sliding window. Default is %default.""")

    parser.add_option("-g","--gap",
                      action = "store",
                      type = "int",
                      dest = "gap",
                      default = 50,
                      help = """Size of the gap between the reads. Default is %default.""")

    parser.add_option("-l","--length",
                      action = "store",
                      type = "int",
                      dest = "length",
                      default = 100,
                      help = """Length of the reads. Default is %default.""")

    parser.add_option("--snp",
                      action = "store",
                      type = "float",
                      dest = "snp",
                      default = 0,
                      help = """Percentage of SNPs to be generated. If set to 1 then it will generate a SNP every 100 nucleotides. If set to 0 then no SNPs are generated. Default is %default.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output1_filename and
            options.output2_filename
            ):
        parser.print_help()
        sys.exit(1)


    # running
    fasta2reads(
        options.input_filename,
        options.output1_filename,
        options.output2_filename,
        options.output3_filename,
        options.step,
        options.gap,
        options.length,
        options.snp
        )
    #
