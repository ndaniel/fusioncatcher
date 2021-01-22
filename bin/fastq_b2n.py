#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes a FASTQ file and converts all nucleotides with low quality score 'B' to ambigous nucleotide N.



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

#
"""
Illumina says:

The Read Segment Quality Control Indicator:

    * At the ends of some reads, quality scores are unreliable. Illumina has an algorithm for identifying these unreliable runs of quality scores, and we use a special indicator to flag these portions of reads
    * A quality score of 2, encoded as a 'B', is used as a special indicator. A quality score of 2 does not imply a specific error rate, but rather implies that the marked region of the read should not be used for downstream analysis.
    * Some reads will end with a run of B (or Q2) basecalls, but there will never be an isolated Q2 basecall.
"""
import os
import sys
import optparse
import string
import gc
import gzip

def reads_from_fastq_file(f_name,size_read_buffer=10**8):
    fid = None
    if f_name == '-':
        fid = sys.stdin
    elif f_name.lower().endswith('.gz'):
        fid = gzip.open(f_name,'r')
    else:
        fid = open(f_name,'r')
    j = 0
    p1 = None
    p2 = None
    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            j = j + 1
            if j == 1:
                p1 = a_line
            elif j == 2:
                p2 = a_line
            elif j == 4:
                yield (p1,p2,a_line)
                p1 = None
                p2 = None
                j = 0
    fid.close()



class lines_to_file:
    def __init__(self,file_name,size_buffer=10**8):
        self.file_name = file_name
        if file_name == '-':
            self.file_handle = sys.stdout
        elif file_name.lower().endswith('.gz'):
            self.file_handle = gzip.open(file_name,'w')
        else:
            self.file_handle = open(file_name,'w')
        self.size_buffer=size_buffer
        self.data=[]
        self.size=0

    def addline(self,line):
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size = self.size + len(line)
        if self.size > self.size_buffer:
            self.__write_buffer()

    def add_line(self,line):
        line=line.rstrip('\r\n')+'\n'
        self.data.append(line)
        self.size=self.size+len(line)
        if self.size>self.size_buffer:
            self.__write_buffer()

    def add_lines(self,lines):
        lines=[line.rstrip('\r\n')+'\n' for line in lines]
        self.data.extend(lines)
        self.size=self.size+sum([len(line) for line in lines])
        if self.size>self.size_buffer:
            self.__write_buffer()
    def __write_buffer(self):
        self.file_handle.writelines(self.data)
        self.size=0
        self.data=[]
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

if __name__=='__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It takes a FASTQ file and converts all nucleotides with low quality score 'B' to ambigous nucleotide N."""
    version="%prog 0.14 beta              Author: Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file (in FASTQ format) containing the short reads to be processed.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTQ file containing the short reads where the nucleotides with low quality score 'B' have been converted to ambigous nucleotide N.""")

    parser.add_option("--find","-f",
                      action="store_true",
                      default = False,
                      dest="find",
                      help="""By default the only tails of B quality are removed. If this is used then the first occurence of B is used to trim the read.""")

    parser.add_option("--ambiguous","-a",
                      action="store_true",
                      default = False,
                      dest="ambiguous",
                      help="""By default the ambigous nucleotides are not convert to As. If this is set then all ambigous nucleotides will be converted to As.""")

    parser.add_option("--sanger","-s",
                      action="store_true",
                      default = False,
                      dest="sanger",
                      help="""By default read qualities are in Illumina v1.5 format and the character 'B' is used to search qualities. If this is used than the character '#' (that is B in Sanger format) is used to search the qualities.""")

    parser.add_option("--replacement","-r",
                      action="store",
                      type="string",
                      dest="replacement",
                      default="N",
                      help="""The character to be use for replacing the nucleotide which have the quality score Q2. Default is '%default'.""")

    parser.add_option("--threshold","-t",
                      action="store",
                      type="int",
                      dest="threshold",
                      default = 0,
                      help="""If if this is larger than zero then all the short reads strictly shorter than this threshold (if the trimming of Bs would be done) are removed. Default is '%default'.""")


    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")

    #print >>sys.stderr,"Starting..."
    data = lines_to_file(options.output_filename)
    threshold = options.threshold
    b = 'B'
    if options.sanger:
        b = '#'
    c = options.replacement
    #ttable = string.maketrans('URYMKWSBDHVN-.TGCA','TACAGACCAAAAAATGCA')
    if options.ambiguous:
        if options.find:
            for reads in reads_from_fastq_file(options.input_filename):
                i = reads[2].find(b) # ORIGINALLY
                if i != -1:
                    n = len(reads[1]) - 1
                    if threshold != 0 and i < threshold:
                        continue
                    r1 = reads[1][:i]+c*(n-i)
                    r2 = reads[2][:i]+b*(n-i)
                    data.addline("%s%s\n+\n%s\n" % (reads[0],r1.upper(),r2))
                else:
                    data.addline("%s%s+\n%s" % (reads[0],reads[1].upper(),reads[2]))
        else:
            for reads in reads_from_fastq_file(options.input_filename):
                i = len(reads[2].rstrip('\r\n').rstrip(b)) # only the trailing B's
                n = len(reads[1].rstrip('\r\n')) - 1
                if i != n:
                    if threshold != 0 and i < threshold:
                        continue
                    r1 = reads[1][:i]+c*(n-i)
                    r2 = reads[2][:i]+b*(n-i)
                    data.addline("%s%s\n+\n%s\n" % (reads[0],r1.upper(),r2))
                else:
                    data.addline("%s%s+\n%s" % (reads[0],reads[1].upper(),reads[2]))
    else:
        if options.find:
            for reads in reads_from_fastq_file(options.input_filename):
                i = reads[2].find(b) # ORIGINALLY
                if i != -1:
                    n = len(reads[1]) - 1
                    if threshold != 0 and i < threshold:
                        continue
                    r1 = reads[1][:i]+c*(n-i)
                    r2 = reads[2][:i]+b*(n-i)
                    data.addline("%s%s\n+\n%s\n" % (reads[0],r1,r2))
                else:
                    data.addline("%s%s+\n%s" % (reads[0],reads[1],reads[2]))
        else:
            for reads in reads_from_fastq_file(options.input_filename):
                i = len(reads[2].rstrip('\r\n').rstrip(b)) # only the trailing B's
                n = len(reads[1].rstrip('\r\n')) - 1
                if i != n:
                    if threshold != 0 and i < threshold:
                        continue
                    r1 = reads[1][:i]+c*(n-i)
                    r2 = reads[2][:i]+b*(n-i)
                    data.addline("%s%s\n+\n%s\n" % (reads[0],r1,r2))
                else:
                    data.addline("%s%s+\n%s" % (reads[0],reads[1],reads[2]))

    #print >>sys.stderr,"Done."
    #
