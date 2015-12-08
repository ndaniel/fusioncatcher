#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given a list of short read names and cut position, it will extracts them from
an input FASTQ file and split each read into two (paired) reads in two separate FASTQ files.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

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
import gc
import shutil
import gzip
import string
#import tempfile


ttable = string.maketrans("ACGTYRSWKMBDHV-.","TGCARYSWMKVHDB-.") # global
#
#
#
def reversecomplement(seq):
    #seq = seq.upper()
    #seq = seq.rstrip('\r\n').translate(ttable)
    return seq.rstrip('\r\n').translate(ttable)[::-1]

#
#
#
def reverse(seq):
    #seq = seq.upper()
    return seq.rstrip('\r\n')[::-1]

#
#
#
def int2str(x,n=2):
    x = str(x)
    return '0' * int(n - len(x)) + x

#
#
#
def reads_from_fastq_file(f_name,size_read_buffer=10**8):
    fid = None
    if f_name == '-':
        fid = sys.stdin
    elif f_name.lower().endswith('.gz'):
        fid = gzip.open(f_name,'r')
    else:
        fid = open(f_name,'r')
    piece = []
    i = 0
    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            i = i + 1
            piece.append(a_line)
            if i == 4:
                yield (piece[0],piece[1],piece[3])
                piece = []
                i = 0
    fid.close()

class lines_to_file:
    def __init__(self,file_name,size_buffer=10**8):
        self.file_name = file_name
        if file_name:
            if file_name == "-":
                self.file_handle = sys.stdout
            elif file_name.lower().endswith('.gz'):
                self.file_handle = gzip.open(file_name,'w')
            else:
                self.file_handle = open(file_name,'w')
        self.size_buffer=size_buffer
        self.data=[]
        self.size=0
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
        lines=[line.rstrip('\r\n')+'\n' for line in lines]
        self.data.extend(lines)
        self.size=self.size+sum([len(line) for line in lines])
        gc.enable()
        if self.size>self.size_buffer:
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

##################
def split_reads(f_in, f_list, f_out_1, f_out_2, wiggle = 0, gap = 0, anchor = 15, replace_solexa_ids = "", rc = False, size_buffer = 2*(10**9)):

    data1 = lines_to_file(f_out_1)
    data2 = lines_to_file(f_out_2)

    fid = open(f_list,'r')
    reads = []
    while True:

        p = fid.tell()
        err = True
        sb = size_buffer
        if sb == 0 :
            sb = 2*(10**9)
        while err:
            gc.disable()
            try:
                lines = fid.readlines(sb)
            except MemoryError:
                print >>sys.stderr,"Warning: Not enough free memory (it needed %d)!!! Trying again with a 50% smaller buffer..." % (sb,)
                sb = int(sb / 2)
                if sb < 10000000:
                    print >>sys.stderr,"Error: Not enough free memory (it needed %d)!!! Giving up..." % (sb,)
                    os.system("free -m")
                    sys.exit(1)
                err = True
                fid.seek(p)
            else:
                err = False
            gc.enable()

        if not lines:
            break

        gc.disable()
        reads = [line.rstrip('\r\n').partition("\t") for line in lines]
        gc.enable()

        gc.disable()
        r = dict()
        wiggle_range = range(-wiggle,wiggle+1)

        for line in reads:
            k = line[0]
            if not r.has_key(k):
                r[k] = set()
            w = int(line[2])
            for wig in wiggle_range:
                r[k].add(w+wig)
        reads = r
        gc.enable()
        am1 = anchor - 1
        am2 = anchor - 2

        for read in reads_from_fastq_file(f_in):
            v = reads.get(read[0][1:].rstrip('\r\n'),None)
            if not v:
                continue
            v = list(v)
            i = 0
            unique = set()
            if gap != 0:
                for agap in xrange(1,gap+1):
                    for cut in v:

                        if cut+1-agap > anchor - 1:
                            k1 = cut+1-agap
                            k2 = cut+1
                            if (k1,k2) not in unique:
                                if replace_solexa_ids:
                                    w = read[0][:-1].replace("/",replace_solexa_ids,1)+'__'+int2str(i)
                                else:
                                    w = read[0][:-1]+'__'+int2str(i)
                                r1a = read[1][0:cut+1-agap]
                                r2a = read[2][0:cut+1-agap]
                                r1b = read[1][cut+1:]
                                r2b = read[2][cut+1:]

                                if len(r1a) > am1 and len(r1b) > am2:
                                    data1.add_line("%sa\n%s\n+\n%s\n" % (w,r1a,r2a))
                                    if rc:
                                        data2.add_line("%sb\n%s\n+\n%s\n" % (w,reversecomplement(r1b),reverse(r2b)))
                                    else:
                                        data2.add_line("%sb\n%s+\n%s" % (w,r1b,r2b))
                                    i = i + 1
                                    unique.add((k1,k2))

                        if len(read[1])-(cut+1+agap) > anchor - 1:
                            k1 = cut+1
                            k2 = cut+1+agap
                            if (k1,k2) not in unique:
                                if replace_solexa_ids:
                                    w = read[0][:-1].replace("/",replace_solexa_ids,1)+'__'+int2str(i)
                                else:
                                    w = read[0][:-1]+'__'+int2str(i)
                                r1a = read[1][0:cut+1]
                                r2a = read[2][0:cut+1]
                                r1b = read[1][cut+1+agap:]
                                r2b = read[2][cut+1+agap:]
                                if len(r1a) > am1 and len(r1b) > am2:
                                    data1.add_line("%sa\n%s\n+\n%s\n" % (w,r1a,r2a))
                                    if rc:
                                        data2.add_line("%sb\n%s\n+\n%s\n" % (w,reversecomplement(r1b),reverse(r2b)))
                                    else:
                                        data2.add_line("%sb\n%s+\n%s" % (w,r1b,r2b))
                                    i = i + 1
                                    unique.add((k1,k2))
            else:
                for cut in v:
                    if replace_solexa_ids:
                        w = read[0][:-1].replace("/",replace_solexa_ids,1)+'__'+int2str(i)
                    else:
                        w = read[0][:-1]+'__'+int2str(i)
                    r1a = read[1][0:cut+1]
                    r2a = read[2][0:cut+1]
                    r1b = read[1][cut+1:]
                    r2b = read[2][cut+1:]
                    if len(r1a) > am1 and len(r1b) > am2:
                        data1.add_line("%sa\n%s\n+\n%s\n" % (w,r1a,r2a))
                        if rc:
                            data2.add_line("%sb\n%s\n+\n%s\n" % (w,reversecomplement(r1b),reverse(r2b)))
                        else:
                            data2.add_line("%sb\n%s+\n%s" % (w,r1b,r2b))
                        i = i + 1
    data1.close()
    data2.close()
    fid.close()
    #


if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""Given a list of short read names and cut position, it will extracts them from
an input FASTQ file and split each read into two (paired) reads in two separate FASTQ files."""
    version="%prog 0.19 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ format (Solexa). Can be given as gzipped file too.""")

    parser.add_option("--list",
                      action="store",
                      type="string",
                      dest="input_list_filename",
                      help="""A text file containing on each line a name of short read which should be extracted from the input FASTQ file and its corresponding cut position (e.g. cut at position N (0-offset); part 1 = [0:N] and part 2 = [N+1:end-read]).""")

    parser.add_option("--output-1",
                      action="store",
                      type="string",
                      dest="output_filename_1",
                      help="""The output FASTQ file where is the first part of the reads (on forward strand).""")

    parser.add_option("--output-2",
                      action="store",
                      type="string",
                      dest="output_filename_2",
                      help="""The output FASTQ file where is the second part of the reads (on forward strand).""")


    parser.add_option("--wiggle-size",
                      action = "store",
                      type = "int",
                      default = 0,
                      dest = "wiggle",
                      help="""The size of the wiggle for the cut. If it is 0 then a read is cut into one paired-reads. If it is 1 then a read is cut into 3 paired-reads. Default is %default.""")

    parser.add_option("--gap-size",
                      action = "store",
                      type = "int",
                      default = 0,
                      dest = "gap",
                      help="""The size of the gap for the cut. Default is %default.""")

    parser.add_option("--anchor-size",
                      action = "store",
                      type = "int",
                      default = 15,
                      dest = "anchor",
                      help="""The minimum size of the anchor (for a mapped read which is splited). Default is %default.""")


    parser.add_option("--replace-solexa-ids",
                      action = "store",
                      type = "string",
                      dest = "replace_solexa_ids",
                      help = """In the reads ids the '/' from '/1' and '/2' will be replaced with the string given here.""")


    parser.add_option("--buffer-size",
                      action = "store",
                      type = "int",
                      default = 2*(10**9),
                      dest = "bucket",
                      help="""The size of the buffer used for keeping the list of reads ids (given by --list). Default is %default.""")

    parser.add_option("--output-2-rc",
                      action = "store_true",
                      default = False,
                      dest = "reverse_complement",
                      help="""The Fastq file specified by '--output-2' will be reverse-complemented. Default is %default.""")

#    parser.add_option("--tmp_dir",
#                  action="store",
#                  type="string",
#                  dest="tmp_dir",
#                  default = None,
#                  help = "The directory which should be used as temporary directory. By default is the OS temporary directory.")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.input_list_filename and
            options.output_filename_1 and
            options.output_filename_2
            ):
        parser.print_help()
        sys.exit(1)

    # running

    split_reads(options.input_filename,
                options.input_list_filename,
                options.output_filename_1,
                options.output_filename_2,
                options.wiggle,
                options.gap,
                options.anchor,
                options.replace_solexa_ids if options.replace_solexa_ids else "",
                options.reverse_complement,
                options.bucket)

#
