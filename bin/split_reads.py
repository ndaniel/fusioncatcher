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
#import tempfile


def int2str(x,n=2):
    x = str(x)
    return '0' * int(n - len(x)) + x

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
def split_reads(f_in, f_list, f_out_1, f_out_2, size_buffer = 2*(10**9)):

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
                print "Warning: Not enough free memory (it needed %d)!!! Trying again with a 50% smaller buffer..." % (sb,)
                sb = int(sb / 2)
                if sb < 10000000:
                    print "Error: Not enough free memory (it needed %d)!!! Giving up..." % (sb,)
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
        for line in reads:
            k = line[0]
            if not r.has_key(k):
                r[k] = []
            r[k].append(int(line[2]))
        reads = r
        gc.enable()

        for read in reads_from_fastq_file(f_in):
            v = reads.get(read[0][1:].rstrip('\r\n'),None)
            if len(v) == 1:
                cut = v[0]
                data1.add_line("%s%s\n+\n%s\n" % (read[0],read[1][0:cut+1],read[2][0:cut+1]))
                data2.add_line("%s%s+\n%s" % (read[0],read[1][cut+1:],read[2][cut+1:]))
            else:
                i = 0
                for cut in v:
                    w = read[0][:-1]+'__'+int2str(i)
                    data1.add_line("%s\n%s\n+\n%s\n" % (w,read[1][0:cut+1],read[2][0:cut+1]))
                    data2.add_line("%s\n%s+\n%s" % (w,read[1][cut+1:],read[2][cut+1:]))
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
    version="%prog 0.15 beta"

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


    parser.add_option("--buffer-size",
                      action = "store",
                      type = "int",
                      default = 2*(10**9),
                      dest = "bucket",
                      help="""The size of the buffer used for keeping the list of reads ids (given by --list). Default is %default.""")


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
                options.bucket)
