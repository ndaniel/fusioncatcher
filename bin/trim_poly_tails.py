#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It trims the leading or trailing poly-A/C/G/T tails of read sequences from a input FASTQ file.


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
    piece = [None,None,None,None]
    i = 0
    while True:
        gc.disable()
        lines=fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            i = i + 1
            piece[i-1] = a_line.rstrip('\r\n')
            #piece.append(a_line[:].rstrip('\r\n'))
            if i == 4:
                piece[2] = '+'
                yield piece
                piece = [None,None,None,None]
                i = 0
    fid.close()

class lines_to_file:
    def __init__(self,file_name,size_buffer=10**8):
        self.file_name=file_name
        if file_name == '-':
            self.file_handle = sys.stdout
        elif file_name.lower().endswith('.gz'):
            self.file_handle = gzip.open(file_name,'w')
        else:
            self.file_handle = open(file_name,'w')
        self.size_buffer=size_buffer
        self.data=[]
        self.size=0
    def add_line(self,line):
        line=line+'\n'
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size=self.size+len(line)
        if self.size>self.size_buffer:
            self.__write_buffer()
    def add_simple_line(self,line):
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size=self.size+len(line)
        if self.size>self.size_buffer:
            self.__write_buffer()
    def add_lines(self,lines):
        gc.disable()
        lines=[line+'\n' for line in lines]
        self.data.extend(lines)
        gc.enable()
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

def trim_poly_5_end(r, q, nucleotide, no_repeats = 9,empty=False):
    n = len(r)
    start = -1
    f = False
    if n > 1 and r[0] == nucleotide and r[1] == nucleotide:
        start = -1
        s = r.upper()
        for i in range(1,n):
            if s[i] == s[i-1]:
                start = i
            else:
                break
        if (n - start + 1) >= no_repeats:
            r = r[start+1:]
            q = q[start+1:]
            f = True # trimmed
            if (not empty) and (not r):
                # do not give reads of length zero!
                r = "A"
                q = "I"
    return (r,q,f)

def trim_poly_3_end(r, q, nucleotide, no_repeats = 9,empty=False):
    n = len(r)
    f = False
    if n > 1 and r[-1] == nucleotide and r[-2] == nucleotide:
        end = n
        s = r.upper()
        for i in range(n-2,-1,-1):
            if s[i] == s[i+1]:
                end = i
            else:
                break
        if (n - end + 1) >= no_repeats:
            r = r[:end]
            q = q[:end]
            f = True # trimmed
            if (not empty) and (not r):
                # do not give reads of length zero!
                r = "A"
                q = "I"
    return (r,q,f)


if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It trims the leading or trailing poly-A/C/G/T tails of read sequences from a input FASTQ file."""
    version="%prog 0.11 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ format.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTQ file containing all the trimmed sequences.""")

    parser.add_option("--repeats","-r",
                      action="store",
                      type="int",
                      dest="repeats",
                      default=9,
                      help="""The number of times a nucleotide specified with '--nucleotide' should be repeated in order to be considered a poly Default is %default.""")

    parser.add_option("--skip_reads","-s",
                      dest="skip_reads",
                      action="store_true",
                      default=False,
                      help="""If this is specified then the reads which are having poly tails are filtered out (i.e. not written to the output) instead of trimming. Default is %default.""")

    parser.add_option("--keep-too-short","-k",
                      dest="keep_too_short",
                      action="store_true",
                      default=False,
                      help="""If this is specified then the reads which are less than N bp will be kept, where N is set using '--keep-too-short-length'. Default is %default.""")

    parser.add_option("--keep-too-short-length","-l",
                      dest = "keep_too_short_length",
                      action = "store",
                      type = "int",
                      default = 20,
                      help = """The threshold used to decide when a read is too short. Default is %default.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    #print >>sys.stderr,"Starting..."
    poly = {}
    poly['A'] = 'A' * options.repeats
    poly['T'] = 'T' * options.repeats
    poly['C'] = 'C' * options.repeats
    poly['G'] = 'G' * options.repeats
    poly['N'] = 'N' * options.repeats
    poly['R'] = 'R' * options.repeats
    poly['Y'] = 'Y' * options.repeats
    poly['S'] = 'S' * options.repeats
    poly['W'] = 'W' * options.repeats
    poly['K'] = 'K' * options.repeats
    poly['M'] = 'M' * options.repeats
    poly['B'] = 'B' * options.repeats
    poly['D'] = 'D' * options.repeats
    poly['H'] = 'H' * options.repeats
    poly['V'] = 'V' * options.repeats

    thr = options.keep_too_short_length

    poly_keys = sorted(poly.keys())
    data = lines_to_file(options.output_filename)
    c = 0
    i = 0
    j = 0

    if options.skip_reads:
        for reads in reads_from_fastq_file(options.input_filename):
            i = i + 1
            id = reads[0]
            ss = reads[1].upper()
            qq = reads[3]

            ts = ss
            tq = qq
            h = False
            caracter1 = ss[-1:]
            caracter2 = ss[0:1]
            if caracter2 and ss.startswith(poly.get(caracter2,'---')):
                h = True
            elif caracter1 and ss.endswith(poly.get(caracter1,'---')):
                h = True
#            for k in poly_keys:
#                if ss.startswith(poly[k]) or ss.endswith(poly[k]):
#                    h = True
#                    break
            if h:
                c = c + 1
            else:
                #data.add_lines([id,ts,'+',tq])
                data.add_simple_line("%s\n%s\n+\n%s\n" % (id,ts,tq))
                j = j + 1
    else:
        for reads in reads_from_fastq_file(options.input_filename):
            i = i + 1
            id = reads[0]
            ss = reads[1].upper()
            qq = reads[3]

            ts = ss
            tq = qq
            h = False
            h1 = False
            h2 = False
            caracter1 = ss[-1:]
            caracter2 = ss[0:1]
            if caracter2 and ss.startswith(poly.get(caracter2,'---')):
                (ts,tq,h1) = trim_poly_5_end(ts, tq, caracter2, no_repeats = options.repeats)
            if caracter1 and ss.endswith(poly.get(caracter1,'---')):
                (ts,tq,h2) = trim_poly_3_end(ts, tq, caracter1, no_repeats = options.repeats)
            if h1 or h2:
                h = True

#            for k in poly_keys:
#                h1 = False
#                h2 = False
#                if ss.startswith(poly[k]):
#                    (ts,tq,h1) = trim_poly_5_end(ts, tq, k, no_repeats = options.repeats)
#                if ss.endswith(poly[k]):
#                    (ts,tq,h2) = trim_poly_3_end(ts, tq, k, no_repeats = options.repeats)
#                if h1 or h2:
#                    h = True
#                    break
            if options.keep_too_short or len(ts) >= thr:
                #data.add_lines([id,ts,'+',tq])
                data.add_simple_line("%s\n%s\n+\n%s\n" % (id,ts,tq))
                j = j + 1
            if h:
                c = c + 1
        data.close()

    print >>sys.stderr,"Found %d reads with poly-A/C/G/T/N tails (equal or more %s repeat nucleotides)" % (c,options.repeats)
    print >>sys.stderr,"Total number of input reads = %d" % (i,)
    print >>sys.stderr,"Total number of reads written in the output = %d" % (j,)
    #
