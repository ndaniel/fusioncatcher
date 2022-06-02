#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It fragments the input paired-end reads in short reads.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2022 Daniel Nicorici

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
import datetime
import optparse
import itertools
import string
import itertools
import re
import gc
import shutil
import errno
import gzip
import math

#ttable = string.maketrans("ACGTYRSWKMBDHV-.","TGCARYSWMKVHDB-N") # global
ttable = string.maketrans("ACGTN","TGCAA") # global
empty_read = ['@N123\n','N\n','+\n','I\n'] # global

#
#
#

#
#
#
def dnaReverseComplement(seq):
    #seq = seq.upper()
    seq = seq.translate(ttable)
    return seq[::-1]

#
#
#
def read_fastq(file_name, size_buffer = 10**8):
    fid = None
    if file_name.lower().endswith('.gz'):
        fid = gzip.open(file_name,'r')
    else:
        fid = open(file_name,'r')
    while True:
        gc.disable()
        lines = fid.readlines(10**8)
        gc.enable()
        if not lines:
            break
        for line in lines:
            if line[:-1]:
                yield line
    fid.close()

#
#
#
def reads_from_paired_fastq_file(file_name_1, file_name_2, size_read_buffer = 10**8):
    i1 = read_fastq(file_name_1, size_buffer = size_read_buffer/2)
    i2 = read_fastq(file_name_2, size_buffer = size_read_buffer/2)
    it = 0
    empty_piece = [(None,None),(None,None),(None,None),(None,None)]
    piece = empty_piece[:]
    for (line_1, line_2) in itertools.izip(i1,i2):
        piece[it] = (line_1,line_2)
        it = it + 1
        if it == 4:
            bucket = (
                piece[0][0], # read 1 - id
                piece[1][0], # read 1 - seq
                piece[3][0], # read 1 - qual

                piece[0][1], # mate read 1 - id
                piece[1][1], # mate read 1 - seq
                piece[3][1]  # mate read 1 - qual
                )
            yield bucket
            piece = empty_piece[:]
            it = 0

    if piece and len(piece) != 4:
        print >>sys.stderr,"WARNING: Found unexpected ending of FASTQ files '%s' and '%s' but still continuing..." % (file_name_1,file_name_2)



#
#
#
def reads_from_single_fastq_file(file_name_1, anchor_size = 20, size_read_buffer = 10**8):
    i1 = read_fastq(file_name_1, size_buffer = size_read_buffer)
    it = 0
    empty_piece = [(None,None),(None,None),(None,None),(None,None)]
    piece = empty_piece[:]
    h = 0
    for line_1 in i1:

        it = it + 1
        if it == 1:
            if line_1[:-1].endswith("/1"):
                piece[it-1] = (line_1,line_1[:-2]+"2\n")
            elif line_1[:-1].endswith("/2"):
                piece[it-1] = (line_1,line_1[:-2]+"1\n")
            else:
                piece[it-1] = (line_1,line_1)
        elif it == 2:
            h = int(math.ceil(float(len(line_1) - 1)/float(2))) + anchor_size
            line_2 = dnaReverseComplement(line_1[:-1])
            piece[it-1] = (line_1[0:h]+"\n",line_2[0:h]+"\n")
        elif it == 4:
            line_2 = line_1[0:-1][::-1]
            piece[it-1] = (line_1[0:h]+"\n",line_2[0:h]+"\n")

            bucket = (
                piece[0][0], # read 1 - id
                piece[1][0], # read 1 - seq
                piece[3][0], # read 1 - qual

                piece[0][1], # mate read 1 - id
                piece[1][1], # mate read 1 - seq
                piece[3][1]  # mate read 1 - qual
                )
            yield bucket
            piece = empty_piece[:]
            it = 0

def trim_tail_n(s,q):
    # trim tails of N from a read
    if s:
        n = len(s)
        if n != 1:
            ts = s.rstrip('N')
            m = len(ts)
            r = n - m
            ts = ts.lstrip('N')
            l = m - len(ts)
            if n != 0 and n!= 1 and l+r>0:
                s = s[l:n-r]
                q = q[l:n-r]
    if not s:
        s = "N"
        q = "I"

    return (s,q)
    
    

#
#
#
def int2str(x,n=2):
    x = str(x)
    return '0' * int(n - len(x)) + x


#
#
#
def fragment_fastq(
        input_file_1,
        input_file_2,
        output_file_1,
        output_file_2,
        log_file = None,
        window_size = 82,
        step_size = 60,
        threshold_size_read = 0,
        anchors = 1,
        wiggle_end = 20,
        skip = 0,
        trim_n = False,
        verbose = False,
            ):
    #
    #
    # finding automatically for adapters
    # - read only the first few millions pairs of reads
    
    threshold = threshold_size_read
    window = window_size
    step = step_size
    
    if threshold == 0:
        threshold = window
    
    limit = 10**5
    

    
    fq1 = open(output_file_1,"w")
    fq2 = None
    if output_file_2 and output_file_2 != '-':
        fq2 = open(output_file_2,"w")
    
    t1 = []
    t2 = []
    zn = 0
    z = []
    
    digits = 2
    limit_digits = 10**digits - 1
    
    get_reads = None
    single = False
    if input_file_2 and input_file_2 == "-":
        get_reads = reads_from_single_fastq_file(input_file_1,anchor_size=window_size-step)
        single = True
    else:
        get_reads = reads_from_paired_fastq_file(input_file_1,input_file_2)
        
        
    for bucket in get_reads:
    
        r1 = bucket[0]
        s1 = bucket[1].rstrip("\r\n")
        q1 = bucket[2].rstrip("\r\n")
        n1 = len(s1)
        
        r2 = bucket[3]
        s2 = bucket[4].rstrip("\r\n")
        q2 = bucket[5].rstrip("\r\n")
        n2 = len(s2)


        if trim_n:
            if n1 != 1 and s1.startswith("N") or s1.endswith("N"):
                (s1,q1) = trim_tail_n(s1,q1)
            if n2 != 1 and s2.startswith("N") or s2.endswith("N"):
                (s2,q2) = trim_tail_n(s2,q2)
            n2 = len(s2)
            n1 = len(s1)

        if n1<skip or n2<skip:
            continue
            

        rr1 = []
        if n1 > threshold:
            if n1 != zn:
                zn = n1
                y1 = range(window,n1,step)
                if (not y1) or single or n1 - y1[-1] >= wiggle_end: #y1[-1] != n1:
                    y1.append(n1)
                x1 = [i-window if i-window > -1 else 0 for i in y1]
                z = zip(x1,y1)
                
            rr1 = [(s1[i:j],q1[i:j]) for i,j in z]
        else:
            rr1 = [(s1,q1)]

        rr2 = []
        if n2 > threshold:
            if n2 != zn:
                zn = n2
                y2 = range(window,n2,step)
                if (not y2) or single or n2 - y2[-1] >= wiggle_end: #y2[-1] != n2:
                    y2.append(n2)
                x2 = [i-window if i-window > -1 else 0 for i in y2]
                z = zip(x2,y2)
            rr2 = [(s2[i:j],q2[i:j]) for i,j in z]
        else:
            rr2 = [(s2,q2)]

        if trim_n:
            rr1 = [trim_tail_n(rra,rrb) for rra,rrb in rr1]
            rr2 = [trim_tail_n(rra,rrb) for rra,rrb in rr2]


        i = -1
        u = set()
        for j in xrange(min(len(rr2),anchors)):
            as2 = rr2[j][0]
            aq2 = rr2[j][1]
            k = -1
            for (as1,aq1) in rr1:
                k = k + 1
                u.add((k,j))
                if len(as1) < skip or len(as2) < skip:
                    continue

                i = i + 1
                if i > limit_digits:
                    digits = digits + 1
                    limit_digits = 10**digits - 1
                ids = int2str(i,digits)

                if fq2:
                    t1.append("@%s_%s%s\n+\n%s\n" % (ids,r1[1:],as1,aq1))
                    t2.append("@%s_%s%s\n+\n%s\n" % (ids,r2[1:],as2,aq2))
                else:
                    t1.append("@%s_%s%s\n+\n%s\n" % (ids,r1[1:],as1,aq1))
                    t1.append("@%s_%s%s\n+\n%s\n" % (ids,r2[1:],as2,aq2))


        for j in xrange(min(len(rr1),anchors)):
            as1 = rr1[j][0]
            aq1 = rr1[j][1]
            #rr2.pop(0) # this is already done
            k = -1
            for (as2,aq2) in rr2:
                k = k + 1
                if (j,k) in u:
                    continue

                if len(as1) < skip or len(as2) < skip:
                    continue

                i = i + 1
                if i > limit_digits:
                    digits = digits + 1
                    limit_digits = 10**digits - 1
                ids = int2str(i,digits)

                if fq2:
                    t1.append("@%s_%s%s\n+\n%s\n" % (ids,r1[1:],as1,aq1))
                    t2.append("@%s_%s%s\n+\n%s\n" % (ids,r2[1:],as2,aq2))
                else:
                    t1.append("@%s_%s%s\n+\n%s\n" % (ids,r1[1:],as1,aq1))
                    t1.append("@%s_%s%s\n+\n%s\n" % (ids,r2[1:],as2,aq2))
        
        if len(t1) > limit:
            fq1.writelines(t1)
            t1 = []
        if fq2 and len(t2) > limit:
            fq2.writelines(t2)
            t2 = []

    if t1:
        fq1.writelines(t1)
    fq1.close()

    if fq2:
        if t2:
            fq2.writelines(t2)
        fq2.close()


######################################################################
######################################################################
######################################################################
######################################################################



def main():

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options] --input_1 <fastq_file_read_1> --input_2 <fastq_file_read_2> --output_2 <fastq_file_1> --output_2 <fastq_file_2> "

    description = """It fragments (i.e. breaks up) the input single-end/paired-end reads in reads which are shorter."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2022 Daniel Nicorici

"""

    version = "%prog 0.97 beta"

    parser = MyOptionParser(
        usage       = usage,
        epilog      = epilog,
        description = description,
        version     = version
        )



    parser.add_option("-1","--input_1",
                      action = "store",
                      type = "string",
                      dest = "input_1_filename",
                      help = """The input FASTQ file containing the reads from 5' fragment end (i.e. 5'-3' orientation for read 1 and 3'-5' for read 2 which needs to be reversed-complemented).""")

    parser.add_option("-2","--input_2",
                      action = "store",
                      type = "string",
                      dest = "input_2_filename",
                      default = "-",
                      help = """The input FASTQ file containing the reads from 3' fragment end (i.e. 5'-3' orientation for read 1 and 3'-5' for read 2 which needs to be reversed-complemented). If this is not used then the input FASTQ file is assume to contain single-end reads.""")

    parser.add_option("-f","--output_1",
                      action = "store",
                      type = "string",
                      dest = "output_1_filename",
                      help = """The output FASTQ file where the reads are trimmed.""")

    parser.add_option("-r","--output_2",
                      action = "store",
                      type = "string",
                      dest = "output_2_filename",
                      default = "-",
                      help = """The output FASTQ file where the reads are trimmed. If this is not used then all the output reads are written in '-f' as interleaved.""")

    parser.add_option("-l","--log",
                      action = "store",
                      type = "string",
                      dest = "log_filename",
                      help = """It outputs a detalied statistics of the trimming.""")

    parser.add_option("-w","--window-size",
                      action = "store",
                      type = "int",
                      dest = "window_size",
                      default = 82,
                      help = """The size of the new reads. This should be shorter than the size of input reads. Default is %default.""")

    parser.add_option("-s","--step-size",
                      action = "store",
                      type = "int",
                      dest = "step_size",
                      default = 60,
                      help = """The size of step for the sliding window used to fragment the input reads. This should be shorter than the size of input reads. Default is %default.""")

    parser.add_option("-t","--threshold-read",
                      action = "store",
                      type = "int",
                      dest = "threshold_size_read",
                      default = 0,
                      help = """The reads shorter than this will not be fragmented. If it is set then it should be higher than window size. By default reads shorter than the window size will not be fragmented.""")

    parser.add_option("-a","--anchors",
                      action = "store",
                      type = "int",
                      dest = "anchors",
                      default = 4,
                      help = """Number of anchors. Default is %default.""")

    parser.add_option("-k","--skip-short",
                      action = "store",
                      type = "int",
                      dest = "skip_short",
                      default = 0,
                      help = """The paired-end where both reads are shorter than this will be filtered out. Default is '%default'.""")

    parser.add_option("-e","--wiggle-end",
                      action = "store",
                      type = "int",
                      dest = "wiggle_end",
                      default = 16,
                      help = """The last fragment will not be generated if it overlaps with the previous generated fragment and the non-overlapping segment is strictly shorter than this threshold. This applies only for paired-end reads. Default is '%default'.""")

    parser.add_option("-n","--trim-n",
                      action = "store_true",
                      dest = "trim_n",
                      default = False,
                      help = """If it is set then Ns from both ends of the read are trimmed.""")

    parser.add_option("-q", "--quiet",
                      action = "store_false",
                      dest = "verbose",
                      default = True,
                      help = "Do not print status messages to stdout.")

    parser.add_option("-u","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.input_1_filename and
            options.input_2_filename and
            options.output_1_filename and
            options.output_2_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    #
    fragment_fastq(
        options.input_1_filename,
        options.input_2_filename,
        options.output_1_filename,
        options.output_2_filename,
        options.log_filename,
        options.window_size,
        options.step_size,
        options.threshold_size_read,
        options.anchors,
        options.wiggle_end,
        options.skip_short,
        options.trim_n,
        options.verbose
        )


if __name__ == '__main__':
    main()
    
    
    
#
