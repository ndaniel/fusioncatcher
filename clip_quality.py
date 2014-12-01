#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It clips low quality 3' end of the reads based on a quality threshold.



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
import gc
import shutil
import multiprocessing
import itertools


#
#
#
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] == '@': # fastq header line
                    last = l[:-1] # save this line
                    break
        if not last:
            break
        name = last.partition(" ")[0]
        seqs = []
        last = None
        for l in fp: # read the sequence
            if l[0] == '+':
                break
            seqs.append(l[:-1])
        seq = ''.join(seqs)
        lenseq = len(seq)
        lenq = 0
        seqs = []
        for l in fp: # read the quality
            seqs.append(l[:-1])
            lenq = lenq + len(l) - 1
            if lenq >= lenseq: # have read enough quality
                last = None
                yield name, seq, ''.join(seqs) # yield a fastq record
                break
        if last: # reach EOF before reading enough quality
            yield name, seq, None # yield a fasta record instead
            break

#
#
#
def fastq(file_name, size_buffer = 10**8):
    fid = open(file_name,'r')
    while True:
        lines = fid.readlines(10**8)
        if not lines:
            break
        for line in lines:
            yield line
    fid.close()


#
#
#
class tofastq:
    #
    def __init__(self, file_name, size_buffer = 10**8):
        self.file_name = file_name
        if file_name:
            if os.path.isfile(file_name) or os.path.islink(file_name):
                os.remove(file_name)
            elif os.path.isdir(file_name):
                os.rmtree(file_name)
            self.file_handle = open(file_name,'w')
        self.size_buffer = size_buffer
        self.data = []
        self.size = 0
    #
    def add(self, name, seq, qual):
        line = "%s\n%s\n+\n%s\n" % (name,seq,qual)
        self.data.append(line)
        self.size = self.size + len(line)
        if self.size > self.size_buffer:
            self.__write_buffer()
    #
    def __write_buffer(self):
        self.file_handle.writelines(self.data)
        self.size = 0
        self.data = []
    #
    def close(self):
        if self.file_name:
            if self.data:
                self.__write_buffer()
            self.file_handle.close()
            self.file_name = None

    def __del__(self):
        self.close()

def low(quality,score,window_length):
    """
    A = Q1 in FASTQ-SOLEXA
    B = Q2 in FASTQ-SOLEXA
    C = Q3 in FASTQ-SOLEXA
    D = Q4 in FASTQ-SOLEXA
    E = Q5 in FASTQ-SOLEXA => p = 0.32
    F = Q6 in FASTQ-SOLEXA
    J = Q10 in FASTQ-SOLEXA => p = 0.10
    O = Q15 in FASTQ-SOLEXA
    T = Q20 in FASTQ-SOLEXA
    ^ = Q30 in FASTQ-SOLEXA
    a = Q33 in FASTQ-SOLEXA
    ...
    " = Q1 in FASTQ-SANGER
    # = Q2 in FASTQ-SANGER
    & = Q5 in FASTQ-SANGER
    + = Q10 in FASTQ-SANGER
    5 = Q20 in FASTQ-SANGER
    ? = Q30 in FASTQ-SANGER
    B = Q33 in FASTQ-SANGER
    """
    n = -1
    q = [(1 if e <= score else 0) for e in quality]
    for i in xrange(0,len(quality)-window_length):
        m = sum(q[i:i+window_length])
        if m >= round(float(window_length)/2):
            for j in xrange(i,i+window_length):
                if q[j] == 1:
                    n = j
                    break
            break
    return n

#
#
#
def shred(stuff):
    aread = stuff[0]
    par = stuff[1]
    name = aread[0]
    seq = aread[1]
    qual = aread[2]
    cut = low(qual,par.score,par.window)
    f = False
    if cut != -1:
        if cut == 0:
            cut = 1
        seq = seq[:cut]
        qual = qual[:cut]
        f = True
    return (name,seq,qual,f)

#
#
#
class param:
    def __init__(self):
        self.window = None
        self.score = None

#
#
#
def clip(
         file_input,
         file_output,
         file_log = None,
         window_length = 4,
         quality_score = 5, # 10 for Q10
         basescore = 'solexa',
         cpus = 0,
         verbose = True):
    """
    It clips the low quality nucleotides from 3' end of the reads.
    """

    p = 0
    f = 0

    if quality_score < 0 or quality_score > 100:
        print >> sys.stderr,"Error: Unknown Solexa Quality score!"
        sys.exit(1)
    score = '!'
    if basescore.lower() == 'solexa':
        score = chr(ord('@')+quality_score)
    elif basescore.lower() == 'sanger':
        score = chr(ord('!')+quality_score)

    #
    if cpus == 0:
        cpus = multiprocessing.cpu_count()
    if verbose:
        print "Using",cpus,"process(es)..."

    # in case the pool.imap does not work use itertools.imap instead

    pool = multiprocessing.Pool(processes=cpus)
    fq = tofastq(file_output)
#    for (name,seq,qual) in pool.imap(shred, readfq(fastq(file_input)), chunksize=100):
    p = 0
    f = 0
    para = param()
    para.window = window_length
    para.score = score
    for w in pool.imap(shred,
                       itertools.izip_longest(readfq(fastq(file_input)),
                                              [],
                                              fillvalue = para),
                                              chunksize = 100
                                       ):
        name = w[0]
        seq = w[1]
        qual = w[2]
        c = w[3]

        p = p + 1
        if c:
            f = f + 1
        fq.add(name,seq,qual)
    fq.close()
    t = "Empty input file!"
    if p != 0:
        t ="%.5f %% reads clipped due to low quality (less or equal than Q%d) at 3' end (%d out of %d)!" % ((100*float(f)/float(p)),quality_score,f,p)
    if verbose:
        print t
    if file_log:
        file(file_log,'w').write(t+'\n')


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
    usage = "%prog [options] --input <fastq_file> --output <fastq_file> "

    description = """It clips 3' end of the reads based on a quality threshold."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
(c) 2013-2014 Daniel Nicorici.

"""

    version = "%prog 0.11 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
)



    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input FASTQ file.""")

    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output FASTQ file contained the clipped reads.""")

    parser.add_option("-l","--log",
                      action = "store",
                      type = "string",
                      dest = "log_filename",
                      help = """It outputs a detailed log.""")

    parser.add_option("-n","--length",
                      action = "store",
                      type = "int",
                      dest = "window_length",
                      default = 4,
                      help = """Number of consecutive nucleotides with the quality scores below or equal the given threshold. Default is %default.""")

    parser.add_option("-t","--threshold",
                      action = "store",
                      type = "int",
                      dest = "quality_score",
                      default = 5, # original was 10
                      help = """The quality score below (or equal) the nucleotides are considered low quality and will be trimmed (for example 10 for Q10). Default is %default.""")

    parser.add_option("-s","--score-type",
                      action = "store",
                      type = "string",
                      dest = "score",
                      default = 'solexa', # original was 10
                      help = """The quality score system used. The choices are SOLEXA or SANGER. Default is %default.""")

    parser.add_option("-a","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")

    parser.add_option("-q", "--quiet",
                      action = "store_false",
                      dest = "verbose",
                      default = True,
                      help = "Do not print status messages to console.")

    parser.add_option("-p", "--processes",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 0,
                      help = """Number of parallel processes/CPUs to be used for computations. In case of value 0 then the program will use all the CPUs which are found. The default value is %default.""")

    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    #
    clip(options.input_filename,
         options.output_filename,
         file_log = options.log_filename,
         window_length = options.window_length,
         quality_score = options.quality_score,
         basescore = options.score,
         cpus = options.processes,
         verbose = options.verbose)

if __name__ == '__main__':
    main()
    #
