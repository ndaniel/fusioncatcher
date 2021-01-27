#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It removes the reads which contain short tandem repeats (STR), as shown here:
Gymrek M. et al., "lobSTR: A short tandem repeat profiler for personal genomes",
Genome Res. 2012 22: 1154-1162, doi:10.1101/gr.135780.111



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2021 Daniel Nicorici

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
import math
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
    fid = None
    if file_name == '-':
        fid = sys.stdin
    else:
        fid = open(file_name,'r')
    while True:
        gc.disable()
        lines = fid.readlines(size_buffer)
        gc.enable()
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
                shutil.rmtree(file_name)
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




def counter(sequence,nucleotide = 2):
    freq = {}
    len_sequence = len(sequence) - nucleotide + 1
    for i in xrange(0,len_sequence):
        e = sequence[i:i+nucleotide]
        freq[e] = freq.get(e,0) + 1
    return freq

def plus(a,b):
    ks = set(a.keys()+b.keys())
    d = {}
    for k in ks:
        d[k] = a.get(k,0) + b.get(k,0)
    return d

def minus(a,b):
    ks = set(a.keys()+b.keys())
    d = {}
    for k in ks:
        x = a.get(k,0) - b.get(k,0)
        if x > 0:
            d[k] = x
    return d

def plusminus(a,b,c):
    ks = set(a.keys()+b.keys()+c.keys())
    d = {}
    for k in ks:
        x = a.get(k,0) + b.get(k,0) - c.get(k,0)
        if x > 0:
            d[k] = x
    return d

def bits(d):
    v = d.values()
    n = sum(v)
    if n != 0:
        v = sum([-(float(e)/float(n))*math.log(float(e)/float(n),2) for e in v if e != 0])
    else:
        v = 0
    return v


def codelength(s,w=24,o=12,kmer=2):
    # w = window length
    # o = window overlap
    # n = nucleotides in kmer
#    print "--------------------",s,w,o,kmer
    x = s.upper()
    if x.find('N') !=-1:
        x = x.replace('N','A')
    step = w-o
    c1 = {}
#    print "c1 empty"
    c2 = counter(x[0:o],kmer)
#    print "0->",x[0:o],len(x[0:o]),"c2"
    len_s = len(s) - w
    m = len(s)*100000
    #print "===>",s
    r = [0]
    if len_s != 0:
        r = range(0,len_s,step)
    if r and r[-1] != len_s:
        r.append(len_s)
    for i in r:
#        print "1->",x[i+w-step-kmer+1:i+w],len(x[i+w-step-kmer+1:i+w]),i,'---',i+w-step-kmer+1,i+w,"c3"
        c3 = counter(x[i+w-step-kmer+1:i+w],kmer)
        #c2 = plus(minus(c2,c1),c3)
#        print "c2 before:",c2,"c1",c1
        c2 = plusminus(c2,c3,c1)
#        print "c2: no of kmers",sum(c2.values()),c2
        b = bits(c2)
        if b < m:
            m = b
        c1 = counter(x[i:i+step+kmer-1],kmer)
#        print "2->",x[i:i+step+kmer-1],len(x[i:i+step+kmer-1]),i,'---',i,i+step,"c1"
        #print "       ",b,sorted(c2.items()),'[',i,'-',i+w,']',s[i:i+w+1]
        #print "       ",b,sorted(counter(x[i:i+w+1]).items()),"(validate)"
    #x = (stuff[0][0],stuff[0][1],stuff[0][2],m)
    return m
#def wrap_codelength(stuff):
#    m = codelength(stuff[0][1],stuff[1].window_length,stuff[1].window_overlap)
#    return (stuff[0][0],stuff[0][1],stuff[0][2],m)

#
#
#
def wrap_codelength(stuff):
#def codelength(s,w=24,o=12):
    x = None
    if stuff and type(stuff[0]).__name__ != 'instance':
        s = stuff[0][1]
        w = stuff[1].window_length
        o = stuff[1].window_overlap
        nuc = stuff[1].nucleotide
        #
        m = codelength(s,w,o,nuc)
        x = (stuff[0][0],stuff[0][1],stuff[0][2],m)
    return x


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
def str_filter(
               file_input,
               file_output,
               file_str = None,
               file_log = None,
               window_length = 24,
               window_overlap = 12,
               kmer = 2,
               threshold = 1.4,
               cpus = 0,
               verbose = True):
    """
    'Short Tandem Repeat' Filter
    """


    #
    if cpus == 0:
        cpus = multiprocessing.cpu_count()
    if verbose:
        print >>sys.stderr,"Using",cpus,"process(es)..."

    # in case the pool.imap does not work use itertools.imap instead

    pool = multiprocessing.Pool(processes=cpus)
    p = 0
    f = 0
    fq = tofastq(file_output)
    fs = tofastq(file_str) if file_str else None

    para = param()
    para.window_length = window_length
    para.window_overlap = window_overlap
    para.nucleotide = kmer
    for w in pool.imap_unordered(wrap_codelength,
                       itertools.izip_longest(readfq(fastq(file_input)),
                                              [],
                                              fillvalue = para),
                       chunksize = 100):
        if not w:
            continue
        name = w[0]
        seq = w[1]
        qual = w[2]
        r = w[3]

        p = p + 1
        if r > threshold:
            fq.add(name,seq,qual)
        else:
            f = f + 1
            if fs:
                fs.add(name,seq,qual)
            #print seq,r
    fq.close()
    if fs:
        fs.close()

    pool.close()
    pool.join()


    t = "Empty input file!"
    if p != 0:
        t ="%.5f %% reads removed due to STR (short tandem repeats) content (%d out of %d)!" % ((100*float(f)/float(p)),f,p)
    if verbose:
        print >>sys.stderr,t
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

    description = """It removes the short reads which contain short tandem repeats (STR)."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
(c) 2013-2014 Daniel Nicorici.

"""

    version = "%prog 0.12 beta"

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
                      help = """The output FASTQ file.""")

    parser.add_option("-s","--str",
                      action = "store",
                      type = "string",
                      dest = "str_filename",
                      help = """The output FASTQ file containing the reads which are removed from the input due to high content of short tandem repeats (STR).""")

    parser.add_option("-l","--log",
                      action = "store",
                      type = "string",
                      dest = "log_filename",
                      help = """It outputs a detailed log.""")

    parser.add_option("-v","--overlap",
                      action = "store",
                      type = "int",
                      dest = "window_overlap",
                      default = 12,
                      help = """The length of region where the two consecutive windows are overlapping. Default is %default.""")

    parser.add_option("-n","--length",
                      action = "store",
                      type = "int",
                      dest = "window_length",
                      default = 24,
                      help = """The length of the sliding window. Default is %default.""")

    parser.add_option("-k","--kmer",
                      action = "store",
                      type = "int",
                      dest = "kmer",
                      default = 2,
                      help = """The length of the kmer used in computing the codelength. Default is %default.""")

    parser.add_option("-t","--threshold",
                      action = "store",
                      type = "float",
                      dest = "codelength_threshold",
                      default = 1.4, # original was 2.2
                      help = """Any window which compresses less this threshold is considered to contain a short tandem repeat and the read will be filtered out. Default is %default.""")

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
    str_filter(options.input_filename,
               options.output_filename,
               options.str_filename,
               options.log_filename,
               options.window_length,
               options.window_overlap,
               options.kmer,
               options.codelength_threshold,
               options.processes,
               options.verbose)

if __name__ == '__main__':
    main()
    #
