#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gives information regarding the overlap of the mate-reads (i.e. library/fragment size) using two FASTQ files as input.



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
import multiprocessing
import itertools
import string
import gzip
import gc

#ttable = string.maketrans("ACGTYRSWKMBDHV-.","TGCARYSWMKVHDB-.")
ttable = string.maketrans("ACGT","TGCA") # global

def dnaReverseComplement(seq):
    #seq = seq.upper()
    seq = seq.translate(ttable)
    return seq[::-1]


def fastq(file_name, size_read_buffer = 10**8):
    fid = None
    if file_name == '-':
        fid = sys.stdin
    elif file_name.lower().endswith('.gz'):
        fid = gzip.open(file_name,'r')#,buffering=16*1024*1024 + 8)
    else:
        fid = open(file_name,'r')#,buffering=16*1024*1024 + 8) # or buffering = (64*1024)+8 => 64KB ; if the 8 bytes are not added, it will cause poor performance
    piece = [None,None,None,None]
    i = 0
    j = 1

    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for line in lines:
            piece[i] = line
            i = i + 1
            j = j + 1
            if i == 1 and not line.startswith('@'):
                print >>sys.stderr,"ERROR: Fastq file '%s' has an unexpected read id that is '%s' ('@' was expected in the beginning) at line %d!" % (file_name,line.rstrip('\r\n'),j)
                yield ('myexit','1')
                piece = None
                break
            elif i == 3 and not line.startswith('+'):
                print >>sys.stderr,"ERROR: Fastq file '%s' has an unexpected read id that is '%s' ('+' was expected in the beginning)! at line %d!" % (file_name,line.rstrip('\r\n'),j)
                yield ('myexit','1')
                piece = None
                break
            elif i == 4:
                bucket = (piece[0][:-1],
                          piece[1][:-1]
                         )
                yield bucket
                i = 0
                piece = [None,None,None,None]
    if piece and len(piece) != 4:
        print >>sys.stderr,"WARNING: Found unexpected ending of FASTQ file '%s' but still continuing..." % (file_name,)
    fid.close()




def reads_from_paired_fastq_file(file_name_1, file_name_2, nn, fail_gracefully = False, size_read_buffer = 10**8):
    #
    for (pie1,pie2) in itertools.izip(fastq(file_name_1,size_read_buffer = 10**8),
                                      fastq(file_name_2,size_read_buffer = 10**8)):
        # check if the read names are matching (do they form a pair)?
        r1 = pie1[0].partition(" ")[0].rstrip('\r\n')
        r2 = pie2[0].partition(" ")[0].rstrip('\r\n')
        if r1 == 'myexit':
            getout = True
            yield ('myexit','1','myexit','1')
            break
        elif not (r1.startswith('@') and r2.startswith('@')):
            print >>sys.stderr, "ERROR: Something wrong with the id of the reads!"
            print >>sys.stderr, "   - read id from '%s' is '%s'" % (file_name_1,pie1[0])
            if r2 == 'myexit':
                print >>sys.stderr, "   - mate read id from '%s' does not exit (i.e. file ends prematurely)!" % (file_name_2,)
            else:
                print >>sys.stderr, "   - mate read id from '%s' is '%s'" % (file_name_2,pie2[0])
            yield ('myexit','1','myexit','1')
            break

        elif not (r1 == r2 or ((r1[-1] == "1" or r1[-1] == "2") and (r1[-2] == '/' or r1[-2] == '.') and r1[:-2] == r2[:-2]) ):
            print >>sys.stderr, "ERROR: The input FASTQ files do not form a pair!"
            print >>sys.stderr, "   - read id from '%s' is '%s'" % (file_name_1,r1)
            if r2 == 'myexit':
                print >>sys.stderr, "   - mate read id from '%s' does not exit (i.e. file ends prematurely)!" % (file_name_2,)
            else:
                print >>sys.stderr, "   - mate read id from '%s' is '%s'" % (file_name_2,r2)
            yield ('myexit','1','myexit','1')
            break

        elif len(pie1[1]) != nn or len(pie2[1]) != nn:
            print >>sys.stderr, "ERROR: The input reads from the input FASTQ have different lengths compared to the rest of reads (it expects that all input reads have exactly the same lengths)!"
            print >>sys.stderr, "   - read id from '%s' is '%s'" % (file_name_1,r1)
            print >>sys.stderr, "   - mate read id from '%s' is '%s'" % (file_name_2,r2)
            if fail_gracefully:
                yield ('myexit','0','myexit','0')
            else:
                yield ('myexit','1','myexit','1')
            break

        bucket = (pie1[0],
                  pie1[1],
                  pie2[0],
                  pie2[1]
                  )
        yield bucket


class lines_to_file:
    #
    def __init__(self, file_name, size_buffer=10**8):
        self.file_name = file_name
        if file_name:
            self.file_handle = open(file_name,'w')
        self.size_buffer = size_buffer
        self.data = []
        self.size = 0
    #
    def add_line(self,line):
        line = line.rstrip('\r\n')+'\n'
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size = self.size + len(line)
        if self.size > self.size_buffer:
            self.__write_buffer()
    #
    def add_lines(self,lines):
        gc.disable()
        lines = [line.rstrip('\r\n')+'\n' for line in lines]
        self.data.extend(lines)
        gc.enable()
        self.size = self.size + sum([len(line) for line in lines])
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




def fast_alignment5(sa, sb, n, positions, wiggle = 2):
    lib = -1
    xa = sa # ""
    xb = ' '* len(sa) + sb # ""
    mis = -1
    misp = -1
    for (pa,pb) in positions:
        z = sa[pa:pb]
        nz = len(z)
        if (not z) or (pb-pa != nz) or z.find('N') != -1 or z.find('.') != -1 or z[0]*nz == z:
            continue
        #p = sb.find(z,wiggle,-wiggle)
        p = sb.find(z,0,-wiggle)
        if p != -1:
            if pa > p:
                lib = pa + n - p
                pap = pa - p
                t = ' ' * pap
                xa = sa
                xb = "%s%s" % (t,sb)
                lxa = len(xa)
                mis = len([1 for ix in xrange(pap,lxa) if (xa[ix] != xb[ix] or (xa[ix] == 'N' and xb[ix] == 'N'))])
                misp = float(mis)/float(lxa-pap)
            else:
                lib = pa + n - p
                ppa = p - pa
                t = ' ' * ppa
                xa = "%s%s" % (t,sa)
                xb = sb
                lxb = len(xb)
                mis = len([1 for ix in xrange(ppa,lxb) if (xa[ix] != xb[ix] or (xa[ix] == 'N' and xb[ix] == 'N'))])
                misp = float(mis)/float(lxb-ppa)

            #print xa
            #print xb
            #print lib, len("%s%s" % (t,sb))
            #print ""
            break

    if misp != -1 and misp > 0.3:
        # too many mismatches not good alignment
        lib = -1
        xa = sa # ""
        xb = ' '* len(sa) + sb # ""
        mis = -1

    return (lib,xa,xb,mis)

def fast_alignment3(sa, sb, n, positions, wiggle = 2):
    lib = -1
    xa = sa #""
    xb = ' '* len(sa) + sb # ""
    mis = -1
    misp = -1
    for (pa,pb) in positions:
        z = sb[pa:pb]
        nz = len(z)
        if (not z) or (pb-pa != nz) or z.find('N') != -1 or z.find('.') != -1 or z[0]*nz == z:
            continue
        #p = sa.find(z,wiggle,-wiggle)
        p = sa.find(z,0,-wiggle)
        if p != -1:
            if p < pa:
                lib =  n - pa + p
                pap = pa - p
                t = ' ' * pap
                xa = "%s%s" % (t,sa)
                xb = sb
                lxb = len(xb)
                mis = len([1 for ix in xrange(pap,lxb) if (xa[ix] != xb[ix] or (xa[ix] == 'N' and xb[ix] == 'N'))])
                misp = float(mis)/float(lxb-pap)
            else:
                lib =  p + n - pa
                ppa = p - pa
                t = ' ' * (p - pa)
                xa = sa
                xb = "%s%s" % (t,sb)
                lxa = len(xa)
                mis = len([1 for ix in xrange(ppa,lxa) if (xa[ix] != xb[ix] or (xa[ix] == 'N' and xb[ix] == 'N'))])
                misp = float(mis)/float(lxa-ppa)
#            print xa
#            print xb
#            print lib, len("%s%s" % (t,sb))
#            print ""
            break
            
    if misp != -1 and misp > 0.3:
        # too many mismatches not good alignment
        lib = -1
        xa = sa #""
        xb = ' '* len(sa) + sb # ""
        mis = -1
            
    return (lib,xa,xb,mis)


#
#
#
def compute(stuff):
#def compute(mate,o,na,nb):
    mate = stuff[0]
    o = stuff[1].o
    na = stuff[1].na
    nb = stuff[1].nb
    a = mate[1]
    b = dnaReverseComplement(mate[3])
    id1 = mate[0]
    id2 = mate[2]
    mis = -1

    if id1 == 'myexit':
        f = na
        x = a
        y = b
        mis = 0
    else:
        f = -1
        if a == b:
            f = na
            x = a
            y = b
            mis = 0
        else:
            (f,x,y,mis) = fast_alignment5(a, b, na, [(na - o - 1, na - 1), (na - o - 10, na - 10), (na - o - 20, na - 20), (na - o - 30, na - 30)])
            if f == -1:
                (f,x,y,mis) = fast_alignment3(a, b, na, [(nb - o - 1, nb - 1), (nb - o - 10, nb - 10), (nb - o - 20, nb - 20), (nb - o - 30, nb - 30)])

    return (f,x,y,id1,id2,mis)


#
#
#
class param:
    def __init__(self):
        self.o = None
        self.na = None
        self.nb = None

#
#
#
def overlap(input_1_filename,
            input_2_filename,
            output_stat_filename,
            output_alignment_filename = None,
            merged = False,
            size_overlap = 15,
            cpus = 0,
            verbose = True,
            fail_gracefully = False):
			
    o = size_overlap
    print >>sys.stderr,"overlap =", o

    print >>sys.stderr,"Reading..."
    print >>sys.stderr," - ",input_1_filename
    print >>sys.stderr," - ",input_2_filename

    k = 0
    i = 0
    # find fast the length of the read
    d = file(input_1_filename,"r").readlines(50000)
    nax = set([len(el.rstrip('\r\n')) for i,el in enumerate(d) if i%4 == 1])
    na = nax.pop()
    d = file(input_2_filename,"r").readlines(50000)
    nbx = set([len(el.rstrip('\r\n')) for i,el in enumerate(d) if i%4 == 1])
    nb = nbx.pop()

    if len(nax) != 0 or len(nbx)!= 0 or na != nb:
        print >>sys.stderr, "ERROR: The input reads from the input FASTQ have different lengths compared to the rest of reads (it expects that all input reads have exatcly the same lengths)!"
        t = "Analysis skipped because the reads have different lengths in the two input FASTQ files!"
        file(options.output_stat_filename,"w").writelines(t)
        if fail_gracefully:
            sys.exit(0)
        else:
            sys.exit(1)

    library = dict()

    para = param()
    para.o = o
    para.na = na
    para.nb = nb
    nn = na

    #
    if cpus == 0:
        cpus = multiprocessing.cpu_count()
    if verbose:
        print >>sys.stderr,"Using",cpus,"process(es)..."

    pool = multiprocessing.Pool(processes=cpus)

    getout = -1
    if not output_alignment_filename:
        for w in pool.imap_unordered(compute,
                                     itertools.izip_longest(
                                        reads_from_paired_fastq_file(
                                            input_1_filename, 
                                            input_2_filename, 
                                            nn, 
                                            fail_gracefully = fail_gracefully),
                                        [],
                                        fillvalue = para),
                                     chunksize = 100
                                     ):
            f = w[0]
            x = w[1]
            #y = w[2]


            if w[3] == 'myexit':
                getout = int(x)
                break

            if f != -1:
                library[f] = library.get(f,0) + 1
            
            i = i + 1
            if i % 10000000 == 0:
                print >>sys.stderr,"Reading... %d reads" % (i,)
    else: # flag_log
        log = lines_to_file(output_alignment_filename)

        for w in pool.imap_unordered(compute,
                                     itertools.izip_longest(
                                         reads_from_paired_fastq_file(
                                            input_1_filename, 
                                            input_2_filename, 
                                            nn, 
                                            fail_gracefully = fail_gracefully),
                                         [],
                                         fillvalue = para),
                                     chunksize = 100
                                     ):
            f = w[0]
            x = w[1]
            y = w[2]

            if w[3] == 'myexit':
                getout = int(x)
                break

            if f != -1:
                library[f] = library.get(f,0) + 1
            if merged:
                if w[5] == -1:
                    log.add_lines([";"])
                else:
                    log.add_lines(["%s;%s;%s" % (len(x),len(y),str(w[5]))] )
            else:
                log.add_lines([w[3],x,y,w[4],"mismatches = "+str(w[5]) ,"",""]) # read 1 id; read seq 1; read seq 2; read id 2; mismatches

            i = i + 1
            if i % 10000000 == 0:
                print >>sys.stderr,"Reading... %d reads" % (i,)
        log.close()

    pool.close()
    pool.join()


    if getout != -1:
        sys.exit(getout)


    print >>sys.stderr,"Total count mate reads = %i" % (i,)
    print >>sys.stderr,"Read length = %i" % (na,)

    library = sorted(library.items())
    print >>sys.stderr,"Writing the statistics...", output_stat_filename
    # compute the mate reads which overlap
    lib = zip(*library)
    total = 0
    if lib:
        total = sum(lib[1])
    # compute the wasted sequenced nucleotides due to overlapping
    w = 0
    for k,v in library:
        w = w + ( 2 * na - k ) * v
    # most common fragment size
    k_max = 0
    v_max = 0
    for k,v in library:
        if v_max <= v:
            v_max = v
            k_max = k

    print >>sys.stderr,"Count overlapping pair-reads = %i [%.3f%%]" % (total,total/float(i)*100)
    print >>sys.stderr,"Count wasted nucleotides = %i [%.3f%%]" %  (w,100*w/float(2*total*na) if total else 0)
    print >>sys.stderr,"Most common fragment size = %i [%i counts]" %  (k_max,v_max)

    data = []
    data.append("Input FASTQ file 1: %s\n" % (input_1_filename,))
    data.append("Input FASTQ file 2: %s\n" % (input_2_filename,))
    data.append("Length read 1: %i\n" % (na,))
    data.append("Length read 2: %i\n" % (nb,))
    data.append("Minimum overlapping size considered by the program: %i\n" % (o,))
    data.append("Total count pair-reads: %i [%.3f%%]\n" % (i,100.000))
    data.append("Total count reads: %i [%.3f%%]\n" % (2*i,200.000))
    data.append("Count pair-reads overlapping: %i [%.3f%%]\n" % (total,total/float(i)*100))
    data.append("Count total wasted nucleotides (due to overlappings): %i [%.3f%%]\n" % (w,100*w/float(2*i*na)))
    data.append("Most common fragment size amont overlapping pair-reads: %i\n" % (k_max,))
    data.append("Count of the most common fragment size among overlapping pair-reads: %i\n" % (v_max,))

    data.append("\n")
    data.append("Fragment_length\tpair-reads_counts\tpair-reads_percentage\tCumulative_pair-reads_percentage\n")
    da = [ [str(line[0]), str(line[1]), float(100*line[1])/float(i)] for line in library]
    if library:
        x = i - total
        da.append([">%s" % (library[-1][0],),str(x), float(100*x)/float(i)])
    # add the cumulative column
    dax = []
    z = 0
    for line in da:
        z = line[2] + z
        l = (line[0], line[1], "%.3f%%" % (line[2],), "%.3f%%" % (z,))
        dax.append('\t'.join(l)+'\n')
    dax.append("\n")
    data.extend(dax)

    file(output_stat_filename,"w").writelines(data)

######################################################################
######################################################################
######################################################################
######################################################################


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Gives information regarding the overlap of the mate-reads (i.e. library/fragment size) using two FASTQ files as input."""
    version = "%prog 0.14 beta"

    parser = optparse.OptionParser(usage       = usage,
                                   description = description,
                                   version     = version)

    parser.add_option("-1","--input_1",
                      action = "store",
                      type = "string",
                      dest = "input_1_filename",
                      help = """The input FASTQ file containing the reads from 5' fragment end.""")

    parser.add_option("-2","--input_2",
                      action = "store",
                      type = "string",
                      dest = "input_2_filename",
                      help = """The input FASTQ file containing the reads from 3' fragment end.""")

    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_stat_filename",
                      help = """It gives a statistics regarding the amount of mate reads which are overlapping.""")

    parser.add_option("-a","--alignment",
                      action = "store",
                      type = "string",
                      dest = "output_alignment_filename",
                      help = """It outputs also the alignment for each found overlapping.""")

    parser.add_option("-m","--merged",
                      action = "store_true",
                      dest = "merged",
                      default = False,
                      help = """It outputs in the alignment file, the info for merging the paired-end reads (len_read_1,len_read_2,mismatches_in_overlapping_region).""")
                      
    parser.add_option("-s","--fragment-size",
                      action = "store",
                      type = "string",
                      dest = "output_fragment_size_filename",
                      help = """It outputs the fragment size for paired reads which are found to overlap.""")

    parser.add_option("-v","--overlap",
                      action = "store",
                      type = "int",
                      dest = "overlap",
                      default = 13,
                      help = """The minimum length of the region which is considered an overlap. Default is %default.""")

    parser.add_option("-f","--fail-gracefully",
                      action = "store_true",
                      dest = "fail_gracefully",
                      default = False,
                      help = """Fail gracefully in case the the input reads have different lengths. Default is %default.""")

    parser.add_option("-p", "--processes",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 0,
                      help = """Number of parallel processes/CPUs to be used for computations. In case of value 0 then the program will use all the CPUs which are found. The default value is %default.""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.input_1_filename and
            options.input_2_filename and
            options.output_stat_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")

    overlap(options.input_1_filename,
            options.input_2_filename,
            options.output_stat_filename,
            options.output_alignment_filename,
            options.merged,
            size_overlap = options.overlap,
            cpus = options.processes,
            verbose = True, # verbose
            fail_gracefully = options.fail_gracefully
            )
    #
