#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It trims (finds and/or removes by trimming) the paired-reads (FASTQ file
produced by Illumina Solexa) which overlap and contain the adapter. The adapter
is found automatically. Also the partial overlapping between a short read and
the adapter sequence is handled.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2017 Daniel Nicorici

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
import itertools
import re
import gc
import shutil
import errno
import gzip

ttable = string.maketrans("ACGTYRSWKMBDHV-.","TGCARYSWMKVHDB-.") # global
empty_read = ['@N123\n','N\n','+\n','I\n'] # global

#
#
#
def remove_file(a_file):
    if os.path.isfile(a_file) or os.path.islink(a_file):
        os.remove(a_file)
    elif os.path.isdir(a_file):
        shutil.rmtree(a_file)
#
#
#
def linkit(file_input, file_output, kind ='soft'):
    #
    remove_file(file_output)
    if os.path.islink(file_input):
        linkto = os.readlink(file_input)
        if kind == 'soft':
            os.symlink(linkto, file_output)
        elif kind == 'hard':
            try:
                os.link(linkto, file_output)
            except OSError as er:
                print >>sys.stderr,"WARNING: Cannot do hard links ('%s' and '%s')!" % (linkto,file_output)
                shutil.copyfile(linkto, file_output)
#                if er.errno == errno.EXDEV:
#                    # they are on different partitions
#                    # [Errno 18] Invalid cross-device link
#                    shutil.copyfile(linkto, file_output)
#                else:
#                    print >>sys.stderr,"ERROR: Cannot do hard links ('%s' and '%s')!" % (linkto, file_output)
#                    print >>sys.stderr,er
#                    sys.exit(1)

    else:
        if kind == 'soft':
            os.symlink(file_input, file_output)
        elif kind == 'hard':
            try:
                os.link(file_input, file_output)
            except OSError as er:
                print >>sys.stderr,"WARNING: Cannot do hard links ('%s' and '%s')!" % (linkto,file_output)
                shutil.copyfile(linkto, file_output)
#                if er.errno == errno.EXDEV:
#                    # they are on different partitions
#                    # [Errno 18] Invalid cross-device link
#                    shutil.copyfile(file_input, file_output)
#                else:
#                    print >>sys.stderr,"ERROR: Cannot do hard links ('%s' and '%s')!" % (file_input, file_output)
#                    print >>sys.stderr,er
#                    sys.exit(1)


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
def read_first_fastq(file_name, first = 2000000, size_buffer = 10**8):
    fid = None
    if file_name.lower().endswith('.gz'):
        fid = gzip.open(file_name,'r')
    else:
        fid = open(file_name,'r')
    i = 0
    while True:
        gc.disable()
        lines = fid.readlines(10**8)
        gc.enable()
        if (not lines) or i > first:
            break
        for line in lines:
            if line[:-1]:
                i = i + 1
                if i > first:
                    break
                else:
                    yield line
    fid.close()

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
def first_reads_from_paired_fastq_file(file_name_1, file_name_2, first = 2000000, size_read_buffer = 10**8):
    # it provides only the first N mate-reads
    # if first == 0 then all mate-reads are provided
    i1 = read_first_fastq(file_name_1, first = first*4+1, size_buffer = size_read_buffer / 2)
    i2 = read_first_fastq(file_name_2, first = first*4+1, size_buffer = size_read_buffer / 2)
    piece = [(None,None),(None,None),(None,None),(None,None)]
    bucket = []
    it = 0

    global empty_read

    #cut = 2000000 # first number of paired reads which are read
    i = 0
    x = True # get a quality score for the empty read
    for (line_1, line_2) in itertools.izip(i1,i2):
        piece[i] = (line_1,line_2)
        i = i + 1
        if first != 0 and it > first - 1:
            continue
        if i == 4:
            bucket = [piece[0][0],
                      piece[1][0].rstrip('\r\n'),
                      "+\n",
                      piece[3][0],

                      piece[0][1],
                      piece[1][1].rstrip('\r\n'),
                      "+\n",
                      piece[3][1]
                      ]
            yield bucket

            if x:
                empty_read[3] = piece[3][0][0]+'\n'
                x = False

            piece = [(None,None),(None,None),(None,None),(None,None)]
            i = 0

            it = it + 1


#
#
#
def reads_from_paired_fastq_file(file_name_1, file_name_2, size_read_buffer = 10**8):
    i1 = read_fastq(file_name_1, size_buffer = size_read_buffer/2)
    i2 = read_fastq(file_name_2, size_buffer = size_read_buffer/2)
    it = 0
    piece = [(None,None),(None,None),(None,None),(None,None)]
    for (line_1, line_2) in itertools.izip(i1,i2):
        it = it + 1
        piece[it-1] = (line_1,line_2)
        #piece.append((line_1,line_2))
        if it == 4:
            bucket = [piece[0][0],
                      piece[1][0].rstrip('\r\n'),
                      "+\n",
                      piece[3][0],

                      piece[0][1],
                      piece[1][1].rstrip('\r\n'),
                      "+\n",
                      piece[3][1]
                      ]
            yield bucket
            piece = [(None,None),(None,None),(None,None),(None,None)]
            it = 0

    if piece and len(piece) != 4:
        print >>sys.stderr,"WARNING: Found unexpected ending of FASTQ files '%s' and '%s' but still continuing..." % (file_name_1,file_name_2)
#
#
#
class lines_to_file:
    #
    def __init__(self, file_name, size_buffer = 10**8):
        self.file_name = file_name
        if file_name:
            remove_file(file_name)
            self.file_handle = open(file_name,'w')
        self.size_buffer = size_buffer
        self.data = []
        self.size = 0
    #
    def add_line(self, line):
#        line = line.rstrip('\n')+'\n'
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size = self.size + len(line)
        if self.size > self.size_buffer:
            self.__write_buffer()
    #
    def add_lines(self,lines):
#        lines=[line.rstrip('\r\n')+'\n' for line in lines]
        gc.disable()
        self.data.extend(lines)
        gc.enable()
        self.size = self.size + sum([len(line) for line in lines])
        if self.size > self.size_buffer:
            self.__write_buffer()
    #
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

#
#
#
def fast_alignment_adapter(sa, sb, len_adapter = 13, overlap = 13):
    lib = -1
    na = len(sa)
    nb = len(sb)
    adapter5 = ''
    adapter3 = ''
    positions = [ (nb - overlap - x, nb - x) for x in range(3,nb/4,7) ]
    #
    for (pa,pb) in positions:
        if pb - pa < overlap - 1 or pb < 0 or pa < 0:
            break
        p = sa.find(sb[pa:pb], 2)
        if p != -1 and p < pa:
            lib =  p + nb - pa
            if lib > na:
                lib = na
            s = pa - p
            #t = ' ' * s
            #xa = "%s%s" % (t,sa)
            #xb = sb
            #print xa
            #print xb
            # count the mismatches
            mis = len([1 for ix in xrange(lib) if s+ix<nb and (sa[ix] != sb[s+ix] or (sa[ix] == 'N' and sb[s+ix] == 'N') )])
            n_notn = len([1 for ix in xrange(lib) if s+ix<nb and (sa[ix] == 'N' or sb[s+ix] == 'N')])
            #print "mismatches",mis
            if lib > 0 and float(n_notn) / float(lib) > 0.3:
                mis = mis + n_notn
            if (mis > 0 and lib > 0 and float(mis) / float(lib) < 0.2 and lib != na):
                continue
            #print xa
            #print xb
            #print lib, len("%s%s" % (t,sb)),mis
            #print sa[lib:lib+len_adapter]
            #print sb[s-len_adapter:s]
            #print ""
            adapter5 = sa[lib:lib+len_adapter]
            adapter3 = sb[s-len_adapter:s]
            if len(adapter5) < len_adapter:
                adapter5 = ""
            if len(adapter3) < len_adapter:
                adapter3 = ""
            #print adapter5,adapter3
            break
    return (adapter5, adapter3)

#
#
#
def find_hard(ss,adp):
    # try hard finding the adapter with one mismatch
    p = -1
    if adp:
        adp_re = re.compile('|'.join(adp[:i] + '.' + adp[i+1:] for i in xrange(len(adp))))
        for elx in adp_re.finditer(ss):
            p = elx.start()
            break
    return p

#
#
#
def fast_alignment(sa, sb, overlap = 13, wiggle = 2, adpt5 = "", adpt3 = ""):
    # align a read on top of its mate read if this is possible
    na = len(sa)
    nb = len(sb)
    trim_a = na
    trim_b = 0
    len_adpt5 = 0
    if adpt5:
        len_adpt5 = len(adpt5)
    len_adpt3 = 0
    if adpt3:
        len_adpt3 = len(adpt3)
    positions = [ (nb - overlap - x, nb - x) for x in range(3,nb/4,7) ]
    positions.insert(0,[-1,-1]) # this is I need it for allowing first to look for adapter
    # do not forget! check that the positions are all positive
    cut_mis_adapt = 0.3
    cut_mis = 0.3
    p = -1
    p5 = -1
    p3 = -1
    lib = 0
    s = -1
    common = 0
    mis = -1
    mis5 = -1
    mis3 = -1
    first = True
    for (pa, pb) in positions:
        if first:
            first = False
            p5 = -1
            if adpt5:
                p5 = sa.find(adpt5)
            p3 = -1
            if adpt3:
                p3 = sb.find(adpt3)
            if p3 == -1:
                if p5 == -1:
                    continue
                else:
                    mis5 = 0
                    # try harder
                    p3 = find_hard(sb,adpt3)
                    if p3 == -1:
                        continue
                    else:
                        mis3 = 1
            elif p5 == -1:
                mis3 = 0
                # try harder
                p5 = find_hard(sa,adpt5)
                if p5 == -1:
                    continue
                else:
                    mis5 = 1
            #do
            lib = p5 # position on first read
            s = p3 + len_adpt3 # position on second read (reverse-complemented)
            trim_a = lib
            trim_b = s
            common = -1
        elif pb - pa < overlap - 1 or pa < 0 or pb < 0:
            break
        else:
            p = sa.find(sb[pa:pb], wiggle)
            if p == -1 or p > pa:
                continue
            else:
                lib =  p + nb - pa
                s = pa - p
                if lib > na:
                    lib = na

        #if lib == 0:
        #    print "lib=0"
        #print "lib",lib
        #print sa,p
        #print sb,s
        #    t = ' ' * s
        #    xa = "%s%s" % (t,sa)
        #    xb = sb
        #    print xa
        #    print xb
        # count the mismatches in the overlap; N is not considered a mismatch in the overlapping part
        xlib = range(lib)
        mis = len([1 for ix in xlib if s+ix<nb and ((sa[ix] != sb[s+ix] and sa[ix] != 'N' and sb[s+ix] != 'N') or ( sa[ix] == 'N' and sb[s+ix] == 'N'))])
        n_notn = len([1 for ix in xlib if s+ix<nb and (sa[ix] == 'N' or sb[s+ix] == 'N')])
        if lib > 0 and float(n_notn) / float(lib) > 0.3:
            mis = mis + n_notn
        if mis > 0 and lib > 0 and ((mis / float(lib) > cut_mis) or (lib == na and mis / float(lib) > 0.05)):
            continue


        mis5 = -1
        if adpt5:
            mis5 = len([1 for ix in xrange(len_adpt5) if lib+ix<na and adpt5[ix] != sa[lib+ix] and sa[lib+ix] != 'N' and sa[lib+ix] != '.'])
            if float(float(mis5) / float(len_adpt5)) <= float(cut_mis_adapt):
                # trim the read
                trim_a = lib

        mis3 = -1
        if adpt3:
            mis3 = len([1 for ix in xrange(len_adpt3) if s-len_adpt3 + ix>-1 and adpt3[ix] != sb[s-len_adpt3 + ix] and sb[s-len_adpt3 + ix] != 'N' and sb[s-len_adpt3 + ix] != '.'])
            if float(mis3) / float(len_adpt3) <= float(cut_mis_adapt):
                # trim the read
                trim_b = s

        #print xa
        #print xb
        #print lib, len("%s%s" % (t,sb)), mis, mis5, mis3
        #print sa[lib:lib+len_adpt5]
        #print sb[s-len_adpt3:s] if s-len_adpt3 > -1 else sb[0:s]
        #print sa[:trim_a]
        #print sb[trim_b:]
        #print ""
        #raw_input()
        common = lib
        break

    if trim_a == na and trim_b == 0 and common != na:
        common = -1
        trim_a = p5
        trim_b = -1
        if p3 != -1:
            trim_b = p3 + len_adpt3
    #print ":::",trim_a,trim_b
    return (trim_a, trim_b, common, mis, mis5, mis3)

#
#
#
def fix_N_in_overlap(ya, ybr, yb, la, lb, shift = -1):
    # replaces N with its corresponding nucleotide from the mate read
    # if the mate read overlaps over N
    fixed_Ns = 0
    if shift != -1 and (ya.find('N') != -1 or ybr.find('N') != -1):
        #print "-----------------------"
        #print la,lb,shift
        #print ' '*shift+ya
        #print ybr
        e = lb - shift
        if e > la:
            e = la
        f = lb
        if la + shift < lb:
            f = la + shift
        ya1 = ya[0:e]
        ya2 = ya[e:]
        yb1 = ybr[0:shift]
        yb2 = ybr[shift:f]
        yb3 = ybr[f:lb]
        c = []
        #print "before"
        #print ya1+'-'+ya2
        #print yb2+'-'+yb3
        for (ia,ib) in zip(ya1,yb2):
            if ia == 'N' and ib != 'N':
                c.append([ib,ib])
                fixed_Ns = fixed_Ns + 1
            elif ia != 'N'and ib == 'N':
                c.append([ia,ia])
                fixed_Ns = fixed_Ns + 1
            else:
                c.append([ia,ib])
        c = zip(*c)
        ya1 = ''.join(c[0])
        yb2 = ''.join(c[1])
        #print "after"
        #print ya1
        #print yb2
        #print "final"
        #print ya1 + ya2
        #print dnaReverseComplement(yb1+yb2+yb3)
        #print "-----------------------"
        #
        return (ya1 + ya2, dnaReverseComplement(yb1 + yb2 + yb3), fixed_Ns)
    else:
        return (ya, yb, fixed_Ns)


#
#
#
def compute(stuff):
#def compute(mate, reads_overlap, wiggle, adapter5, adapter3, flag_log):
    mate = stuff[0]
    w = stuff[1]
    reads_overlap = w.reads_overlap
    wiggle = w.wiggle
    adapter5 = w.adapter5
    adapter3 = w.adapter3
    flag_log = w.flag_log

    a = mate[1]
    b = dnaReverseComplement(mate[5])
    bb = mate[5]
    na = len(a)
    nb = len(b)
    (qa, qb, com, ml, m5, m3) = fast_alignment(a, b, overlap = reads_overlap, wiggle = wiggle, adpt5 = adapter5, adpt3 = adapter3)

    st1 = -1
    st2 = -1
    stn = -1
    jj = 0
    #print qa,qb,na,nb
    # trimming
    fixedns = 0
    if com > -1:
        if com == na: # full/perfect overlap
            jj = jj + 1
            (mate[1], mate[5], fixed_Ns) = fix_N_in_overlap(a, b, bb, na, nb, nb - com)
            fixedns = fixedns + fixed_Ns
            #mate[5] = 'N' * nb
            #mate[5] = "N" ## orig
            #mate[7] = mate[7][0] + '\n' ## orig
            mate[5] = dnaReverseComplement(mate[1])
            mate[7] = mate[3][:-1][::-1] + '\n'
            stn = nb
        elif qa != na: # partial overlap
            if qb != 0:
                jj = jj + 2
                if qa == 0: # both have the adapter and no overlap
                    #mate[1] = 'N'
                    #mate[5] = 'N'
                    mate[1] = 'N'
                    mate[3] = mate[3][0] + '\n'
                    mate[5] = 'N'
                    mate[7] = mate[7][0] + '\n'
                    st1 = 0
                    st2 = 0
                else: # qa != na & qa != 0 & qb != 0 => overlap and two adapters
                    (mate[1], mate[5], fixed_Ns) = fix_N_in_overlap(a, b, bb, na, nb, nb - com)
                    fixedns = fixedns + fixed_Ns
                    mate[1] = mate[1][0:qa]
                    mate[3] = mate[3][0:qa] + '\n'
                    #mate[5] = 'N' ## orig
                    #mate[7] = mate[7][0] + '\n' ## orig
                    mate[5] = dnaReverseComplement(mate[1])
                    mate[7] = mate[3][:-1][::-1] + '\n'
                    st1 = qa
                    st2 = nb-qb
                    stn = qa
            else: # qa != na & qb == 0 => overlapping and only one adapter is found (on first read) instead of two adapters
                jj = jj + 1
                (mate[1], mate[5], fixed_Ns) = fix_N_in_overlap(a, b, bb, na, nb, nb - com)
                fixedns = fixedns + fixed_Ns
                mate[1] = mate[1][0:qa] # just trimming. should I set the other read to N?
                mate[3] = mate[3][0:qa] + '\n'
                #mate[5] = 'N'
                #mate[7] = mate[7][0] + '\n'
                st1 = qa
        elif qb != 0: # qa == na & qb != 0 => overlapping and only one adapter is found (on second read) instead of two adapters
                jj = jj + 1
                (mate[1], mate[5], fixed_Ns) = fix_N_in_overlap(a, b, bb, na, nb, nb - com)
                fixedns = fixedns + fixed_Ns
                #mate[1] = 'N'
                #mate[3] = mate[3][0]+'\n'
                mate[5] = mate[5][0:nb-qb] # just trimming. should I set the other read to N?
                mate[7] = mate[7][0:nb-qb] + '\n'
                st1 = nb-qb
                #statn[nb] = statn.get(nb,0) + 1
    else: # if there is no overlap but still I find the adapter => trim
        # brute force
        # searching directly for the adapter
        if qa != -1:
            jj = jj + 1
            m5 = 0
            if qa == 0:
                mate[1] = 'N'
                mate[3] = mate[3][0] + '\n'
                st1 = 0
            else:
                mate[1] = mate[1][0:qa]
                mate[3] = mate[3][0:qa] + '\n'
                st1 = qa

        if qb != -1:
            jj = jj + 1
            m3 = 0
            if qb == nb:
                mate[5] = 'N'
                mate[7] = mate[7][0] + '\n'
                st1 = 0
            else:
                mate[5] = mate[5][:nb-qb]
                mate[7] = mate[7][:nb-qb] + '\n'
                st1 = nb-qb



    mate[1] = mate[1] + '\n'
    mate[5] = mate[5] + '\n'

    x = ''
    if flag_log:
        x1 = ' ' * (nb - com) + ' ' * qa + '*' * (na-qa)+'\n' if com != -1 else ' ' * qa + '*' * (na-qa)+'\n' if m5 != -1 else '\n'
        x2 = '-' * (nb - com) + a + '\n' if com != -1 else a+'\n'
        x3 = b + '-' * (na - com)+'\n' if com != -1 else b + '\n'
        x4 = '*' * qb+'\n'
        x5 = "mismatches in overlapping part = %d \n" % (ml,)
        x6 = "mismatches in 3-end adapter = %d \n" % (m5,)
        x7 = "mismatches in 5-end adapter = %d \n" % (m3,)
        x8 = "length overlapping = %d \n" % (com)
        x = ''.join(["ID READ 1:        "+mate[0],
                       "ID READ 2:        "+mate[4],
                       "SEQ READ 1:       "+a+'\n',
                       "SEQ READ 2:       "+bb+'\n',
                       "FIXED SEQ READ 1: "+mate[1],
                       "FIXED SEQ READ 2: "+mate[5],
                       '\n',
                       'Alignment:\n',
                       x1,
                       x2,
                       x3,
                       x4,
                       x5,
                       x6,
                       x7,
                       x8,
                       '\n',
                       '\n',
                       '\n'
                      ])

    return (mate, st1, st2, stn, jj, x, fixedns)


def norepeats(x, t = 0.30):
    cc = dict()
    tt = float(len(x)-1) * t
    xx = (x[i:i+2] for i in range(0,len(x)-1))
    r = True
    for u in xx:
        cc[u] = cc.get(u,0) + 1
        if cc[u] >= tt:
            r = False
            break
    return r


#
#
#
class param:
    def __init__(self):
        self.reads_overlap = None
        self.wiggle = None
        self.adapter5 = None
        self.adapter3 = None
        self.flag_log = None


#
#
#
def trim_tail_n(s,q,count = 1):
    # trim tails of N from a read

    if s:
        n = len(s)
        if n != 1:
            ts = s.rstrip('N').rstrip('.')
            m = len(ts)
            r = n - m
            ts = ts.lstrip('N').rstrip('.')
            l = m - len(ts)

            if l+r >= count:
                s = s[l:n-r]
                q = q[l:n-r]

    if not s:
        s = "N"
        q = "I"


    return (s+'\n',q+'\n')

#
#
#
def trim_adapter(input_file_1,
                 input_file_2,
                 output_file_1,
                 output_file_2,
                 log_file,
                 align_file,
                 len_adapter,
                 reads_overlap,
                 reads_infer_adapter,
                 threshold_infer_adapter,
                 verbose = False,
                 link = 'soft',
                 shortest_read = 20,
                 trim_n = 3,
                 cpus = 0):
    #
    #
    # finding automatically for adapters
    # - read only the first few millions pairs of reads
    if verbose:
        print >>sys.stderr,"Reading the files for automated finding of adapters..."
        print >>sys.stderr," - ",input_file_1
        print >>sys.stderr," - ",input_file_2
    i = 0
    # looking for adapter
    adapt5 = dict()
    adapt3 = dict()
    #
    for mate in first_reads_from_paired_fastq_file(
        input_file_1,
        input_file_2,
        reads_infer_adapter):
        
        a = mate[1]
        b = dnaReverseComplement(mate[5])
        i = i + 1
        (a5,a3) = fast_alignment_adapter(a, b, len_adapter = len_adapter)
        if a5:
            adapt5[a5] = adapt5.get(a5,0) + 1
        if a3:
            adapt3[a3] = adapt3.get(a3,0) + 1

    #
    # deal with the adapter
    # sort the most common flanking sequences
    #
    a5 = [(v,k) for (k,v) in adapt5.items()]
    a3 = [(v,k) for (k,v) in adapt3.items()]
    # remove the adapters which contain homopolymers!
    a5 = [(v,k) for (v,k) in a5 if len(set(k)) > 2]
    a3 = [(v,k) for (v,k) in a3 if len(set(k)) > 2]
    # remove the adapter which have repeats
    a5 = [(v,k) for (v,k) in a5 if norepeats(k,0.3)]
    a3 = [(v,k) for (v,k) in a3 if norepeats(k,0.3)]
    # sort based on count
    a5 = sorted(a5, reverse = True)
    a3 = sorted(a3, reverse = True)

    adapter5 = ""
    adapter3 = ""
    adapter3reverse = ""
    if a5:
        count_adapter5 = a5[0][0]
        adapter5 = a5[0][1]
        if a5[0][0] > threshold_infer_adapter * float(reads_infer_adapter):
            if verbose:
                print >>sys.stderr,"Adapter 3' end found!\n   [%s...]\n   [reverse-complement:...%s]\n   [count=%d/%d] (%.5f%%)" % (adapter5,dnaReverseComplement(adapter5),a5[0][0],reads_infer_adapter,100*float(a5[0][0])/reads_infer_adapter)
        else:
            if verbose:
                print >>sys.stderr,"Found [%s...]\n   [reverse-complement=...%s]\n   [count=%d/%d (%.5f%%)] but does not look like 3-end adapter! (too low count)" % (adapter5,dnaReverseComplement(adapter5), a5[0][0], reads_infer_adapter,100*float(a5[0][0])/reads_infer_adapter)
            adapter5 = None
    if a3:
        count_adapter3 = a3[0][0]
        adapter3 = a3[0][1]
        adapter3reverse = dnaReverseComplement(adapter3)
        if a3[0][0] > threshold_infer_adapter * float(reads_infer_adapter):
            if verbose:
                print >>sys.stderr,"Adapter 5' end found!\n   [%s...]\n   [reverse-complement:...%s]\n   [count=%d/%d (%.5f%%)]" % (adapter3reverse, adapter3,a3[0][0],reads_infer_adapter,100*float(a3[0][0])/reads_infer_adapter)
        else:
            if verbose:
                print >>sys.stderr,"Found [%s...]\n   [reverse-complement=...%s]\n   [count=%d/%d (%.5f%%)] but does not look like 5-end adapter! (too low count)" % (adapter3reverse, adapter3, a3[0][0], reads_infer_adapter,100*float(a3[0][0])/reads_infer_adapter)
            adapter3 = None

    #
    # remove the adapters
    #
    threshold = 10 # if more than 10% of reads have adapter then use more than one process, else use only one process
    if (reads_infer_adapter != 0) and adapter5 and adapter3 and ((100*float(a5[0][0])/reads_infer_adapter >= threshold) or (100*float(a3[0][0])/reads_infer_adapter >= threshold)):
        if cpus == 0:
            cpus = multiprocessing.cpu_count()
    else:
        if verbose:
            print >>sys.stderr,"NOTE: Too few adapters found in order to use several processes! Only one CPU will be used!"
        cpus = 1
    if verbose:
        print >>sys.stderr,"Using",cpus,"process(es)..."
    #
    flag_log = False
    log = None
    if align_file:
        flag_log = True
        log = lines_to_file(align_file)

    log_stat = None
    if log_file:
        log_stat = lines_to_file(log_file)

    stat = dict()
    statn = dict()
    if adapter5 and adapter3:
        # I found both adapters
        #
        if verbose:
            print >>sys.stderr,"Scanning for adapters..."
        i = 0
        j = 0
        last_j = j

        out_1 = lines_to_file(output_file_1)
        out_2 = lines_to_file(output_file_2)


        para = param()
        para.reads_overlap = reads_overlap
        para.wiggle = 2
        para.adapter5 = adapter5
        para.adapter3 = adapter3
        para.flag_log = flag_log

        #

        all_fixed = 0
        give = None
        if cpus == 1:

            for stuff in itertools.imap(
                                    compute,
                                    itertools.izip_longest(
                                        reads_from_paired_fastq_file(
                                            input_file_1, 
                                            input_file_2
                                        ),
                                        [],
                                        fillvalue = para
                                    )
                                ):

                i = i + 1
                mate = stuff[0]
                st1 = stuff[1]
                st2 = stuff[2]
                stn = stuff[3]
                jj = stuff[4]
                xx = stuff[5]
                fixed = stuff[6]
                all_fixed = all_fixed + fixed

                if st1 != -1:
                    stat[st1] = stat.get(st1,0) + 1
                if st2 != -1:
                    stat[st2] = stat.get(st2,0) + 1
                if stn != -1:
                    statn[stn] = statn.get(stn,0) + 1
                j = j + jj

                if trim_n:
                    # do N trimming from both ends
                    mm1 = mate[1].rstrip('\r\n')
                    if mm1.startswith('N') or mm1.endswith('N'):
                        mm2 = mate[3].rstrip('\r\n')
                        (mate[1], mate[3]) = trim_tail_n(mm1,mm2,trim_n)

                    mm1 = mate[5].rstrip('\r\n')
                    if mm1.startswith('N') or mm1.endswith('N'):
                        mm2 = mate[7].rstrip('\r\n')
                        (mate[5], mate[7]) = trim_tail_n(mm1,mm2,trim_n)


                if len(mate[1]) < shortest_read + 1:
                    out_1.add_line(mate[0])
                    out_1.add_line(empty_read[1])
                    out_1.add_line(empty_read[2])
                    out_1.add_line(empty_read[3])
                else:
                    out_1.add_lines(mate[0:4])

                if len(mate[5]) < shortest_read + 1:
                    out_2.add_line(mate[4])
                    out_2.add_line(empty_read[1])
                    out_2.add_line(empty_read[2])
                    out_2.add_line(empty_read[3])
                else:
                    out_2.add_lines(mate[4:8])


                if flag_log and last_j != j:
                    log.add_line(xx)

                last_j = j
                if verbose and i % 1 == 10000000: # 10000000
                    print >>sys.stderr,"   %d pair-reads [ %d (%f%%) reads trimmed ]" % (i,j,100*float(j)/float(2*i))

            
        else:
            pool = multiprocessing.Pool(processes = cpus)

            for stuff in pool.imap_unordered(
                                            compute,
                                            itertools.izip_longest(
                                                reads_from_paired_fastq_file(
                                                    input_file_1,
                                                    input_file_2
                                                ),
                                                [],
                                                fillvalue = para
                                            ),
                                        chunksize = 100
                                        ):

                i = i + 1
                mate = stuff[0]
                st1 = stuff[1]
                st2 = stuff[2]
                stn = stuff[3]
                jj = stuff[4]
                xx = stuff[5]
                fixed = stuff[6]
                all_fixed = all_fixed + fixed

                if st1 != -1:
                    stat[st1] = stat.get(st1,0) + 1
                if st2 != -1:
                    stat[st2] = stat.get(st2,0) + 1
                if stn != -1:
                    statn[stn] = statn.get(stn,0) + 1
                j = j + jj

                if trim_n:
                    # do N trimming from both ends
                    mm1 = mate[1].rstrip('\r\n')
                    if mm1.startswith('N') or mm1.endswith('N'):
                        mm2 = mate[3].rstrip('\r\n')
                        (mate[1], mate[3]) = trim_tail_n(mm1,mm2,trim_n)

                    mm1 = mate[5].rstrip('\r\n')
                    if mm1.startswith('N') or mm1.endswith('N'):
                        mm2 = mate[7].rstrip('\r\n')
                        (mate[5], mate[7]) = trim_tail_n(mm1,mm2,trim_n)


                if len(mate[1]) < shortest_read + 1:
                    out_1.add_line(mate[0])
                    out_1.add_line(empty_read[1])
                    out_1.add_line(empty_read[2])
                    out_1.add_line(empty_read[3])
                else:
                    out_1.add_lines(mate[0:4])

                if len(mate[5]) < shortest_read + 1:
                    out_2.add_line(mate[4])
                    out_2.add_line(empty_read[1])
                    out_2.add_line(empty_read[2])
                    out_2.add_line(empty_read[3])
                else:
                    out_2.add_lines(mate[4:8])


                if flag_log and last_j != j:
                    log.add_line(xx)

                last_j = j
                if verbose and i % 1 == 10000000: # 10000000
                    print >>sys.stderr,"   %d pair-reads [ %d (%f%%) reads trimmed ]" % (i,j,100*float(j)/float(2*i))


            
            pool.close()
            pool.join()




        out_1.close()
        out_2.close()


        # total number of joined reads
        s = 0
        for (k,v) in statn.items():
            s = s + v

        if verbose:
            print >>sys.stderr,"Total count reads = %i" % (2*i,)
            print >>sys.stderr,"Count trimmed reads = %i [%f%%]" % (j,100*float(j)/float(2*i))
            print >>sys.stderr,"Count joined pair-reads = %i [%f%%]" % (s,100*float(s)/float(i))
            print >>sys.stderr,"Count of fixed Ns = %i" % (all_fixed,)

        if log_stat:
            log_stat.add_line("Input file read 1: %s\n" % (input_file_1,))
            log_stat.add_line("Input file read 2: %s\n" % (input_file_2,))
            log_stat.add_line("Adapter 3' end = [%s...] [reverse-complement=...%s] [count=%d/%d]\n" % (adapter5,dnaReverseComplement(adapter5),a5[0][0],reads_infer_adapter))
            log_stat.add_line("Adapter 5' end = [%s...] [reverse-complement=...%s] [count=%d/%d]\n" % (adapter3reverse, adapter3,a3[0][0],reads_infer_adapter))
            log_stat.add_line("Total count reads = %i\n" % (2*i,))
            log_stat.add_line("Count trimmed reads = %i [%f%%]\n" % (j,100*float(j)/float(2*i)) )
            log_stat.add_line("Count not-trimmed reads = %i [%f%%]\n" % (2*i-j,100*float(2*i-j)/float(2*i)) )
            log_stat.add_line("Count joined pair-reads = %i [%f%%]\n" % (s,100*float(s)/float(i)))
            log_stat.add_line("Count of fixed Ns = %i\n" % (all_fixed,))
            log_stat.add_line("--------------------------------------------------------------------------------\n")
            log_stat.add_line("----------------TRIMMED READS---------------------------------------------------\n")
            log_stat.add_line("size trimmed-read\tcount\tpercentage[%]\tcumulated percentage[%]\n")
            w = 0
            for (k,v) in sorted(stat.items()):
                x = 100*float(v)/float(2*i)
                w = w + x
                log_stat.add_line("%d\t%d\t%f\t%f\n" % (k,v,x,w))
            log_stat.add_line("--------------------------------------------------------------------------------\n")
            log_stat.add_line("----------------JOINED PAIR-READS (overlapping and adapters found)------\n")
            log_stat.add_line("size joined pair-read\tcount\tpercentage[%]\tcumulated percentage[%]\n")
            w = 0
            for (k,v) in sorted(statn.items()):
                x = 100*float(v)/float(i)
                w = w + x
                log_stat.add_line("%d\t%d\t%f\t%f\n" % (k,v,x,w))
            log_stat.add_line("================================================================================\n")
            log_stat.close()

        if flag_log:
            log.close()

    else:
        # create symbolic links and move on
        if verbose:
            print >>sys.stderr,"No trimming was done !"
            print >>sys.stderr,"NOTE: In case that it is known a priori that there are adaptors then "
            print >>sys.stderr,"      try to run it again with a shorter adapter, see '--len_adapter' "
            print >>sys.stderr,"      of 'remove_adapter.py'!"

        linkit(input_file_1,output_file_1,kind=link)
        linkit(input_file_2,output_file_2,kind=link)

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

    description = """It trims (finds and/or removes by trimming) the paired-reads (FASTQ file
produced by Illumina Solexa) which overlap and contain the adapter. The adapter
is found automatically. Also the partial overlapping between a short read and
the adapter sequence is handled."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2017 Daniel Nicorici

"""

    version = "%prog 0.98 beta"

    parser = MyOptionParser(usage       = usage,
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
                      help = """The input FASTQ file containing the reads from 3' fragment end (i.e. 5'-3' orientation for read 1 and 3'-5' for read 2 which needs to be reversed-complemented).""")

    parser.add_option("-f","--output_1",
                      action = "store",
                      type = "string",
                      dest = "output_1_filename",
                      help = """The output FASTQ file where the reads are trimmed.""")

    parser.add_option("-r","--output_2",
                      action = "store",
                      type = "string",
                      dest = "output_2_filename",
                      help = """The output FASTQ file where the reads are trimmed.""")

    parser.add_option("-l","--log",
                      action = "store",
                      type = "string",
                      dest = "log_filename",
                      help = """It outputs a detalied statistics of the trimming.""")

    parser.add_option("-g","--alignment_log",
                      action = "store",
                      type = "string",
                      dest = "alignment_log_filename",
                      help = """It outputs also the alignment for each found overlapping.""")

#    parser.add_option("--mismatches",
#                      action = "store",
#                      type = "int",
#                      dest = "mismatches",
#                      default = 3,
#                      help = """Maximum number or mismatches to be allowed in the overlapping region. Default is %default.""")

    parser.add_option("-o","--reads_overlap",
                      action = "store",
                      type = "int",
                      dest = "reads_overlap",
                      default = 10,
                      help = """The minimum length of the region which is considered an overlap. Default is %default.""")


    parser.add_option("-n","--len_adapter",
                      action = "store",
                      type = "int",
                      dest = "len_adapter",
                      default = 13,
                      help = """The length of the adapter which is found automaticaly and further used for trimming. Default is %default.""")

    parser.add_option("-s","--shortest_read",
                      action = "store",
                      type = "int",
                      dest = "shortest_read",
                      default = 20,
                      help = """The reads stricly shorter than %default after the adapter trimming will be removed (i.e. replaced with a read of length one containing the sequence 'N'). Default is %default.""")

    parser.add_option("-t","--threshold_infer_adapter",
                      action = "store",
                      type = "float",
                      dest = "threshold_infer_adapter",
                      default = 0.0001,
                      help = """The percentage of reads which should contain the found candidate-adapter (during the automatic adapter-finding step) in order to be considered a real adapter and used further for trimming. The range is [0..1]. Default is %default.""")


    parser.add_option("-i","--reads_infer_adapter",
                      action = "store",
                      type = "int",
                      dest = "reads_infer_adapter",
                      default = 3000000,
                      help = """The number of first reads which are used for finding automatically the adapter. If it is set to 0 then all the reads from the files are used. Default is %default.""")


    parser.add_option("-x","--trim-n",
                      action = "store",
                      type = "int",
                      dest = "trim_n",
                      default = 0,
                      help = """The number of Ns needed to be found in order the trimming of Ns (from both reads ends) is done. Default is %default.""")

    parser.add_option("-k","--link",
                      action = "store",
                      type = "string",
                      dest = "link",
                      default = 'soft',
                      help = """The type of link between the input and output file when there are no changes done. The choices are ['soft','hard']. Default is '%default'.""")





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

    parser.add_option("-p", "--processes",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 0,
                      help = """Maximum number of parallel processes/CPUs to be used for computations. In case of value 0 then the program will try to use, if it see fit, all the CPUs which are found. The default value is %default.""")

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
    trim_adapter(options.input_1_filename,
                 options.input_2_filename,
                 options.output_1_filename,
                 options.output_2_filename,
                 options.log_filename,
                 options.alignment_log_filename,
                 options.len_adapter,
                 options.reads_overlap,
                 options.reads_infer_adapter,
                 options.threshold_infer_adapter,
                 options.verbose,
                 options.link,
                 options.shortest_read,
                 options.trim_n,
                 options.processes)


if __name__ == '__main__':
    main()
    
    
    
#
