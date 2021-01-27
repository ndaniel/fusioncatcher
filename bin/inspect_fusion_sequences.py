#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a list of chromosomal coordinates of fusion genes and it labels them if they have a very low complexity fusion sequence.

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

import os
import sys
import optparse
import gzip
import math


#
#
#
polyA = "A"
polyC = "C"
polyG = "G"
polyT = "T"

#
#
#
def counter(sequence,nucleotide = 2):
    freq = {}
    len_sequence = len(sequence) - nucleotide + 1
    for i in xrange(0,len_sequence):
        e = sequence[i:i+nucleotide]
        freq[e] = freq.get(e,0) + 1
    return freq

#
#
#
def plus(a,b):
    ks = set(a.keys()+b.keys())
    d = {}
    for k in ks:
        d[k] = a.get(k,0) + b.get(k,0)
    return d

#
#
#
def minus(a,b):
    ks = set(a.keys()+b.keys())
    d = {}
    for k in ks:
        x = a.get(k,0) - b.get(k,0)
        if x > 0:
            d[k] = x
    return d

#
#
#
def plusminus(a,b,c):
    ks = set(a.keys()+b.keys()+c.keys())
    d = {}
    for k in ks:
        x = a.get(k,0) + b.get(k,0) - c.get(k,0)
        if x > 0:
            d[k] = x
    return d

#
#
#
def bits(d):
    v = d.values()
    n = sum(v)
    if n != 0:
        v = sum([-(float(e)/float(n))*math.log(float(e)/float(n),2) for e in v if e != 0])
    else:
        v = 0
    return v

#
#
#
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



#
#
#
def evaluate_fusion_sequence(
             input_file,
             output_file,
             window_length = 24,
             window_overlap = 12,
             kmer = 2,
             threshold = 2.0,
             threshold2 = 2.0,
             poly = 15,
             remove_poly_filename = '',
             verbose = True):

    polyA = poly * "A"
    polyC = poly * "C"
    polyG = poly * "G"
    polyT = poly * "T"

    if verbose:
        print >>sys.stderr,"Reading the list of fusion genes..."
    # read the list of fusion genes and their chromosomal positions
    data = [line.rstrip('\r\n').split('\t') for line in file(input_file,'r').readlines() if line.rstrip('')]
    header = data.pop(0)

    if remove_poly_filename:
        file(remove_poly_filename,'a').write("\n\n\nCandidate fusions removed due to PolyA/C/G/T and long repeats:\n===============================================\n")


    if data:
        res = []
        for line in data:
            s = line[14].replace(" ","").replace("*","").replace("N","")
            label = line[2]
            l1 = codelength(s,window_length,window_overlap,kmer)
            if l1 < threshold:
                label = label + "," if label else label
                label = label + 'short_repeats,sr%.2f' % (l1,)
                #print "l1",s,l1
            if len(s) > 50:
                l2 = codelength(s,w=50,o=40,kmer=9)
                if l2 < threshold2:
                    label = label + "," if label else label
                    label = label + 'long_repeats,lr%.2f' % (l2,)
                #print "l2",s,l2
            poly_found = False
            if s.find(polyA) != -1:
                label = label + "," if label else label
                label = label + 'polyA'
                poly_found = True
            if s.find(polyC) != -1:
                label = label + "," if label else label
                label = label + 'polyC'
                poly_found = True
            if s.find(polyG) != -1:
                label = label + "," if label else label
                label = label + 'polyG'
                poly_found = True
            if s.find(polyT) != -1:
                label = label + "," if label else label
                label = label + 'polyT'
                poly_found = True
            
            if remove_poly_filename:
                if poly_found or label.find("long_repeats") != -1:
                    uf = line[0:2]+[label]+line[3:]
                    file(remove_poly_filename,'a').write('\t'.join(uf)+'\n')
                else:
                    res.append(line[0:2]+[label]+line[3:])
            else:
                res.append(line[0:2]+[label]+line[3:])

        if verbose:
            print >>sys.stderr,"Parsing the GTF file..."

        res.insert(0,header)
        file(output_file,'w').writelines(['\t'.join(line)+'\n' for line in res])
    else:
        file(output_file,'w').write('\t'.join(header)+'\n')



################################################################################
################################################################################
################################################################################
def main():

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options]"

    description = """It takes as input a list of chromosomal coordinates of fusion genes and it labels them if they have a very low complexity fusion sequence."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2021 Daniel Nicorici

"""



    version = "%prog 0.04 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
                            )




    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file containing the fusion (chromosomal) coordinates for each fusion genes.""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file where the frame predictions are written. """)

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

    parser.add_option("-p","--poly",
                      action = "store",
                      type = "int",
                      dest = "poly",
                      default = 15,
                      help = """The minimum length of the polyN. Default is %default.""")


    parser.add_option("-t","--threshold",
                      action = "store",
                      type = "float",
                      dest = "codelength_threshold",
                      default = 1.4, # original was 2.2
                      help = """Any window which compresses less this threshold is considered to contain a short tandem repeat and the read will be filtered out. Default is %default.""")

    parser.add_option("-s","--threshold2",
                      action = "store",
                      type = "float",
                      dest = "codelength_threshold2",
                      default = 4.5, # original was 2.2
                      help = """Any window which compresses less this threshold is considered to contain a short tandem repeat and the read will be filtered out. Default is %default.""")


    parser.add_option("-x","--remove-poly",
                      action = "store",
                      type = "string",
                      dest = "remove_poly_filename",
                      help = """The fusions that are found to contains polyA/C/G/T will be filtered out and written in this file. Default is %default.""")


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


    ( options , args ) = parser.parse_args()

    # validate options
    if not (
            options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")


    #
    evaluate_fusion_sequence(
                    options.input_filename,
                    options.output_filename,
                    options.window_length,
                    options.window_overlap,
                    options.kmer,
                    options.codelength_threshold,
                    options.codelength_threshold2,
                    options.poly,
                    options.remove_poly_filename,
                    options.verbose)


if __name__ == '__main__':
    main()
    
    
#    
