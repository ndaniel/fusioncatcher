#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It splits a FASTA file into several equally sized FASTA files.



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


import sys
import os
import optparse
import math


buffer_size = 10**8

def sequence(afile):
    fid = open(afile,"r")
    box = []
    m = 0
    while True:
        lines = fid.readlines(buffer_size)
        if not lines:
            if box:
                yield (box,m)
            break
        for line in lines:
            if line.startswith('>'):
                if box and m:
                    yield (box,m)
                if line.rstrip("\r\n"):
                    box = [line]
                else:
                    box = []
                m = 0
            else:
                if line.rstrip('\r\n'):
                    box.append(line)
                    m = m + len(line) - 1
    fid.close()


def int2str(x,n=3):
    if n < 1:
        n = 1
    x = str(x)
    return '0' * int(n - len(x)) + x


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It splits a FASTA file into several equally sized FASTA files such that all are below a given size threshold."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file in FASTA format.""")

    parser.add_option("--threshold","-t",
                      action = "store",
                      type = "int",
                      dest = "threshold",
                      default = 2**31-1000,
                      help = """Any splitted file should not have the size (in bytes) over this threshold if that is possible. If this is set to 1 then one sequence will be written in each output FASTA file. Default is %default.""")

    parser.add_option("--seq_per_fasta","-f",
                      action = "store",
                      type = "int",
                      dest = "seq_per_fasta",
                      help = """This is another way to split the input FASTA file by giving the number of sequences per output FASTA file. If this option is used it will over ride the '--threshold'.""")

    parser.add_option("--size","-s",
                      action = "store",
                      type = "string",
                      dest = "size_input_fasta_filename",
                      help = """A file containing the number of nucleotides which are in the input FASTA file. If it is not given then the file size will be used.""")

    parser.add_option("--seqs","-q",
                      action = "store",
                      type = "string",
                      dest = "seqs_input_fasta_filename",
                      help = """A file containing the number of sequences which are in the input FASTA file. If it is not given then the file size divided by 50 will be used. It is used/needed by '--seq_per_fasta'""")

    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file which contains the paths to the splitted FASTA files.""")


    parser.add_option("--output-max-lens","-x",
                      action = "store",
                      type = "string",
                      dest = "output_max_lens_filename",
                      help = """The output file which contains the paths to the maximum lengths of the splitted FASTA files.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    # start

    n = 0
    bucket = 0

    # reading
    (xdir,xfile) = os.path.splitext(options.output_filename)


    # find size of input FASTA file
    size = os.path.getsize(options.input_filename)

    #
    if options.threshold < 10:
        options.seq_per_fasta = 1

    if options.seq_per_fasta:
        n = int(math.ceil(float(size) / float(30)))
        if options.seqs_input_fasta_filename:
            n = int(file(options.seqs_input_fasta_filename,'r').readline().strip())
    else:
        # threshold size
        threshold = options.threshold

        if options.size_input_fasta_filename:
            size = int(file(options.size_input_fasta_filename,"r").readline().strip())

        # compute/estimate the approximate size of the chunk (chunks should be approximately equally)
        n = int(math.ceil(float(size) / float(threshold)))
        bucket = int(math.ceil(float(size) / float(n)))


    truck = []
    level = 0
    xid = 0

    no_digits = int(math.ceil(math.log(n+1,10)))

    f = options.output_filename+'.'+int2str(xid,no_digits)
    xlist = [f]
    fod = open(f,"w")
    
    max_lens = False
    if options.output_max_lens_filename:
        max_lens = True
    maxlist = []
    
    temp = 0
    seq_count = 0
    written = False
    max_x = 0
    if options.seq_per_fasta:
        for (seq,x) in sequence(options.input_filename):
            if x > max_x:
                max_x = x
            truck.extend(seq)
            seq_count = seq_count + 1
            temp = temp + x
            if seq_count == options.seq_per_fasta:
                fod.writelines(truck)
                truck = []
                seq_count = 0
                temp = 0
                #
                fod.close()
                if max_lens:
                    fmax = options.output_max_lens_filename+'.'+int2str(xid,no_digits)
                    file(fmax,'w').write("%d" %(max_x,))
                    maxlist.append(fmax)
                    max_x = 0
                    
                xid = xid + 1
                f = options.output_filename+'.'+int2str(xid,no_digits)
                fod = open(f,'w')
                written = False
                xlist.append(f)
            elif temp >= buffer_size:
                fod.writelines(truck)
                truck = []
                temp = 0
                written = True
    else:
        for (seq,x) in sequence(options.input_filename):
            if x > max_x:
                max_x = x
            truck.extend(seq)
            level = level + x
            temp = temp + x
            if level >= bucket:
                fod.writelines(truck)
                truck = []
                level = 0
                temp = 0
                #
                fod.close()
                if max_lens:
                    fmax = options.output_max_lens_filename+'.'+int2str(xid,no_digits)
                    file(fmax,'w').write("%d" %(max_x,))
                    maxlist.append(fmax)
                    max_x = 0
                
                xid = xid + 1
                f = options.output_filename+'.'+int2str(xid,no_digits)
                fod = open(f,'w')
                written = False
                xlist.append(f)
            elif temp >= buffer_size:
                fod.writelines(truck)
                truck = []
                temp = 0
                written = True
    if truck:
        fod.writelines(truck)
        written = True

    if not fod.closed:
        fod.close()
        if max_lens:
            fmax = options.output_max_lens_filename+'.'+int2str(xid,no_digits)
            file(fmax,'w').write("%d" %(max_x,))
            maxlist.append(fmax)
            max_x = 0


    if written == False: # the last opened file might be empty
        last = xlist.pop()
        os.remove(last)
        if max_lens:
            last = maxlist.pop()
            os.remove(last)

    # EXTRA CHECKING
    if xlist:
        if os.path.isfile(xlist[-1]) and os.path.getsize(xlist[-1]) == 0:
            last = xlist.pop()
            os.remove(last)
            if max_lens:
                last = maxlist.pop()
                os.remove(last)

    file(options.output_filename,"w").writelines([el+'\n' for el in xlist])
    
    if max_lens:
        file(options.output_max_lens_filename,"w").writelines([el+'\n' for el in maxlist])
    #
    # end
    #
