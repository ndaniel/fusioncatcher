#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It compresses using lossy compression the ids of all reads from a input FASTQ file (using the read index/count).
The reads ids have all the ids in alphabetically order.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

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
import gzip
import string
import math
import itertools

def generate_id(t, lowercase = False, interleaved = True, no12 = False):
    digits = string.digits + string.ascii_uppercase
    if lowercase:
        digits = digits + string.ascii_lowercase
    l = len(digits)
    if interleaved:
        r = int(math.ceil(math.log(float(t+10)/float(2),l)))
        if no12:
            for el in itertools.product(digits,repeat = r):
                x = ''.join(el)
                yield "@%s\n" % (x,)
                yield "@%s\n" % (x,)
        else:
            for el in itertools.product(digits,repeat = r):
                x = ''.join(el)
                yield "@%s/1\n" % (x,)
                yield "@%s/2\n" % (x,)
    else:
        r = int(math.ceil(math.log(float(t+5),l)))
        for el in itertools.product(digits,repeat = r):
            x = ''.join(el)
            yield "@%s\n" % (x,)




if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It compresses using lossy compression the ids of all reads from a input FASTQ file (using the read index/count). The compressed reads ids have all the ids in alphabetically order."""
    version = "%prog 0.15 beta"

    parser = optparse.OptionParser(
                usage = usage,
                description = description,
                version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ Solexa file (also given thru STDOUT or as gzipped file). By default, it is assumed that the input reads are shuffled/interleaved (that is read id 1 is followed by read id 2 where read 1 and read 2 form a pair). """)

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text file containg the changed reads ids (also given thru STDOUT or as gzipped file).""")

    parser.add_option("--count-reads","-n",
                      action = "store",
                      type = "string",
                      dest = "count",
                      help="""The total number of reads in the input file. This is used in order to compress the best the reads ids.""")

    parser.add_option("--no12","-s",
                      action = "store_true",
                      dest = "no12",
                      default = False,
                      help="""If this is set than no /1 and /2 will be added to the compressed reads ids. It has an effect only on interleaved inputs.""")

    parser.add_option("--not-interleaved","-x",
                      action = "store_true",
                      dest = "not_interleaved",
                      default = False,
                      help="""If it is set then the input reads from the input FASTQ files are not interleaved. Also no /1 or /2 is added to the reads ids.""")

    parser.add_option("--lowercase","-l",
                      action = "store_true",
                      dest = "lowercase",
                      default = False,
                      help="""If this is set then also lowercase charcaters will be used for read ids in FASTQ files.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("Input and output files should be specified!")
        sys.exit(1)


    fin = None
    if options.input_filename == '-':
        fin = sys.stdin
    elif options.input_filename.lower().endswith('.gz'):
        fin = gzip.open(options.input_filename,'r')
    else:
        fin = open(options.input_filename,'r')

    fou = None
    if options.output_filename == '-':
        fou = sys.stdout
    elif options.output_filename.lower().endswith('.gz'):
        fou = gzip.open(options.output_filename,'w')
    else:
        fou = open(options.output_filename,'w')

    n = 0
    if options.count:
        if os.path.isfile(options.count):
            n = sum([int(line.strip()) for line in file(options.count,'r').readlines() if line.strip()])
        else:
            n = int(options.count)
    if not n:
        print >>sys.stderr,"Warning: Cannot read/use the '--count-reads' option! Therefore the counting will be done by the script!"
        # then use the file size as a proxy for the number of reads
        if options.input_filename != '-':
            n = 0
            fid = open(options.input_filename,'r')
            sb = 10**8
            while True:
                gc.disable()
                lines = fid.readlines(sb)
                gc.enable()
                if not lines:
                    break
                n = n + len(lines)
            fid.close()
            n = n / 4
        else:
            print >>sys.stderr,"ERROR: '--count-reads' option is needed to be specified!"
            sys.exit(1)
    i = 0
    ids = generate_id(n, lowercase = options.lowercase, interleaved = (not options.not_interleaved), no12=options.no12)
    sb = 10**8
    while True:
        gc.disable()
        lines = fin.readlines(sb)
        gc.enable()
        if not lines:
            break
        gc.disable()
        lines = [ids.next() if (j+i)%4 == 0 else '+\n' if (j+i)%4 == 2 else line for (j,line) in enumerate(lines)]
        gc.enable()
        i = i + len(lines)
        fou.writelines(lines)
    fin.close()
    fou.close()

    #
