#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It converts a FASTQ file to FASTA format.



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
import optparse
import Bio.SeqIO

def fq2fa(f_in,f_out):
    if hasattr(Bio.SeqIO,'convert'):
        counts=Bio.SeqIO.convert(f_in, "fastq-solexa", f_out, "fasta")
    else:
        print "Bio.SeqIO.convert() not supported!"
        print "Trying to go around it!"
        #record=list(Bio.SeqIO.parse(file(f_in,'r'), "fastq-solexa"))
        #Bio.SeqIO.write(record,file(f_out, "w"),"fasta")
        input_handle=open(f_in, "rU")
        output_handle=open(f_out, "w")
        sequences=Bio.SeqIO.parse(input_handle,"fastq-solexa")
        counts=Bio.SeqIO.write(sequences,output_handle,"fasta")
        output_handle.close()
        input_handle.close()
    print "Converted %i records" % counts

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It converts a FASTQ file to a FASTA file."""
    version="%prog 0.10 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ format (Solexa).""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTA file.""")

    (options,args)=parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    # running
    print "Converting..."
    fq2fa(options.input_filename,options.output_filename)

    print "The end."
