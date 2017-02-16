#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a fasta file containing the transcriptome and another file with
genes ids and it gives the transcript sequences of the given genes as a fasta file.



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


import os
import sys
import Bio
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet
import optparse



if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It takes as input a fasta file containing the transcriptome and another file with genes ids and it gives the transcript sequences of the given genes as a fasta file."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

    parser.add_option(
        "--input_genes","-g",
        action = "store",
        type = "string",
        dest = "input_genes_ids",
        help = """Input file with Ensembl genes ids.""")

    parser.add_option(
        "--input_transcriptome","-i",
        action = "store",
        type = "string",
        dest = "input_transcriptome_fasta",
        help = """Input FASTA genome containing the genome sequences.""")


    parser.add_option(
        "--output",
        action = "store",
        type = "string",
        dest = "output_filename",
        help = """The output FASTA file where the genes sequences are written.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_genes_ids and
            options.input_transcriptome_fasta and
            options.output_filename
            ):
        parser.print_help()
        parser.error("Missing argument(s)!")
        sys.exit(1)

    #
    # Process the input list of genes ids
    # first column is "gene ids" and second column (which is optional) is the "gene symbol"
    ids = [line.rstrip('\r\n').split('\t') for line in file(options.input_genes_ids,'r').readlines() if line.rstrip('\r\n')]
    ids = dict([(line[0],'') if len(line) <2 else (line[0],line[1]) for line in ids])

    transcripts = []
    output_handle = open(os.path.join(options.output_filename), "w")
    for record in Bio.SeqIO.parse(open(options.input_transcriptome_fasta, "rU"), "fasta") :
        t = record.id.partition(';')[2]
        if t in ids:
            g = ids[t]
            if g:
                g = "%s;%s" % (g,record.id)
            else:
                g = record.id
            transcripts.append(Bio.SeqRecord.SeqRecord(record.seq,id=g,name="",description=""))
    Bio.SeqIO.write(transcripts, output_handle, "fasta")
    output_handle.close()
    #
