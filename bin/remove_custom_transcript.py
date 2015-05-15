#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Remove human transcripts which are wrongly annotated.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

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
import symbols

def remove(outdir,
           transcript_id = ''):

    headex = dict([(line.rstrip('\r\n'),i) for i,line in enumerate(file(os.path.join(outdir,'exons_header.txt'),'r').readlines()) if line.rstrip('\r\n')])
    exons = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(outdir,'exons.txt'),'r').readlines() if line.rstrip('\r\n')]
    tid = headex['ensembl_transcript_id']
    gid = headex['ensembl_gene_id']

    headge = dict([(line.rstrip('\r\n'),i) for i,line in enumerate(file(os.path.join(outdir,'genes_header.txt'),'r').readlines()) if line.rstrip('\r\n')])
    genes = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(outdir,'genes.txt'),'r').readlines() if line.rstrip('\r\n')]


    data = []
    fix = []
    gene_id = [line[gid] for line in exons if line[tid] == transcript_id]
    falg = False
    if len(set(gene_id)) == 1:
        gene_id = gene_id.pop()
        for line in exons:
            if line[gid] == gene_id:
                if line[tid] != transcript_id:
                    fix.append(line)
            else:
                data.append(line)
        if fix:
            # find new start and end position for gene after the transcript has been removed
            p1id = headex['start_position']
            p2id = headex['end_position']
            t1id = headex['transcript_start']
            t2id = headex['transcript_end']
            t1 = fix[0][t1id]
            t2 = fix[0][t2id]
            for line in fix:
                if line[t1id] < t1:
                    t1 = line[t1id]
                if line[t2id] > t2:
                    t2 = line[t2id]
            for line in fix:
                line[p1id] = t1
                line[p2id] = t2
            # join everything
            data.extend(fix)
            # fix the genes.txt
            p1id = headge['start_position']
            p2id = headge['end_position']
            gid = headge['ensembl_gene_id']
            flag = False
            for line in genes:
                if line[gid] == gene_id:
                    line[p1id] = t1
                    line[p2id] = t2
                    flag = True
            if flag:
                file(os.path.join(outdir,'genes.txt'),'w').writelines(['\t'.join(line)+'\n' for line in genes])
        file(os.path.join(outdir,'exons.txt'),'w').writelines(['\t'.join(line)+'\n' for line in data])

        # remove it from the GTF file also
        gtf = [line for line in file(os.path.join(outdir,'organism.gtf'),'r').readlines() if line.find(transcript_id) == -1]
        file(os.path.join(outdir,'organism.gtf'),'w').writelines(gtf)
    else:
        print >>sys.stderr, "WARNING: %s transcript not found!" % (transcript_id,)




if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Remove human transcripts which are wrongly annotated."""
    version = "%prog 0.13 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the list of allowed candidate fusion genes is generated, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the list of allowed candidate fusion genes is generated. Default is '%default'.""")

    parser.add_option("--transcripts",
                      action="store",
                      type="string",
                      dest="transcripts",
                      help="""List of transcripts which should be removed.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    print "Remove transcripts which are wrongly annotated..."

    if options.transcripts:
        gl = [line.rstrip("\r\n") for line in file(options.transcripts,'r').readlines() if line.rstrip("\r\n")]
        for g in gl:
            remove(options.output_directory,g)

    if options.organism.lower() == 'homo_sapiens':
        remove(options.output_directory,"ENST00000467125") # ENST00000467125 is listed as GOPC, but is actually the GOPC/ROS1
        remove(options.output_directory,"ENST00000507166") # ENST00000507166 is listed as FIP1L1 but is really the FIP1L1/PDGFRA
        remove(options.output_directory,"ENST00000621209") # ENST00000621209 is listed as CEL but is really the CEL/CELP

    #
