#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the fusion genes found in CCLE from depmap.org portal (from: https://depmap.org/portal/download/).



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
import urllib2
import gzip
import socket
import optparse
import concatenate
import symbols
import shutil
import datetime


def fix(g):
    u = g
    if g.upper() == "IGL-@":
        u = "IGL@"
    elif g.upper() == "IGH-@":
        u = "IGH@"
    return u

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the human fusion genes found in 675 cell lines from article: C. Klijn et al., A comprehensive transcriptional portrait of human cancer cell lines, Nature Biotechnology, Dec. 2014, doi:10.1038/nbt.3080"""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--organism","-g",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the known fusion genes are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the known fusion genes are stored. Default is '%default'.""")

    parser.add_option("--data",
                      action = "store",
                      type = "string",
                      dest = "data_filename",
                      help="""The input CSV file containg.""")

    parser.add_option("--skip-filter-overlap",
                      action="store_true",
                      dest="skip_filter_overlap",
                      default = False,
                      help="""If set then it filters out the known fusion genes where the (i) genes are fully overlapping, or (ii) the genes are partially overlapping and are on the same strand. Default is '%default'.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    # NEW:
    #url = 'https://ndownloader.figshare.com/files/25494407'
    url = 'https://ndownloader.figshare.com/files/27902109' # 5-2021
    tmp_file = os.path.join(options.output_directory,'temp_ccle.csv')

    if options.organism.lower() == 'homo_sapiens':
        today = datetime.date.today()
        data = []
        sem = True
        if options.data_filename:
            print "Using the local file..."
        else:
            print "Downloading the known fusion genes from the DEPMAP..."
            try:
                da = urllib2.urlopen('%s' % (url,))
                file(tmp_file,'w').write(da.read())
            except:
                print >>sys.stderr, "Warning: Cannot access '%s'! The output file will be empty!" % (url,)
                sem = False


        data = set()
        if sem:
            print "Parsing..."
            # parsing
            dudu = [e.rstrip("\r\n").split(",") for e in file(tmp_file,"r").readlines() if e.rstrip("\r\n")]
            for line in dudu:
                if line and len(line)>1 and line[1].find("--") != -1:
                    gg = line[1].upper().split("--")
                    g1 = fix(gg[0])
                    g2 = fix(gg[1])
                    (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                    data.add((g1,g2))
            print " - found",len(data)

            # save version of
            txt = ['DepMap Portal: %s\n' % (today.strftime("%Y-%m-%d"),)]
            file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)

    #
            # read the gene symbols
            file_symbols = os.path.join(options.output_directory,'synonyms.txt')
            loci = symbols.generate_loci(file_symbols)

            genes = symbols.read_genes_symbols(file_symbols)

            d = []
            for (g1,g2) in sorted(data):
                if g1.upper() != g2.upper():
                    if g1 and g2:
                        ens1 = symbols.ensembl(g1.upper(),genes,loci)
                        ens2 = symbols.ensembl(g2.upper(),genes,loci)
                        if ens1 and ens2:
                            for e1 in ens1:
                                for e2 in ens2:
                                    if e1 != e2:
                                        d.append([e1,e2])

            data = ['\t'.join(sorted(line)) + '\n' for line in d]
            data = sorted(set(data))

            print "%d known fusion genes found" % (len(data),)

            if not options.skip_filter_overlap:
                ens2hugo = dict([tuple(line.rstrip('\r\n').split('\t')) for line in file(os.path.join(options.output_directory,'genes_symbols.txt'),'r').readlines() if line.rstrip('\r\n')])

                d1 = []
                overlappings = ['ensembl_fully_overlapping_genes.txt',
                                'ensembl_same_strand_overlapping_genes.txt',
                                'refseq_fully_overlapping_genes.txt',
                                'refseq_same_strand_overlapping_genes.txt',
                                "ucsc_fully_overlapping_genes.txt",
                                "ucsc_same_strand_overlapping_genes.txt",
                                'pairs_pseudogenes.txt',
                                'banned.txt',
                                'paralogs.txt']
                ensembls = set(['ensembl_fully_overlapping_genes.txt',
                                'ensembl_same_strand_overlapping_genes.txt'])
                ens = []
                for ov in overlappings:
                    p = os.path.join(options.output_directory,ov)
                    print "Parsing file:",p
                    d2 = []
                    if os.path.isfile(p):
                        d2 = sorted(set([tuple(sorted(line.rstrip('\r\n').split('\t'))) for line in file(p,'r').readlines() if line.rstrip('\r\n')]))
                        d3 = set(['\t'.join(line)+'\n' for line in d2])
                        d4 = sorted([line for line in data if line in d3])
                        d4 = [line.rstrip('\r\n').split('\t') for line in d4]
                        d4 = [line+[ens2hugo.get(line[0],'')]+[ens2hugo.get(line[1],'')] for line in d4]
                        d4 = ['\t'.join(line)+'\n' for line in d4]
                        file(os.path.join(options.output_directory,"ccle___%s"%(ov,)),'w').writelines(d4)
                    if ov in ensembls:
                        ens.extend(d2)
                    d1.extend(d2)
                d = set()
                for line in d1:
                    (a,b) = (line[0],line[1])
                    if a > b:
                        (a,b) = (b,a)
                    d.add("%s\t%s\n" % (a,b))
                ens = set(['\t'.join(line)+'\n' for line in ens])
                ensembl = [line for line in data if line in ens]
                file(os.path.join(options.output_directory,'ccle___ensembl.txt'),'w').writelines(sorted(ensembl))
                skipped = [line for line in data if line in d]
                data = [line for line in data if line not in d]
                file(os.path.join(options.output_directory,'ccle___all.txt'),'w').writelines(sorted(skipped))

                print "%d known fusion genes left after removing the overlappings" % (len(data),)

        file(os.path.join(options.output_directory,'ccle.txt'),'w').writelines(data)

        if os.path.exists(tmp_file):
            os.remove(tmp_file)

    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'ccle.txt'),'w').write('')
#
