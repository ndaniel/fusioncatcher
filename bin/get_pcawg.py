#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the lastest known fusion genes from PCWAG project
<https://dcc.icgc.org/releases/PCAWG/transcriptome/fusion>.


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



import os
import sys
import socket
import optparse
import shutil
import urllib2
import gzip
import symbols

def check(s):
    r = ''
    if isinstance(s, basestring):
        r = s
    return r
        


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known fusion genes from PCWAG project."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the known fusion genes are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the known fusion genes are stored. Default is '%default'.""")

    parser.add_option("--server",
                      action="store",
                      type="string",
                      dest="server",
                      default = "https://dcc.icgc.org",
                      help="""The ChimerDB 3.0 server from where the known fusion genes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    file(os.path.join(options.output_directory,'pcawg.txt'),'w').write('')


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)


    #http://ercsb.ewha.ac.kr/FusionGene/document/PO_down.xls
    url = '/api/v1/download?fn=/PCAWG/transcriptome/fusion/gene.fusions.V1.tsv.gz'

    headers = { 'User-Agent' : 'Mozilla/5.0' }

    if options.organism.lower() == 'homo_sapiens':

        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(['PCAWG database version: 1.0\n'])

        print "Downloading the known fusion genes from PCAWG database!"
        data = []

        tmp_file = os.path.join(options.output_directory,'temp_fusions.txt.gz')
        sem = True
        try:
            req = urllib2.Request('%s%s' % (options.server,url), headers=headers)
            d = urllib2.urlopen(req)
            file(tmp_file,'w').write(d.read())
        except:
            sem = False
            print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty!" % (options.server,url)

        if sem:
            print "Parsing...",tmp_file
            # parse the ChimerDB file with the known fusion genes
            fusions = set()
            fuse = set()

            # parse the PubMed sheet
            fi = [e.rstrip("\r\n").split("\t") for e in gzip.open(tmp_file,"r").readlines() if e.rstrip("\r\n")]
            fi.pop(0)
            for row in fi:
                tg = row[0].split("->")
                g1 = tg[0]
                g2 = tg[1]
                (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)

                fusions.add((g1,g2))

                fg1 = row[4].partition(".")[0]
                fg2 = row[5].partition(".")[0]
                (fg1,fg2) = (fg2,fg1) if fg2 < fg1 else (fg1,fg2)
                fuse.add((fg1,fg2))

            fusions = list(fusions)
            print "%d known gene fusions found!" % (len(fusions),)

            # read the gene symbols
            file_symbols = os.path.join(options.output_directory,'synonyms.txt')
            genes = symbols.read_genes_symbols(file_symbols)

            banned = set()
            loci = symbols.generate_loci(file_symbols)
            #for v in symbols.locus.values():
            for v in loci.values():
                if v:
                    n = len(v)
                    if n > 1:
                        for i in xrange(n-1):
                            for j in xrange(i+1,n):
                                if v[i].upper() != v[j].upper():
                                    ens1 = symbols.ensembl(v[i].upper(),genes,loci)
                                    ens2 = symbols.ensembl(v[j].upper(),genes,loci)
                                    if ens1 and ens2:
                                        for e1 in ens1:
                                            for e2 in ens2:
                                                if e1 != e2:
                                                    (e1,e2) = (e2,e1) if e2 < e1 else (e1,e2)
                                                    banned.add((e1,e2))


            d = []
            for (g1,g2) in fusions:
                if ( g1.upper() == g2.upper() or ((g1.endswith('@') and g2.endswith('@')) and g1.upper()[:2] == g2.upper()[:2])):
                    print "%s-%s skipped!" % (g1,g2)
                    continue
                ens1 = symbols.ensembl(g1,genes,loci)
                ens2 = symbols.ensembl(g2,genes,loci)

                if ens1 and ens2:
                    for e1 in ens1:
                        for e2 in ens2:
                            if e1 != e2 and ((e1,e2) not in banned) and ((e2,e1) not in banned):
                                d.append([e1,e2])
            for (ens1,ens2) in fuse:
                if ens1 and ens2:
                    if ens1 != ens2 and ((ens1,ens2) not in banned) and ((ens2,ens1) not in banned):
                        d.append([ens1,ens2])

            data = ['\t'.join(sorted(line)) + '\n' for line in d]
            data = sorted(set(data))

            print "%d known fusion genes converted succesfully to Ensembl Gene ids!" % (len(data),)

        file(os.path.join(options.output_directory,"pcawg.txt"),'w').writelines(data)
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)


#
