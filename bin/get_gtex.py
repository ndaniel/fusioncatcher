#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the lastest known candidate fusion genes found in silico in 
The Genotype-Tissue Expression Project (GTEx) 
<http://commonfund.nih.gov/GTEx> thru FusionAnnotator."



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
import symbols
import datetime

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known candidate fusion genes found in silico in The Genotype-Tissue Expression Project (GTEx) <http://commonfund.nih.gov/GTEx> thru FusionAnnotator."""
    version = "%prog 0.20 beta"

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
                      default = "http://raw.githubusercontent.com",
                      help="""The GETx server (thru FusionAnnotator) from where the known fusion genes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)
    tmp_file = 'temp_gtex.txt'

    url1 = '/ndaniel/FusionAnnotator/master/Hg19_CTAT_fusion_annotator_lib/GTEx_v2.txt'
    url2 = '/ndaniel/FusionAnnotator/master/Hg19_CTAT_fusion_annotator_lib/Stransky2014_GTEx_normals_pairs.txt'
    url3 = '/ndaniel/FusionAnnotator/master/Hg19_CTAT_fusion_annotator_lib/GTEx_Recurrent_Blacklist_July222016.txt'

    if options.organism.lower() == 'homo_sapiens':
        today = datetime.date.today()
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(['GTEx database version (thru FusionAnnotator): %s\n' % (today.strftime("%Y-%m-%d"),)])

        print "Downloading the known fusion genes from GTEx database thru FusionAnnotator!"
        sem = True
        try:
            d1 = urllib2.urlopen('%s%s' % (options.server,url1))
            file(tmp_file,'w').write(d1.read())
            file(tmp_file,'a').write('\n')
            d2 = urllib2.urlopen('%s%s' % (options.server,url2))
            file(tmp_file,'a').write(d2.read())
#            d3 = urllib2.urlopen('%s%s' % (options.server,url3))
#            file(tmp_file,'a').write(d3.read())
        except:
            sem = False
            print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty!" % (options.server,url)

        removed = []
        if sem:
            print "Parsing..."
            fusions = set()
            for line in file(tmp_file,'r').readlines():
                li = line.rstrip('\r\n').split('\t')
                if (not li) or li[0].startswith('#'):
                    continue
                counts = 0
                if len(li) > 1:
                    counts = int(li[1])
                li = li[0].split('--')
                if len(li) != 2:
                    continue

                if counts < 5:
                    removed.append('%s--%s\t%d\n' % (li[0],li[1],counts))
                    continue
                u1 = li[0].upper().split('|')
                u2 = li[1].upper().split('|')
                for g1 in u1:
                    for g2 in u2:
                        (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                        if (g2 == "TMPRSS2" and g1 == "ERG") or (g2 == "BCR" and g1 == "ABL1"):
                            pass
                        else:
                            fusions.add((g1,g2))

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


            data = ['\t'.join(sorted(line)) + '\n' for line in d]
            data = sorted(set(data))

            print "%d known fusion genes converted succesfully to Ensembl Gene ids!" % (len(data),)
        else:
            data = []
        file(os.path.join(options.output_directory,'gtex.txt'),'w').writelines(data)
        file(os.path.join(options.output_directory,'gtex_removed.txt'),'w').writelines(removed)
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)
    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'gtex.txt'),'w').write('')
#
