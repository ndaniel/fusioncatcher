#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the lastest known fusion genes from OncoKB repository 
<https://github.com/oncokb/oncokb> (Chakravarty et al., JCO PO 2017).



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


def is_url_valid(u):
    flag = True
    try:
        urllib2.urlopen(u)
    except urllib2.HTTPError, e:
        flag = False
    except urllib2.URLError, e:
        flag = False
    return flag



if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known fusion genes from OncoKB repository <https://github.com/oncokb/oncokb> (Chakravarty et al., JCO PO 2017)."""
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
                      default = "http://github.com",
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
    tmp_file = 'temp_oncokb.txt'

    url = '/oncokb/oncokb-public/raw/master/data/v1.XYZ/allAnnotatedVariants.txt'

    if options.organism.lower() == 'homo_sapiens':
        today = datetime.date.today()
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(['OncoKB (Chakravarty et al., JCO PO 2017): %s\n' % (today.strftime("%Y-%m-%d"),)])

        # find the latest version by starting from release 23
        print "Find latest release..."
        k = 23
        while True:
            if is_url_valid(options.server+url.replace("XYZ",str(k))):
                print " - testing:",k,"-> ok"
                k = k + 1
            else:
                print " - testing:",k,"-> not ok"
                break
        k = k - 1
        url = url.replace("XYZ",str(k))
        print " - found: 1."+str(k)

        print "Downloading the known fusion genes..."
        sem = True
        try:
            d = urllib2.urlopen('%s%s' % (options.server,url))
            file(tmp_file,'w').write(d.read())
        except:
            sem = False
            print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty!" % (options.server,url)

        monkey = set(["IGH","IGL","IGK","TRA","TRB"])

        removed = []
        if sem:
            print "Parsing..."
            fusions = set()
            for line in file(tmp_file,'r').readlines():
                li = line.rstrip('\r\n').split('\t')
                fus = li[4]
                if not fus.endswith("usion"):
                    continue
                fus = fus.upper().partition(" FUSION")[0]
                fus = fus.replace("?","-").split("-")
                if len(fus) == 2:
                    g1 = fus[0].strip()
                    g2 = fus[1].strip()
                elif len(fus) == 3:
                    if fus[0] == "HLA":
                        g1 = fus[0].strip()+"-"+fus[1].strip()
                        g2 = fus[2].strip()
                    elif fus[1] == "HLA":
                        g1 = fus[1].strip()+"-"+fus[2].strip()
                        g2 = fus[0].strip()
                    if fus[0] == "NKX2":
                        g1 = fus[0].strip()+"-"+fus[1].strip()
                        g2 = fus[2].strip()
                    elif fus[1] == "NKX2":
                        g1 = fus[1].strip()+"-"+fus[2].strip()
                        g2 = fus[0].strip()
                else:
                    continue

                if g1 in monkey:
                    g1 = g1 + "@"
                if g2 in monkey:
                    g2 = g2 + "@"

                # fix bug in OncoKB
                if g1 == "SEC16A1":
                    g1 = "SEC16A"
                if g2 == "SEC16A1":
                    g2 = "SEC16A"
                if g1 == "NKX2":
                    g1 = "NKX2-1"
                if g2 == "NKX2":
                    g2 = "NKX2-1"

                (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                print g1+"/"+g2
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
        file(os.path.join(options.output_directory,'oncokb.txt'),'w').writelines(data)
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)
    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'oncokb.txt'),'w').write('')
#
