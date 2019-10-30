#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest known fusion genes from COSMIC database
<http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>.



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

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known fusion genes from COSMIC database."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

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

    parser.add_option("--server","-s",
                      action="store",
                      type="string",
                      dest="server",
                      default = "http://cancer.sanger.ac.uk",
                      help="""The COSMIC server from where the known fusion genes are downloaded. Default is '%default'.""")

    parser.add_option("--data",
                      action="store",
                      type="string",
                      dest="data_filename",
                      help="""The input TSV.GZ file containg the data from the COSMIC database. It should be used when the COSMIC website cannot be reached. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    # OLD
    #ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicFusionExport_v63_300113.tsv.gz
    #ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/ # not working anymore
    #
    # NEW:
    # http://cancer.sanger.ac.uk/files/cosmic/v69/CosmicFusionExport.tsv.gz
    # http://cancer.sanger.ac.uk/files/cosmic/current_release/CosmicFusionExport.tsv.gz
    #


    url = '/files/cosmic/current_release/CosmicFusionExport.tsv.gz'
    tmp_file_gz = os.path.join(options.output_directory,'temp_cosmic.tsv.gz')

    headers = { 'User-Agent' : 'Mozilla/5.0' }

    if options.organism.lower() == 'homo_sapiens':
        today = datetime.date.today()
        data = []
        sem = True
        if options.data_filename:
            tmp_file_gz = options.data_filename
            print "Using the local COSMIC database file..."
        else:
            print "Downloading the known fusion genes from the COSMIC database..."
            try:
                req = urllib2.Request('%s%s' % (options.server,url), headers=headers)
                da1 = urllib2.urlopen(req)
                file(tmp_file_gz,'w').write(da1.read())
            except:
                print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty! Hint: Please, login to http://cancer.sanger.ac.uk using in your browser with your login and password and then re-run this script." % (options.server,url)
                sem = False

        if sem:
            try:
                f = gzip.open(tmp_file_gz, 'rb')
                file_content = f.read()
                f.close()
            except IOError:
                print >>sys.stderr,"Warning: The downloaded file is not a gzipped file! Hint: Please, login to http://cancer.sanger.ac.uk using in your browser with your login and password and then re-run this script."
                if not options.data_filename:
                    os.remove(tmp_file_gz)
                sem = False

        if sem:
            fn = tmp_file_gz[:-3] # remove '.gz'
            fod = open(fn,'wb')
            fod.write(file_content)
            fod.close()
            if not options.data_filename:
                os.remove(tmp_file_gz)
            tmp_file = fn

            # save version of COSMIC
            txt = ['COSMIC database version: %s\n' % (today.strftime("%Y-%m-%d"),)]
            file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)

            # parse the COSMIC file with the known fusion genes
            data = [line.upper().rstrip('\r\n').split('\t') for line in file(tmp_file,'r') if line.rstrip('\r\n')]
            os.remove(tmp_file)
            header = dict([(v.lower(),k) for (k,v) in enumerate(data.pop(0))])
            #h = header['fusion description']
            h = header['translocation name']
            # 'Fusion description' column contains:
            #   TMPRSS2{NM_005656.2}:r.1_142_ERG{NM_004449.3}:r.444_3097
            #   FGFR1{NM_023110.2}:r.1_2919_ZNF703{NM_025069.1}:r.441_2174
            def clean(s):
                # given r.1_2919_ZNF703{NM_025069.1}
                # it gives ZNF703---NM_025069.1
                s = s[:-1]
                s = s.split('{')
                sb = s[-1]
                sa = s[0].split('_')[-1]
                return "%s---%s" % (sa,sb)

            data = [list(set([clean(el) for el in line[h].split(':') if el.endswith('}')])) for line in data if line and len(line)>= h and line[h]]
            fusions = set()
            for line in data:
                n = len(line)
                if n == 1:
                    pass
                elif n == 2:
                    a = line[0].split('---')[0]
                    b = line[1].split('---')[0]
                    (a,b) = (b,a) if a > b else (a,b)
                    if a.lower() != b.lower():
                        fusions.add((a,b))
                else:
                    for i in xrange(n-1):
                        a = line[i].split('---')
                        for j in xrange(i+1,n):
                            b = line[j].split('---')
                            if a[1] != b[1]:
                                aa = a[0]
                                bb = b[0]
                                (aa,bb) = (bb,aa) if aa > bb else (aa,bb)
                                if aa.lower() != bb.lower():
                                    fusions.add((aa,bb))

    #
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

            print "%d known fusion genes found in COSMIC database" % (len(data),)

        file(os.path.join(options.output_directory,'cosmic.txt'),'w').writelines(data)

    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'cosmic.txt'),'w').write('')
#
