#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the cancer associated genes from: http://www.bushmanlab.org/links/genelists


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
    description = """It downloads the cancer associated genes from: http://www.bushmanlab.org/links/genelists"""
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


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)



    url = 'http://www.bushmanlab.org/assets/doc/allOnco_Feb2017.tsv'
    tmp_file = os.path.join(options.output_directory,'temp_cancer.tsv')

    headers = {     'User-agent': 'Mozilla/5.0 (Windows; U; Windows NT 6.1; en-GB; rv:1.9.2.3) Gecko/20100401 Firefox/3.6.3',
    'Accept' : 'text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*;q=0.5',
    'Accept-Language' : 'en-gb,en;q=0.5'
    }

    file(os.path.join(options.output_directory,'cancer_genes.txt'),'w').write('')

    # save version of
    today = datetime.date.today()
    txt = ['Cancer Gene List database version: %s\n' % (today.strftime("%Y-%m-%d"),)]
    file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)


    if options.organism.lower() == 'homo_sapiens':

        data = []
        sem = True
        print "Downloading the cancer genes ..."
        try:
            req = urllib2.Request('%s' % (url,), headers=headers)
            da = urllib2.urlopen(req)
            file(tmp_file,'w').write(da.read())
        except:
            print >>sys.stderr, "Warning: Cannot access '%s'! The output file will be empty!" % (url,)
            sem = False

        data = []
        if sem:
            print "Parsing..."
            # parse the file with the known fusion genes
            data = set()
    
            d = [line.upper().rstrip("\r\n").split("\t")[0] for line in file(tmp_file,"r")]
            d.pop(0) # remove the header
            data = set(d)

            print "%d known genes found (using gene symbols)" % (len(data),)

            # read the gene symbols
            file_symbols = os.path.join(options.output_directory,'synonyms.txt')
            loci = symbols.generate_loci(file_symbols)

            genes = symbols.read_genes_symbols(file_symbols)

            d = []
            for g in data:
                ens = symbols.ensembl(g.upper(),genes,loci)
                if ens:
                    d.extend(ens)

            data = [line + '\n' for line in d]
            data = sorted(set(data))

            print "%d known genes found (after conversion to Ensembl ids)" % (len(data),)



        file(os.path.join(options.output_directory,'cancer_genes.txt'),'w').writelines(data)

        if os.path.exists(tmp_file):
            os.remove(tmp_file)

    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'cancer_genes.txt'),'w').write('')
#
