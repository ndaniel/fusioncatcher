#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the lastest known duplicated genes from the DGD database
<http://dgd.genouest.org/>.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2014 Daniel Nicorici

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
import itertools
import gzip

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known fusion genes from the duplicat genes DGD database."""
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
                      default = "http://dgd.genouest.org",
                      help="""The TICdb server from where the known fusion genes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)
    tmp_file = 'temp_dgd.txt.gz'
    # columns are:
    # chr - 0
    # group_id - 1
    # NB_Genes - 2
    # start - 3
    # end - 4
    # strand - 5
    # ENS_ID - 6
    # Name - 7
    # Description - 8

    # http://dgd.genouest.org/listRegion/rawgz/homo_sapiens/all%3A0..x/
    # http://dgd.genouest.org/listRegion/rawgz/canis_familiaris/all%3A0..x/
    # http://dgd.genouest.org/listRegion/rawgz/mus_musculus/all%3A0..x
    # http://dgd.genouest.org/listRegion/rawgz/rattus_norvegicus/all%3A0..x/
    url = '/listRegion/rawgz/homo_sapiens/all%3A0..x/'

    urls = {'homo_sapiens': '/listRegion/rawgz/homo_sapiens/all%3A0..x/',
            'rattus_norvegicus':'/listRegion/rawgz/rattus_norvegicus/all%3A0..x/',
            'mus_musculus':'/listRegion/rawgz/mus_musculus/all%3A0..x',
            'canis_familiaris':'/listRegion/rawgz/canis_familiaris/all%3A0..x/',
            'bos_taurus':'/listRegion/rawgz/bos_taurus/all%3A0..x/',
            'danio_rerio':'/listRegion/rawgz/danio_rerio/all%3A0..x/',
            'equus_caballus':'/listRegion/rawgz/equus_caballus/all%3A0..x/',
            'gallus_gallus':'/listRegion/rawgz/gallus_gallus/all%3A0..x/',
            'sus_scrofa':'/listRegion/rawgz/sus_scrofa/all%3A0..x/'
           }

    if options.organism.lower() in urls:
        url = urls[options.organism.lower()]
        today = datetime.date.today()
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(['DGD database version: %s\n' % (today.strftime("%Y-%m-%d"),)])

        print "Downloading the known fusion genes from DGD database!"
        sem = True
        try:
            d = urllib2.urlopen('%s%s' % (options.server,url))
            file(tmp_file,'w').write(d.read())
        except:
            sem = False
            print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty!" % (options.server,url)

        if sem:
            try:
                f = gzip.open(tmp_file, 'rb')
                file_content = f.read()
                f.close()
            except IOError:
                print >>sys.stderr,"Warning: The downloaded file is not a gzipped file!"
                os.remove(tmp_file)
                sem = False
                sys.exit(1)


            fn = tmp_file[:-3] # remove '.gz'
            fod = open(fn,'wb')
            fod.write(file_content)
            fod.close()
            os.remove(tmp_file)
            tmp_file = fn

            print "Parsing..."
            fusions = set()
            e_ids= set()
            bucket = set()
            bucket_e = set()
            last = None
            d = file(tmp_file,'r').readlines()
            d.pop(0)
            for line in d:
                if not line:
                    continue
                line = line.rstrip('\r\n').split('\t')
                group = line[1]
                gene = line[7].upper()
                e = line[6]
                if last != group:
                    if bucket:

                        for (g1,g2) in itertools.combinations(list(bucket),2):
                            (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                            fusions.add((g1,g2))
                        for (g1,g2) in itertools.combinations(list(bucket_e),2):
                            (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                            if g1.startswith('ENS') and g2.startswith('ENS'):
                                e_ids.add((g1,g2))

                        bucket = set()
                        bucket_e = set()
                    last = group
                bucket.add(gene)
                bucket_e.add(e)

            if bucket:
                for (g1,g2) in itertools.combinations(list(bucket),2):
                    (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                    fusions.add((g1,g2))
                for (g1,g2) in itertools.combinations(list(bucket_e),2):
                    (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                    if g1.startswith('ENS') and g2.startswith('ENS'):
                        e_ids.add((g1,g2))


            #
            fusions = list(fusions)
            print "%d known duplicated genes found (using gene symbols)!" % (len(fusions),)
            print "%d known duplicated genes found (using Ensemb gene ids)!" % (len(e_ids),)

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
            data.extend('\t'.join(sorted(line)) + '\n' for line in list(e_ids))
            data = sorted(set(data))

            print "%d known duplicated genes converted succesfully to Ensembl Gene ids!" % (len(data),)
        else:
            data = []
        file(os.path.join(options.output_directory,'dgd.txt'),'w').writelines(data)
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)
    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'dgd.txt'),'w').write('')
#
