#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the lastest known fusion genes from the Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer
<https://cgap.nci.nih.gov/Chromosomes/Mitelman>.



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
import socket
import optparse
import shutil
import urllib2
import symbols
import datetime
import zipfile

def mx(e):
    r = e
    if e[0] > e[1]:
        r = (e[1],e[0])
    return r
    
def mysplit(x):
# parses a line from a TSV where the comma is separator
# skip the commas which are between ""
    t = x.rstrip("\r\n")
    t = t.split('"')
    for i in xrange(len(t)):
        if i % 2 == 1:
            t[i] = t[i].replace(',','\t')
    t = ''.join(t)
    t = [e.replace('\t',',') for e in t.split(',')]
    return t

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known fusion genes from the Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer."""
    version = "%prog 0.15 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism","-g",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the known fusion genes are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--url","-u",
                      action="store",
                      type="string",
                      dest="url",
                      #default = 'https://storage.cloud.google.com/mitelman-data-files/prod/mitelman_db.zip',
                      default = 'https://storage.googleapis.com/mitelman-data-files/prod/mitelman_db.zip',
                      help="""The URL for Mitelman database dump file. Default is '%default'.""")

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
    tmp_file = 'temp_mitelman_db.zip'
    tmp2_file = 'temp_mitelman.dat'

    headers = { 'User-Agent' : 'Mozilla/5.0' }

    url = options.url
    #url = "https://storage.cloud.google.com/mitelman-data-files/mitelman_db.zip"
    # https://storage.cloud.google.com/mitelman-data-files/mitelman_db.zip
    data = []
    fusions = []
    file(os.path.join(options.output_directory,'mitelman.txt'),'w').write('')

    if options.organism.lower() == 'homo_sapiens':
        today = datetime.date.today()
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(['Mitelman database version: %s\n' % (today.strftime("%Y-%m-%d"),)])

        sem = True
        if url.startswith("ftp") or url.startswith("http"):
            print "Downloading the known fusion genes from the Mitelman database!"
            try:
                req = urllib2.Request(url, headers=headers)
                da1 = urllib2.urlopen(req)
                file(tmp_file,'w').write(da1.read())
            except:
                print >>sys.stderr, "Warning: Cannot access '%s'! The output file will be empty!" % (url,)
                sem = False
        else:
            shutil.copyfile(url,tmp_file)




        if sem:

            print "Parsing file..."

            ret = None
            try:
                zf = zipfile.ZipFile(tmp_file,'r')
                ret = zf.testzip()
            except:
                print >>sys.stderr, "Warning: File '%s' is a bad ZIP file! The output file will be empty!" % (url,)
                ret = 1

            if ret:
                print >>sys.stderr, "Warning: File '%s' is not a valid ZIP file! The output file will be empty!" % (url,)
            else:
            
                d = zf.read('mitelman_db/MBCA.TXT.DATA').decode('ascii').splitlines()

                if d:
                    #h = d.pop(0) # remove header
                    fusions = set()

                    f = [e.rstrip("\r\n").split("\t")[7] for e in d if e.rstrip("\r\n") and len(e.split("\t"))>7]
                    for e in f:
                        x = e.split(",")
                        for z in x:
                            if z and z.find("::") != -1:
                                uf = tuple(z.upper().replace("+","").split("::"))
                                if uf and len(uf) == 2:
                                    fusions.add(uf)
                                elif len(uf) == 3:
                                    fusions.add((uf[0],uf[1]))
                                    fusions.add((uf[1],uf[2]))
                                    fusions.add((uf[0],uf[2]))

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
                        if g1 in ("IGH","IGL","IGK","IG","TRA","TRB","TRD","TRG"):
                            g1 = g1 + "@"
                        if g2 in ("IGH","IGL","IGK","IG","TRA","TRB","TRD","TRG"):
                            g2 = g2 + "@"
                        if ( g1.upper() == g2.upper() or ((g1.endswith('@') and g2.endswith('@')) and g1.upper()[:2] == g2.upper()[:2])):
                            print "%s-%s skipped!" % (g1,g2)
                            continue
                        ens1 = symbols.ensembl(g1,genes,loci)
                        ens2 = symbols.ensembl(g2,genes,loci)

                        if ens1 and ens2:
                            for e1 in ens1:
                                for e2 in ens2:
                                    if e1 != e2 and ((e1,e2) not in banned) and ((e2,e1) not in banned):
                                        if e1 > e2:
                                            d.append([e2,e1])
                                        else:
                                            d.append([e1,e2])

                    data = ['\t'.join(sorted(line)) + '\n' for line in d]
                    data = sorted(set(data))

                    print "%d known fusion genes converted succesfully to Ensembl Gene ids!" % (len(data),)
        else:
            data = []
        file(os.path.join(options.output_directory,'mitelman.txt'),'w').writelines(data)
        file(os.path.join(options.output_directory,'mitelman_genes.txt'),'w').writelines(sorted(set(["--".join(mx(e))+"\n" for e in fusions])))

    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'mitelman.txt'),'w').write('')
    if os.path.isfile(tmp_file):
        os.remove(tmp_file)
#
