#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available sequences of the viral genomes.



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
import ftplib
import tarfile
import socket
import optparse
import shutil
import datetime
import StringIO

BLOCK_SIZE = 10 ** 8

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest available sequences of the viral genomes."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the chromosomes are stored. Default is '%default'.""")

    parser.add_option("--server",
                      action="store",
                      type="string",
                      dest="server",
                      default = "ftp.ncbi.nlm.nih.gov",
                      help="""The NCBI server from where the viral genomes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)


    #ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
    #ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz


    #ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz
    #
    # FTP => ftp://ftp.ncbi.nih.gov/genomes/
    # web => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Genome
    #
    # FTP => ftp://ftp.ncbi.nih.gov/refseq/
    # web => http://www.ncbi.nlm.nih.gov/LocusLink/refseq.html


    url = 'refseq/release/viral/'
    version = ""
    #mypath = os.path.join(options.output_directory,'viruses_index')
    mypath = options.output_directory
    if not os.path.isdir(mypath):
        os.makedirs(mypath)

    list_files = []
    downloaded = []
    print "Downloading the viral genomes from NCBI server!"
    try:
        ftp = ftplib.FTP(options.server)
        print ftp.login()
        ftp.cwd(url)

        list_files = ftp.nlst()

        list_files = [el for el in list_files if el.lower().endswith('.genomic.fna.gz')]

        if list_files:
            for filename in list_files:
                print "Downloading: %s/%s/%s" % (options.server,url,filename)
                nf = os.path.join(mypath,filename)
                downloaded.append(nf)
                fid = open(nf,'wb')
                ftp.retrbinary("RETR " + filename, fid.write)
                modified_time = ftp.sendcmd('MDTM ' + filename)
                version = datetime.datetime.strptime(modified_time[4:], "%Y%m%d%H%M%S").strftime("%Y-%m-%d")
                fid.close()
        else:
            print >>sys.stderr,"ERROR: File(s) not '.genomic.fna.gz' found!"
            sys.exit(1)
        ftp.close()
    except ftplib.all_errors, e:
        print 'FTP Error = ' + str(e)
        sys.exit(1)
    except Exception, e:
        print "Error: Generic exception!",str(e)
        sys.exit(1)

    file_output = os.path.join(mypath,'viruses.fa')
    if downloaded:
        r = os.system('zcat %s > %s' % (' '.join(downloaded),file_output) )
        if not r:
            for e in downloaded:
                os.remove(e)
        else:
            print >>sys.stderr,"ERROR: somewthing wrong with the .genomic.fna.gz files!"
            sys.exit(1)
         
    else:        
        file(file_output,'w').write('>empty\nAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAACCCCCCCCCCCGGGGGGGTTTTTTTTTT\n')

    txt =["NCBI Viral Genomes version: %s\n" % (version,)]
    file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)
    #
