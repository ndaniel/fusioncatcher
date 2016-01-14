#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available sequences of the viral genomes.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2016 Daniel Nicorici

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
    version = "%prog 0.10 beta"

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

    #ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz

    #
    # FTP => ftp://ftp.ncbi.nih.gov/genomes/
    # web => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Genome
    #
    # FTP => ftp://ftp.ncbi.nih.gov/refseq/
    # web => http://www.ncbi.nlm.nih.gov/LocusLink/refseq.html


    url = 'genomes/Viruses/'
    version = ""
    #mypath = os.path.join(options.output_directory,'viruses_index')
    mypath = options.output_directory
    if not os.path.isdir(mypath):
        os.makedirs(mypath)

    print "Downloading the viral genomes from NCBI server!"
    try:
        ftp = ftplib.FTP(options.server)
        print ftp.login()
        ftp.cwd(url)

        list_files = ftp.nlst()

        list_files = [el for el in list_files if el.lower() == 'all.fna.tar.gz']

        if len(list_files) == 1:
            new_files = []
            for filename in list_files:
                print "Downloading: %s/%s/%s" % (options.server,url,filename)
                nf = os.path.join(mypath,filename)
                new_files.append(nf)
                fid = open(nf,'wb')
                ftp.retrbinary("RETR " + filename, fid.write)
                modified_time = ftp.sendcmd('MDTM ' + filename)
                version = datetime.datetime.strptime(modified_time[4:], "%Y%m%d%H%M%S").strftime("%Y-%m-%d")
                fid.close()
        else:
            print "ERROR: File not 'all.fna.tar.gz' found!"
            sys.exit(1)
        ftp.close()
    except ftplib.all_errors, e:
        print 'FTP Error = ' + str(e)
        sys.exit(1)
    except Exception, e:
        print "Error: Generic exception!",str(e)
        sys.exit(1)

    file_output = os.path.join(mypath,'viruses.fa')
    fout = file(file_output,'w')
    for filename in new_files:
        try:
            f = tarfile.open(filename, 'r:gz')
        except:
            print "ERROR: File '%s' is corrupt!" % (filename,)
            sys.exit(1)

        try:
            members = f.getmembers()
        except:
            print "ERROR: File '%s' is corrupt!" % (filename,)
            sys.exit(1)

        for member in members:
            if member.name.lower().endswith('.fna') or member.name.lower().endswith('.fa'):
                try:
                    fin = f.extractfile(member)
                except:
                    print "ERROR: File '%s' is corrupt!" % (filename,)
                    sys.exit(1)
                while True:
                    data = fin.read(BLOCK_SIZE)
                    if not data:
                        break
                    fout.write(data)
                fin.close()
            elif member.name.lower().endswith('.fna.tgz') or member.name.lower().endswith('.fna.tar.gz') or member.name.lower().endswith('.fa.tgz') or member.name.lower().endswith('.fa.tar.gz'):
                try:
                    fin = f.extractfile(member)
                except:
                    print "ERROR: File '%s' is corrupt!" % (filename,)
                    sys.exit(1)

                compressedstream = StringIO.StringIO(fin.read())

                try:
                    f2 = tarfile.open(fileobj = compressedstream, mode = 'r:gz')
                except:
                    print "ERROR: File '%s' is corrupt!" % (filename,)
                    sys.exit(1)

                try:
                    members2 = f2.getmembers()
                except:
                    print "ERROR: File '%s' is corrupt!" % (filename,)
                    sys.exit(1)
                for member2 in members2:
                    try:
                        fin2 = f2.extractfile(member2)
                    except:
                        print "ERROR: File '%s' is corrupt!" % (filename,)
                        sys.exit(1)
                    while True:
                        data = fin2.read(BLOCK_SIZE)
                        if not data:
                            break
                        fout.write(data)
                fin2.close()
                fin.close()


        #print "--------------"
        #print f.getnames()
        #f.extractall(mypath)
        f.close()
        os.remove(filename)
    fout.close()

    txt =["NCBI Viral Genomes version: %s\n" % (version,)]
    file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)
    #
