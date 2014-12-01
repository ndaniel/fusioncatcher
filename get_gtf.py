#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the annotation in GTF foramt for a given organism from Ensembl.


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
#ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna
#
#Files:
#
#Homo_sapiens.NCBI36.*.dna.chromosome.**.fa.gz
#
#where:
#* is 54 (ensembl version)
#** is chromosome name e.g. 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT

#Hint: use ftplib and gzip

import os
import sys
import ftplib
import gzip
import socket
import optparse
import concatenate
import shutil


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the annotation in GTF format for a given organism from Ensembl."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the chromosomes are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

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
                      default = "ftp.ensembl.org",
                      help="""The Ensembl server from where the chromosomes are downloaded. Default is '%default'.""")

    parser.add_option("--server-path",
                      action="store",
                      type="string",
                      dest="server_path",
                      default = "/pub/current_gtf/", #ftp://ftp.ensembl.org/pub/release-75/gtf/
                      help="""The path of Ensembl server from where the data is downloaded. Default is '%default'.""")
                      
    parser.add_option("--filter-chrom",
                      action="store_true",
                      dest="filter",
                      default = False,
                      help="""All chromsomes which are not 1..99,X,Y,MT, and UN will be filtered out. Default is '%default'.""")


    #command line parsing
    (options, args) = parser.parse_args()

    #
    # validate options
    #
    if (
        (not options.output_directory)
        ):
        parser.print_help()
        sys.exit(1)

    #download file Homo_sapiens.GRCh37.56.gtf.gz
    try:
        ftp = ftplib.FTP(options.server)
        print ftp.login()
        ftp.cwd(options.server_path+options.organism.lower())

        list_files = ftp.nlst()

        list_files = [el for el in list_files if el.lower().startswith(options.organism.lower()) and el.lower().endswith('.gtf.gz') ]
        if len(list_files)!=1:
            print "Too many files or too few were found!"
            print list_files
            sys.exit(1)
        filename = list_files[0]


        fid=open(os.path.join(options.output_directory,filename),'wb')
        ftp.retrbinary("RETR " + filename, fid.write)
        fid.close()

        ftp.close()
    except ftplib.all_errors, e:
        print 'FTP Error = ' + str(e.code)
        sys.exit(1)
    except Exception, e:
        print "Error: Generic exception!",str(e)
        sys.exit(1)


    f=gzip.open(os.path.join(options.output_directory,filename), 'rb')
    file_content=f.read()
    f.close()
    temp_file = os.path.join(options.output_directory,'temp_organism.gtf')
    fod=open(temp_file,'wb')
    fod.write(file_content)
    fod.close()
    os.remove(os.path.join(options.output_directory,filename))

    #keep only the chromosomes 1, 2, 3, ... 99, MT,
    allchr = [str(i) for i in range(1,100)]
    chromosomes = set(allchr+['X','Y','MT','UN'])
    #
    if options.filter:
        data = [line for line in file(temp_file,'r').readlines() if line.startswith("#") or (line.split("\t",1)[0].upper() in chromosomes)]
        file(os.path.join(options.output_directory,"organism.gtf"),"w").writelines(data)
        os.remove(temp_file)
    else:
        os.rename(temp_file,os.path.join(options.output_directory,"organism.gtf"))

    #
