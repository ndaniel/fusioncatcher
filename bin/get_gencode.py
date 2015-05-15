#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available genes annotations (chromosomal positions and gene symbols) from RefSeq NCBI database.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

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
import gzip
import socket
import optparse
import concatenate
import shutil
import StringIO


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest available genes annotations (chromosomal positions and gene symbols) from RefSeq NCBI database."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

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
                      default = "ftp.sanger.ac.uk",
                      help="""The Gencode server from where the gene annotations are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    #ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz



    org = ''
    if options.organism.lower() == 'homo_sapiens':
        org = 'Gencode_human'
    elif options.organism.lower() == 'mus_musculus':
        org = 'Gencode_mouse'



    if org:

        url = 'pub/gencode/%s/' % (org,)
        print "Downloading the GTF file of organism '%s' from Gencode!" % (options.organism.lower(),)
        version = ''
        nf = None
        filename= None
        try:
            ftp = ftplib.FTP(options.server)
            print ftp.login()
            ftp.cwd(url)

            list_files = ftp.nlst()

            list_files = [el.replace("release_","") for el in list_files if el.lower().startswith('release_')]
            if options.organism.lower() == 'homo_sapiens':
                list_files = sorted([int(el) for el in list_files if el.isdigit()])
            else:
                list_files = sorted(list_files)
            version = str(list_files[-1])
            last = "release_"+version
            filename = "gencode.v%s.annotation.gtf.gz" % (version,)
            url = "%s%s" % (url,last)
            ftp.cwd(last)
            print "cd ",last
            print "Downloading: %s%s/%s" % (options.server,url,filename)
            nf = os.path.join(options.output_directory,filename)
            fid = open(nf,'wb')
            ftp.retrbinary("RETR " + filename, fid.write)
            fid.close()

            ftp.close()
        except ftplib.all_errors, e:
            print 'FTP Error = ' + str(e)
            sys.exit(1)
        except Exception, e:
            print "Error: Generic exception!",str(e)
            sys.exit(1)

        print "Decompressing files ..."
        if filename.endswith('.gz'):
            f = gzip.open(nf, 'rb')
            file_content = f.read()
            f.close()
            f = nf[:-3]
            fod = open(f,'wb')
            fod.write(file_content)
            fod.close()
            os.remove(nf)
            nf = f

        print "Parsing the GTF file..."
        d = [line.split("\t") for line in file(nf,'r').readlines() if (not line.startswith("#")) and line]
        d = [(line[0],line[3],line[4],line[6],line[8].partition('gene_name "')[2].partition('"')[0]) for line in d if line[2] == 'gene']
        print "%d genes found!" % (len(d),)

        data = []
        for x in d:
            s = int(x[2])
            e = int(x[1])
            if s > e:
                (s,e) = (e,s)
            data.append([x[4],str(e),str(s),x[3],x[0]])
        data = ['\t'.join(line)+'\n' for line in data]
        file(os.path.join(options.output_directory,'gencode_genes.txt'),'w').writelines(data)
        txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
        file(os.path.join(options.output_directory,'gencode_genes_header.txt'),'w').write(txt)

        txt = ["Gencode database version: %s\n" % (version,)]
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)

        #
        #
        os.remove(nf)
    else:
        data = []
        file(os.path.join(options.output_directory,'gencode_genes.txt'),'w').writelines(data)
        txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
        file(os.path.join(options.output_directory,'gencode_genes_header.txt'),'w').write(txt)
    #
