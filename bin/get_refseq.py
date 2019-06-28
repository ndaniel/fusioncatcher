#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available genes annotations (chromosomal positions and gene symbols) from RefSeq NCBI database.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

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
                      default = "hgdownload.cse.ucsc.edu",
                      help="""The UCSC server from where the RefSeq gene annotations are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
    #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.sql

    #Files:
    #
    # refFlat.txt.gz
    # refFlat.sql
    #

    org = ''
    if options.organism.lower() == 'homo_sapiens':
        org = 'hg19'
        # find genome information
        d = [line for line in file(os.path.join(options.output_directory,'version.txt'),'r') if line.lower().startswith('genome version') ]
        if d:
            if d[0].lower().find('grch38') !=-1:
                org = 'hg38'
    elif options.organism.lower() == 'rattus_norvegicus':
        org = 'rn6'
    elif options.organism.lower() == 'mus_musculus':
        org = 'mm10'
    elif options.organism.lower() == 'canis_familiaris':
        org = 'canFam3'

    if org:
        url = 'goldenPath/%s/database' % (org,)

        files = {'refFlat':'refFlat.txt.gz',
                 'sql':'refFlat.sql'
        }


    if org:
        url = 'goldenPath/%s/database' % (org,)

        files = {'refFlat':'refFlat.txt.gz',
                 'sql':'refFlat.sql'
        }

        print "Downloading the SQL files of organism '%s' from UCSC!" % (options.organism.lower(),)
        try:
            ftp = ftplib.FTP(options.server)
            print ftp.login()
            ftp.cwd(url)

            for k,filename in files.items():
                print "Downloading: %s%s/%s" % (options.server,url,filename)
                nf = os.path.join(options.output_directory,filename)
                files[k] = nf
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
        for k,filename in files.items():
            if filename.endswith('.gz'):
                f = gzip.open(filename, 'rb')
                file_content = f.read()
                f.close()
                fn = filename[:-3]
                files[k] = fn
                fod = open(fn,'wb')
                fod.write(file_content)
                fod.close()
                os.remove(filename)

        print "Parsing the MySQL files..."
        # parse the SQL file for column numbers
        table = dict()
        last_name = None
        version = ''
        for key,filename in files.items():
            if filename.endswith('.sql'):
                mysql = [line.rstrip('\r\n') for line in file(filename,'r').readlines() if line.rstrip('\r\n')]

                # get version UCSC
                x = 'dump completed on '
                v = [line.lower().strip().partition(x)[2].partition(' ')[0] for line in mysql if line.lower().find(x) != -1]
                if v and not version:
                    version = v.pop(0)

                idx = -1
                for line in mysql:
                    if (not line) or line.strip().startswith('/*') or line.strip().startswith('--') or line.lower().strip().startswith('drop table') or line.lower().strip().startswith('key'):
                        continue
                    parts = line.strip().split("`")
                    n = len(parts)
                    if parts[0].find('CREATE TABLE') != -1:
                        last_name = parts[1]
                        if last_name not in table:
                            table[last_name] = {}
                            idx = -1
                    elif n > 2 and parts[0].find('KEY') == -1 and parts[0].find('ENGINE') == -1 and not parts[0].strip().startswith(")"):
                        idx = idx + 1
                        col = parts[1]
                        table[last_name][col] = idx

        # parse "refFlat.txt" file
        gene = [line.rstrip('\r\n').split('\t') for line in file(files['refFlat'],'r').readlines() if line.rstrip('\r\n')]
        gene_symbol = table['refFlat']['geneName']
        gene_chrom = table['refFlat']['chrom']
        gene_strand = table['refFlat']['strand']
        gene_start = table['refFlat']['txStart']
        gene_end = table['refFlat']['txEnd']
        gene = set([(line[gene_symbol].strip().upper().replace(' ','_'),
                     line[gene_end],
                     line[gene_start],
                     line[gene_strand],
                     line[gene_chrom]
                     ) for line in gene if line[gene_symbol]])
        print "%d genes found in 'refFlat.txt'!" % (len(gene),)

        data = []
        for x in gene:
            s = int(x[2])
            e = int(x[1])
            if s > e:
                (s,e) = (e,s)
            data.append([x[0],str(e),str(s),x[3],x[4]])
        data = ['\t'.join(line)+'\n' for line in data]
        file(os.path.join(options.output_directory,'refseq_genes.txt'),'w').writelines(data)
        txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
        file(os.path.join(options.output_directory,'refseq_genes_header.txt'),'w').write(txt)

        txt = ["RefSeq NCBI database version (downloaded thru UCSC database; %s): %s\n" % (org,version)]
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)

        #
        #
        for f in files.values():
            os.remove(f)
    else:
        data = []
        file(os.path.join(options.output_directory,'refseq_genes.txt'),'w').writelines(data)
        txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
        file(os.path.join(options.output_directory,'refseq_genes_header.txt'),'w').write(txt)
    #
