#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available genes annotations (chromosomal positions and gene symbols) from RefSeq NCBI database.



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
import gzip
import socket
import optparse
#import concatenate
import shutil
import StringIO
import time


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
                      default = "ftp.ebi.ac.uk",
                      help="""The Gencode server from where the gene annotations are downloaded. Default is '%default'.""")

    parser.add_option("--release",
                      action="store",
                      type="int",
                      dest="release",
                      default = 0,
                      help="""Pin to a specific GENCODE release number (e.g. 48). Default 0 = use latest available.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    #ftp://ftp.sanger.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz



    org = ''
    if options.organism.lower() == 'homo_sapiens':
        org = 'Gencode_human'
    elif options.organism.lower() == 'mus_musculus':
        org = 'Gencode_mouse'



    if org:

        base_url = 'pub/databases/gencode/%s/' % (org,)
        print "Downloading the GTF file of organism '%s' from Gencode!" % (options.organism.lower(),)
        version = ''
        nf = None
        filename= None
        MAX_RETRIES = 5
        for attempt in xrange(MAX_RETRIES):
            url = base_url
            try:
                ftp = ftplib.FTP(options.server)
                print ftp.login()
                ftp.cwd(url)

                if options.release:
                    version = str(options.release)
                    if options.organism.lower() == 'mus_musculus':
                        last = "release_M"+version
                        filename = "gencode.vM%s.annotation.gtf.gz" % (version,)
                    else:
                        last = "release_"+version
                        filename = "gencode.v%s.annotation.gtf.gz" % (version,)
                else:
                    list_files = ftp.nlst()
                    list_files = [el.replace("release_M","").replace("release_","") for el in list_files if el.lower().startswith('release_')]
                    if options.organism.lower() == 'homo_sapiens':
                        list_files = sorted([int(el) for el in list_files if el.isdigit() ]) # and el != '24' new "and el!='24'" because is something wrong with release_24
                        version = str(list_files[-1])
                        last = "release_"+version
                        filename = "gencode.v%s.annotation.gtf.gz" % (version,)
                    elif options.organism.lower() == 'mus_musculus':
                        list_files = sorted([int(el) for el in list_files if el.isdigit() ])
                        version = str(list_files[-1])
                        last = "release_M"+version
                        filename = "gencode.vM%s.annotation.gtf.gz" % (version,)
                    else:
                        list_files = sorted(list_files)
                        version = str(list_files[-1])
                        last = "release_"+version
                        filename = "gencode.v%s.annotation.gtf.gz" % (version,)
                url = "%s%s" % (url,last)
                ftp.cwd(last)
                print "cd ",last
                print "Downloading: %s/%s/%s" % (options.server,url,filename)
                nf = os.path.join(options.output_directory,filename)
                fid = open(nf,'wb')
                ftp.retrbinary("RETR " + filename, fid.write)
                fid.close()

                ftp.close()
                break
            except ftplib.all_errors, e:
                print '\nWarning: FTP attempt %d/%d failed: %s' % (attempt+1, MAX_RETRIES, str(e))
                try:
                    ftp.close()
                except Exception:
                    pass
                if attempt < MAX_RETRIES - 1:
                    time.sleep(60 * (attempt + 1))
                else:
                    print '\nError: FTP failed after %d attempts, aborting' % (MAX_RETRIES,)
                    sys.exit(1)
            except Exception, e:
                print "Error: Generic exception!",str(e)
                sys.exit(1)

        print "Decompressing files ..."
        if filename.endswith('.gz'):
            f_in = gzip.open(nf, 'rb')
            nf_decompressed = nf[:-3]
            f_out = open(nf_decompressed, 'wb')
            shutil.copyfileobj(f_in, f_out)
            f_in.close()
            f_out.close()
            os.remove(nf)
            nf = nf_decompressed

        print "Parsing the GTF file..."
        d = []
        with open(nf, 'r') as _fparse:
            for _line in _fparse:
                if _line.startswith("#") or not _line.strip():
                    continue
                _parts = _line.split("\t")
                if len(_parts) >= 9 and _parts[2] == 'gene':
                    d.append((_parts[0],_parts[3],_parts[4],_parts[6],_parts[8].partition('gene_name "')[2].partition('"')[0]))
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

