#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available sequences of HLA types (human histocompatibility
complex) from the IMGT/HLA Database.

Credits:
* Robinson J, Halliwell JA, McWilliam H, Lopez R, Parham P, Marsh SGE, The IMGT/HLA Database, Nucleic Acids Research (2013) 41:D1222-7
* Robinson J, Malik A, Parham P, Bodmer JG, Marsh SGE: IMGT/HLA - a sequence database for the human major histocompatibility complex Tissue Antigens (2000), 55:280-287



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2020 Daniel Nicorici

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
    description = """It downloads the lastest available sequences of the HLA types (human histocompatibility complex) from the IMGT/HLA Database."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the sequences are downloaded. Default is '%default'.""")

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
                      help="""The EBI server from where the HLA sequences are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    # ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla
    # http://www.ebi.ac.uk/imgt/hla
    # http://www.ebi.ac.uk/imgt/hla/download.html
    # ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta

    # Robinson J, Halliwell JA, McWilliam H, Lopez R, Parham P, Marsh SGE
    # The IMGT/HLA Database
    # Nucleic Acids Research (2013) 41:D1222-7
    # Robinson J, Malik A, Parham P, Bodmer JG, Marsh SGE:
    # IMGT/HLA - a sequence database for the human major histocompatibility complex
    # Tissue Antigens (2000), 55:280-287



    url = '/pub/databases/imgt/mhc/hla/'
    version = ""
    #mypath = os.path.join(options.output_directory,'viruses_index')
    mypath = options.output_directory
    if not os.path.isdir(mypath):
        os.makedirs(mypath)

    final_file = os.path.join(options.output_directory,"hla.fa")
    version = []

    if options.organism.lower() != 'homo_sapiens':
        file(final_file,'w').write('>fake-hla\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n')
        sys.exit(0)

    print "Downloading the HLA sequences from the EBI server!"
    try:
        ftp = ftplib.FTP(options.server)
        print ftp.login()
        ftp.cwd(url)

        list_files = ftp.nlst()

        #list_files = [el for el in list_files if (el.lower().find('nuc') != -1 or el.lower().find('gen') != -1 ) and el.lower().endswith('.fasta')]
        list_files = [el for el in list_files if el.lower() == 'hla_nuc.fasta']

        fid = open(final_file,'wb')
        new_files = []
        for filename in list_files:
            print "Downloading: %s/%s/%s" % (options.server,url,filename)
            ftp.retrbinary("RETR " + filename, fid.write)
            modified_time = ftp.sendcmd('MDTM ' + filename)
            version.append(datetime.datetime.strptime(modified_time[4:], "%Y%m%d%H%M%S").strftime("%Y-%m-%d"))
        fid.close()
        ftp.close()
    except ftplib.all_errors, e:
        print 'FTP Error = ' + str(e)
        sys.exit(1)
    except Exception, e:
        print "Error: Generic exception!",str(e)
        sys.exit(1)

    if version:
        version = sorted(version,reverse = True)
        version = version.pop(0)
        txt =["IMGT/HLA Database: %s\n" % (version,)]
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)
    #
