#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the positions of all exons from the Ensembl database.



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
import sys
import os
import socket
# timeout in seconds
timeout = 6000 # one hour
socket.setdefaulttimeout(timeout)
import urllib
import urllib2
import optparse


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the positions of all exons from the Ensembl database."""
    version = "%prog 0.10 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the exons positions are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the exons positions are stored. Default is '%default'.""")


    parser.add_option("--threshold_length","-t",
                      action = "store",
                      type = "int",
                      dest = "gene_length",
                      default = 150,
                      help = "Genes shorter than this will be skipped and they will not be in the output. "+
                             "Default is '%default'.")

    parser.add_option("--server",
                      action="store",
                      type="string",
                      dest="server",
                      default = "www.ensembl.org",
                      help="""The Ensembl server from where the exons positions are downloaded, e.g. 'www.ensembl.org', 'uswest.ensembl.org', 'useast.ensembl.org', 'asia.ensembl.org', etc. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    #
    #
    #
    ense = options.organism.lower().split('_',1)
    ensembl_organism = ense[0][0]+ense[1]+'_gene_ensembl'


    CHUNK_SIZE=65536 # 2**20 1 MB

    query = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
        <Dataset name = "%%%organism%%%" interface = "default" >
            <Attribute name = "ensembl_peptide_id" />
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "ensembl_exon_id" />
            <Attribute name = "exon_chrom_start" />
            <Attribute name = "exon_chrom_end" />
            <Attribute name = "rank" />
            <Attribute name = "start_position" />
            <Attribute name = "end_position" />
            <Attribute name = "transcript_start" />
            <Attribute name = "transcript_end" />
            <Attribute name = "strand" />
            <Attribute name = "chromosome_name" />
        </Dataset>
    </Query>
    """.replace('%%%organism%%%',ensembl_organism).replace("\n"," ").strip()

    temp_exons_filename = os.path.join(options.output_directory,'temp_exons.txt')
    exons_filename = os.path.join(options.output_directory,'exons.txt')
    exons_header_filename = os.path.join(options.output_directory,'exons_header.txt')

    genes_filename = os.path.join(options.output_directory,'genes.txt')
    genes_header_filename = os.path.join(options.output_directory,'genes_header.txt')

    mydata=urllib.urlencode( {"query" : query} )
    headers = {
    'User-agent': 'Mozilla/5.0 (Windows; U; Windows NT 6.1; en-GB; rv:1.9.2.3) Gecko/20100401 Firefox/3.6.3',
    'Accept' : 'text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*;q=0.5',
    'Accept-Language' : 'en-gb,en;q=0.5'
    }

    print "Starting..."

    s = ""
    ns = 0
    server = "http://%s/biomart/martservice" % (options.server,)
    try:
        req=urllib2.Request(server,mydata,headers)
        page=urllib2.urlopen(req)


        fid=open(temp_exons_filename,'w')
        size=0
        while True:
            part=page.read(CHUNK_SIZE)
            if not part:
                break
            size=size+len(part)
            fid.write(part)
            amount=size/float(1024*1024)
            sys.stdout.write("\b"*ns)
            sys.stdout.flush()
            s = "Downloaded: %9.2f MB" % amount
            ns = len(s)
            sys.stdout.write(s)
            sys.stdout.flush()
        fid.close()
        print "\nDownloaded: %9.2f MB\n" % amount
    except urllib2.HTTPError, error:
        print '\nHTTPError = ' + str(error)
        sys.exit(1)
    except urllib2.URLError, error:
        print '\nURLError = ' + str(error)
        sys.exit(1)
    except IOError, error:
        print '\nIOError = ' + str(error)
        sys.exit(1)
    except Exception, error:
        print "\nError: Generic exception!",str(error)
        sys.exit(1)

    # print keeping only the chromosomes (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT)
    allchr=[str(i) for i in range(1,100)]
    chromosomes=set(allchr+['X','Y','MT','UN'])

    chromosomes=set(el.upper() for el in chromosomes)

    #print "Keeping only the genes/transcript/exons from the chromosomes :",chromosomes

    fid=open(temp_exons_filename,'rt')
    fod=open(exons_filename,'wt')
    removed = []
    while True:
        lines = fid.readlines(10**8)
        if not lines:
            break
        lines = [line.rstrip('\r\n').split("\t") for line in lines if line.rstrip('\r\n')]
        li = []
        for line in lines:
            if (line[12].upper() not in chromosomes) or abs(int(line[7])-int(line[8])) <= options.gene_length:
                removed.append(line)
            else:
                li.append(line)
        lines = li
#        lines = [line for line in lines if line[12].upper() in chromosomes]
#        lines = [line for line in lines if abs(int(line[7])-int(line[8]))>options.gene_length]
        if lines:
            fod.writelines(['\t'.join(line)+'\n' for line in lines])
    fod.close()
    fid.close()

    file(os.path.join(options.output_directory,'exons_removed.txt'),"w").writelines(['\t'.join(line)+'\n' for line in removed])

    os.remove(temp_exons_filename)

    header = [el.split('"')[1]+'\n' for el in query.split('Attribute name =')[1:]]
    file(exons_header_filename,'w').writelines(header)


    # write the gene positions also if I am here
    data = [line.rstrip('\r\n').split('\t') for line in  file(exons_filename,'r').readlines() if line.rstrip('\r\n')]
    data = set([(line[1],line[8],line[7],line[11],line[12]) for line in data])
    data = sorted(data)
    file(genes_filename,'w').writelines(['\t'.join(line)+'\n' for line in data])

    header = (header[1],header[8],header[7],header[11],header[12])
    file(genes_header_filename,'w').writelines(header)
    print "End."
