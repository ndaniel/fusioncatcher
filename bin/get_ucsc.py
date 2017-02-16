#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available genes annotations (chromosomal positions and gene symbols) from UCSC.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2017 Daniel Nicorici

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

# filtering the human genes which are marked as overlapping but actually are real fusion genes
skip = set(['ELK4','SLC45A3','ETV6','BCL2L14','NOTCH2NL',
'NBPF10','CCDC142','MRPL53','MAP3K3','DDX42','MATK','ZFR2',
'MBTPS2','YY2','NAIP','OCLN','PPCS','CCDC30','CRLF2'])


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest available genes annotations (chromosomal positions and gene symbols) from UCSC."""
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
                      help="""The UCSC server from where the gene annotations are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz
    #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.sql
    #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
    #ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql
    #
    #Files:
    #
    # kgXref.txt.gz
    # kgXref.sql
    # knownGene.txt.gz
    # knownGene.sql
    #

    org = ''
    if options.organism.lower() == 'homo_sapiens':
        org = 'hg19'
        # find genome information
        d = [line for line in file(os.path.join(options.output_directory,'version.txt'),'r') if line.lower().startswith('genome version') ]
        if d:
            if d[0].lower().find('grch38') !=-1:
                org = 'hg38'
#    elif options.organism.lower() == 'rattus_norvegicus':
#        org = 'rn6'
    elif options.organism.lower() == 'mus_musculus':
        org = 'mm10'
    elif options.organism.lower() == 'canis_familiaris':
        org = 'canFam3'

    if org:
        url = '/goldenPath/%s/database' % (org,)

        files = {'kgXref':'kgXref.txt.gz',
                 'kgXref-sql':'kgXref.sql',
                 'knownGene':'knownGene.txt.gz',
                 'knownGene-sql':'knownGene.sql',
        }

        print "Downloading the SQL files of organism '%s' from UCSC!" % (options.organism.lower(),)
        flag = True
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
            flag = False
        except Exception, e:
            print "Error: Generic exception!",str(e)
            flag = False

        if flag:
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


            # parse "knownGene.txt" file
            gene = [line.rstrip('\r\n').split('\t') for line in file(files['knownGene'],'r').readlines() if line.rstrip('\r\n')]
            gene_id = table['knownGene']['name']
            gene_chrom = table['knownGene']['chrom']
            gene_strand = table['knownGene']['strand']
            gene_start = table['knownGene']['txStart']
            gene_end = table['knownGene']['txEnd']
            gene = dict([(line[gene_id].strip(),
                            (line[gene_chrom],
                             line[gene_strand],
                             line[gene_start],
                             line[gene_end])) for line in gene if line[gene_id]])
            print "%d genes found in 'knownGene.txt'!" % (len(gene),)

            # parse "kgXref.txt" file
            ref = [line.rstrip('\r\n').split('\t') for line in file(files['kgXref'],'r').readlines() if line.rstrip('\r\n')]
            ref_id = table['kgXref']['kgID']
            ref_symbol = table['kgXref']['geneSymbol']
            ref = [(line[ref_id],line[ref_symbol].upper().replace(' ','_')) for line in ref if line[ref_symbol] and line[ref_id]]
            r = dict()
            for k,v in ref:
                x = gene.get(k,None)
                if x:
                    if not r.has_key(v):
                        r[v] = set()
                    r[v].add(x)
            ref = r
            print "%d genes found in 'knownGene.txt' and 'kgXref.txt'!" % (len(ref),)

            # save the UCSC data
            data = []
            for k,v in ref.items():
                for x in v:
                    if x[2] and x[3]:
                        s = int(x[2])
                        e = int(x[3])
                        if s > e:
                            (s,e) = (e,s)
                        data.append([k,str(e),str(s),x[1],x[0]]) # gene_id, end, start, strand, chromosome (to be like ensembl for further processing with generate_overlapping_genes.py)
            # filtering the genes which are marked as overlapping but actually are real fusion genes
            if org in ('hg19','hg38'):
                data = [line for line in data if line[0].upper() not in skip]
            data = ['\t'.join(line)+'\n' for line in data]
            file(os.path.join(options.output_directory,'ucsc_genes.txt'),'w').writelines(data)
            txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
            file(os.path.join(options.output_directory,'ucsc_genes_header.txt'),'w').write(txt)

            txt = ["UCSC database version (%s): %s\n" % (org,version)]
            file(os.path.join(options.output_directory,'version.txt'),'a').writelines(txt)

            #
            for f in files.values():
                os.remove(f)
        else:
            data = []
            file(os.path.join(options.output_directory,'ucsc_genes.txt'),'w').writelines(data)
            txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
            file(os.path.join(options.output_directory,'ucsc_genes_header.txt'),'w').write(txt)

    else:
        data = []
        file(os.path.join(options.output_directory,'ucsc_genes.txt'),'w').writelines(data)
        txt = "%s\n%s\n%s\n%s\n%s\n" % ('gene_symbol','end','start','strand','chromosome')
        file(os.path.join(options.output_directory,'ucsc_genes_header.txt'),'w').write(txt)
    #
