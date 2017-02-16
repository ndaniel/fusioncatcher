#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

It downloads the lastest available gene symbols and synonyms of the gene symbols from Ensembl.



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


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest available gene symbols and synonyms of the gene symbols from Ensembl."""
    version = "%prog 0.11 beta"

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

    parser.add_option("--server-path",
                      action="store",
                      type="string",
                      dest="server_path",
                      default = "/pub/current_mysql/", #ftp://ftp.ensembl.org/pub/release-75/mysql/
                      help="""The path of Ensembl server from where the data is downloaded. Default is '%default'.""")

    parser.add_option("--server",
                      action="store",
                      type="string",
                      dest="server",
                      default = "ftp.ensembl.org",
                      help="""The Ensembl server from where the chromosomes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)

    #ftp://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_71_37/
    #
    #Files:
    #
    # xref.txt.gz
    # gene.txt.gz
    # external_synonym.txt.gz
    # external_db.txt.gz
    # homo_sapiens_core_71_37.sql.gz
    # object_xref.txt.gz
    #


    #ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/
    #ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/
    #ftp://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/dna/

    url = options.server_path

    files = {'xref':'xref.txt.gz',
             'gene':'gene.txt.gz',
             'external_synonym':'external_synonym.txt.gz',
             'sql':'%s_core_---MAGIC---.sql.gz' % (options.organism.lower(),),
             'object_xref':'object_xref.txt.gz'
    }

    print "Downloading the MySQL files of organism '%s' from Ensembl!" % (options.organism.lower(),)
    try:
        ftp = ftplib.FTP(options.server)
        print ftp.login()
        ftp.cwd(url)
        key = "%s_core_" % (options.organism.lower(),)
        organism_dir = [el for el in ftp.nlst() if el.lower().startswith(key) ]
        if organism_dir and len(organism_dir) == 1:
            organism_dir = organism_dir[0]
            magic = organism_dir.lower().split(key)[1]
            files['sql'] = files['sql'].replace("---MAGIC---",magic)
            ftp.cwd(organism_dir)
        else:
            print "ERROR: '%s%s' not found!" % (url,key)
            sys.exit(1)

        for k,filename in files.items():
            print "Downloading: %s%s%s/%s" % (options.server,url,organism_dir,filename)
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
        f = gzip.open(filename, 'rb')
        file_content = f.read()
        f.close()
        fn = filename[:-3]
        files[k] = fn
        fod = open(fn,'wb')
        fod.write(file_content)
        fod.close()
        os.remove(filename)

    print "Parsing the tables in the SQL file..."
    # parse the SQL file for column numbers
    sql = [line.rstrip('\r\n') for line in file(files['sql'],'r').readlines()]
    #sql = [line.strip() for line in file('homo_sapiens_core_71_37.sql','r').readlines() if line.strip()]
    table = dict()
    last_name = None
    idx = -1
    for line in sql:
        if not line:
            continue
        parts = line.split("`")
        n = len(parts)
        if parts[0].find('CREATE TABLE') != -1:
            last_name = parts[1]
            if last_name not in table:
                table[last_name] = {}
                idx = -1
        elif n > 2 and parts[0].find('KEY') == -1 and parts[0].find('ENGINE') == -1 and not parts[0].startswith(")"):
            idx = idx + 1
            col = parts[1]
            table[last_name][col] = idx

    # get the id for ensembl_object_type = 'gene' in table "object_xref"
    object_xref = [line.rstrip('\r\n').split('\t') for line in file(files['object_xref'],'r').readlines() if line.rstrip('\r\n')]
    ox_eot = table['object_xref']['ensembl_object_type']
    ox_xi = table['object_xref']['xref_id']
    #ox_ei = table['object_xref']['ensembl_id']
    object_xref = set([line[ox_xi] for line in object_xref if line[ox_eot].lower() == 'gene' or line[ox_eot].lower() == 'transcript' or line[ox_eot].lower() == 'translation'])


    # get the table "external_synonym"
    external_synonym = [line.rstrip('\r\n').split('\t') for line in file(files['external_synonym'],'r').readlines() if line.rstrip('\r\n')]
    es_xi = table['external_synonym']['xref_id']
    es_s = table['external_synonym']['synonym']
    ids = dict()
    for line in external_synonym:
        id = line[es_xi]
        v = line[es_s].upper().replace("'","").replace('"','')
        if (not v) or (not id) or (id not in object_xref):
            continue
        if not ids.has_key(id):
            ids[id] = set()
        ids[id].add(v)
    external_synonym = None

    # get the table "xref"
    xref = [line.rstrip('\r\n').split('\t') for line in file(files['xref'],'r').readlines() if line.rstrip('\r\n')]
    x_xi = table['xref']['xref_id']
    x_dl = table['xref']['display_label']
    x_da = table['xref']['dbprimary_acc']
    mx = max([x_xi,x_dl,x_da])
    for line in xref:
        if len(line)<mx+1:
            print "Warning:",line
            continue
        id = line[x_xi].strip()
        v = line[x_dl].upper().replace("'","").replace('"','').strip()
        vv = line[x_da].upper().replace("'","").replace('"','').strip()
        if ((not v) and (not vv)) or (not id) or (id not in object_xref):
            continue
        if not ids.has_key(id):
            ids[id] = set()
        if v and not v.isdigit():
            ids[id].add(v)
        if vv and not vv.isdigit():
            ids[id].add(vv)
    xref = None

    # get the table "gene"
    gene = [line.rstrip('\r\n').split('\t') for line in file(files['gene'],'r').readlines() if line.rstrip('\r\n')]
    g_dxi = table['gene']['display_xref_id']
    g_si = table['gene']['stable_id']
    gene_test = [line for line in gene if line[g_si].upper().startswith('ENS')]
    if gene_test:
        gene = gene_test
    else: # this is in case that all ensembl gene ids do not start with ENS
        x = options.organism.upper().split('_')
        templ = "ENS"+x[0][0]+x[1][0:2]+"G"
        for hui in xrange(len(gene)):
            gene[hui][g_si] = templ+gene[hui][g_si].upper()
    g = dict()
    for line in gene:
        id = line[g_dxi].replace("'","").replace('"','')
        v = line[g_si].replace("'","").replace('"','')
        if (not v) or (not id) or (id not in object_xref):
            continue
        if not g.has_key(id):
            g[id] = set()
        g[id].add(v)
    gene = None

    print "Joining all tables..."
    synonyms = []
    for (id,symbols) in ids.iteritems():
        if not g.has_key(id):
            continue
        ens = g[id]
        for e in ens:
            for s in symbols:
                synonyms.append("%s\t%s\n" %(e,s.upper().replace("'","")))
    file(os.path.join(options.output_directory,'synonyms.txt'),'w').writelines(synonyms)
    #
    for f in files.values():
        os.remove(f)
    #
