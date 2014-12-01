#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It downloads the lastest known fusion genes the Cancer Genome Project
<http://www.sanger.ac.uk/genetics/CGP/Census/>.



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



import os
import sys
import socket
import optparse
import shutil
import urllib2
import symbols
import datetime

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It downloads the lastest known fusion genes from CGP (Cancer Genome Project) database."""
    version = "%prog 0.12 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism","-g",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the known fusion genes are downloaded, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the known fusion genes are stored. Default is '%default'.""")

    parser.add_option("--excel",
                      action="store",
                      type="string",
                      dest="excel_filename",
                      help="""The input Excel file containg the data from the CGP database. Default is '%default'.""")

    parser.add_option("--tsv",
                      action="store",
                      type="string",
                      dest="tsv_filename",
                      help="""The input TSV file containg the data from the CGP database. Default is '%default'.""")

    parser.add_option("--server",
                      action="store",
                      type="string",
                      dest="server",
                      default = "http://cancer.sanger.ac.uk",
                      help="""The CGP server from where the known fusion genes are downloaded. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    try:
        import xlrd
        print "The Python Excel XLRD parser was found!"
    except:
        print >> sys.stderr,"WARNING: Python XLRD library not found!"
        print >> sys.stderr,"Please, install XLRD python library in order to be able to parse the CGP database!"
        file(os.path.join(options.output_directory,'cgp.txt'),'w').write('')
        sys.exit(0)


    # timeout in seconds
    timeout = 1800
    socket.setdefaulttimeout(timeout)
    tmp_file = 'temp_cgp.xls'
    tmp2_file = 'temp_cgp.tsv'

    # http://www.sanger.ac.uk/genetics/CGP/Census/Table_1_full_2012-03-15.xls
    #url = '/genetics/CGP/Census/Table_1_full_2012-03-15.xls'
    #url = "/cancergenome/assets/cancer_gene_census.xls"
    url = "/files/cosmic/current_release/cancer_gene_census.xls"
    url2 = '/cosmic/census/trans?&sEcho=1&iColumns=8&sColumns=&iDisplayStart=0&iDisplayLength=10&sSearch=&bRegex=false&sSearch_0=&bRegex_0=false&bSearchable_0=true&sSearch_1=&bRegex_1=false&bSearchable_1=true&sSearch_2=&bRegex_2=false&bSearchable_2=true&sSearch_3=&bRegex_3=false&bSearchable_3=true&sSearch_4=&bRegex_4=false&bSearchable_4=true&sSearch_5=&bRegex_5=false&bSearchable_5=true&sSearch_6=&bRegex_6=false&bSearchable_6=true&sSearch_7=&bRegex_7=false&bSearchable_7=true&iSortingCols=1&iSortCol_0=0&sSortDir_0=asc&bSortable_0=true&bSortable_1=true&bSortable_2=true&bSortable_3=true&bSortable_4=true&bSortable_5=true&bSortable_6=true&bSortable_7=true&export=tsv'

    if options.organism.lower() == 'homo_sapiens':
        today = datetime.date.today()
        file(os.path.join(options.output_directory,'version.txt'),'a').writelines(['CGP database version: %s\n' % (today.strftime("%Y-%m-%d"),)])

        sem = True
        if options.excel_filename:
            tmp_file = options.excel_filename
        else:
            print "Downloading the known fusion genes from the CGP database!"
            try:
                da1 = urllib2.urlopen('%s%s' % (options.server,url))
                file(tmp_file,'w').write(da1.read())
            except:
                print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty!" % (options.server,url)
                sem = False

        if options.tsv_filename:
            tmp2_file = options.tsv_filename
        else:
            try:
                da2 = urllib2.urlopen('%s%s' % (options.server,url2))
                file(tmp2_file,'w').write(da2.read())
            except:
                print >>sys.stderr, "Warning: Cannot access '%s%s'! The output file will be empty!" % (options.server,url2)
                sem = False


        if sem:

            print "Parsing Excel file..."
            # parse the Excel CGP file with the known fusion genes

            fusions = set()
            d = [line.rstrip('\r\n') for line in file(tmp_file,'r').readlines() if line.rstrip('\r\n')]
            h = d.pop(0) # remove header
            if h.lower().find('<?xml') != -1:
                print >>sys.stderr,"WARNING: Not able to download the CGP database!"
            else:
                wb = xlrd.open_workbook(tmp_file)
                # parse the PubMed sheet
                sheet = wb.sheet_by_name('List') # u'PubMed'
                for row in xrange(sheet.nrows):
                    if row == 0: # skip first line
                        continue
                    line = sheet.row_values(row)
                    g1 = line[0].upper().encode('ascii','ignore')
                    g2_list = line[13].upper().encode('ascii','ignore')
                    g2_list = g2_list.replace(". "," ").replace(' ',',').replace(",,,,",",").replace(",,,",",").replace(",,",",")
                    for el in g2_list.split(','):
                        g2 = el.strip().upper()
                        if g2 and g2 != '?':
                            (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                            fusions.add((g1,g2))

            print "Parsing TSV file..."
            # parse the tsv CGP file with the known fusion genes
            d = [line.rstrip('\r\n').split('\t') for line in file(tmp2_file,'r').readlines() if line.rstrip('\r\n')]
            h = d.pop(0) # remove header
            if h[0].lower().find('<?xml') != -1:
                print >>sys.stderr,"WARNING: Not able to download the CGP database!"
            else:
                for line in d:
                    g1 = line[0].upper().encode('ascii','ignore')
                    g2_list = line[6].upper().encode('ascii','ignore')
                    g2_list = g2_list.replace(". "," ").replace('"',' ').replace(",",";").replace(' ',';').replace(";;;;",";").replace(";;;",";").replace(";;",";")
                    for el in g2_list.split(';'):
                        g2 = el.strip().upper()
                        if g2 and g2 != '?':
                            (g1,g2) = (g2,g1) if g2 < g1 else (g1,g2)
                            fusions.add((g1,g2))

            fusions = list(fusions)
            print "%d known gene fusions found!" % (len(fusions),)

            # read the gene symbols
            file_symbols = os.path.join(options.output_directory,'synonyms.txt')
            genes = symbols.read_genes_symbols(file_symbols)

            banned = set()
            loci = symbols.generate_loci(file_symbols)
            #for v in symbols.locus.values():
            for v in loci.values():
                if v:
                    n = len(v)
                    if n > 1:
                        for i in xrange(n-1):
                            for j in xrange(i+1,n):
                                if v[i].upper() != v[j].upper():
                                    ens1 = symbols.ensembl(v[i].upper(),genes,loci)
                                    ens2 = symbols.ensembl(v[j].upper(),genes,loci)
                                    if ens1 and ens2:
                                        for e1 in ens1:
                                            for e2 in ens2:
                                                if e1 != e2:
                                                    (e1,e2) = (e2,e1) if e2 < e1 else (e1,e2)
                                                    banned.add((e1,e2))


            d = []
            for (g1,g2) in fusions:
                if ( g1.upper() == g2.upper() or ((g1.endswith('@') and g2.endswith('@')) and g1.upper()[:2] == g2.upper()[:2])):
                    print "%s-%s skipped!" % (g1,g2)
                    continue
                ens1 = symbols.ensembl(g1,genes,loci)
                ens2 = symbols.ensembl(g2,genes,loci)

                if ens1 and ens2:
                    for e1 in ens1:
                        for e2 in ens2:
                            if e1 != e2 and ((e1,e2) not in banned) and ((e2,e1) not in banned):
                                d.append([e1,e2])

            data = ['\t'.join(sorted(line)) + '\n' for line in d]
            data = sorted(set(data))

            print "%d known fusion genes converted succesfully to Ensembl Gene ids!" % (len(data),)
        else:
            data = []
        file(os.path.join(options.output_directory,'cgp.txt'),'w').writelines(data)
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)
        if os.path.isfile(tmp2_file):
            os.remove(tmp2_file)
    else:
        # write an empty file for other organisms than human
        file(os.path.join(options.output_directory,'cgp.txt'),'w').write('')
#
