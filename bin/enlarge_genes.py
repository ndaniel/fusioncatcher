#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It enlarges the genes slightly in the Ensembl annotation database.


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
import optparse


def give_gene(d, g = 1):
    w = ''
    bucket = []
    for l in d:
        lg = l[g]
        if w != lg:
            if bucket:
                yield bucket
            bucket = [l]
            w = lg
        else:
            bucket.append(l)
    if bucket:
        yield bucket


def int2str(x,n=9):
    x = str(x)
    return '0' * int(n - len(x)) + x



if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It enlarges the genes slightly in the Ensembl annotation database. The enlargement is done only where is possible, for example the genes which overlaps each other will be skipped from enlarging)."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(
        usage = usage,
        description = description,
        version = version)

    parser.add_option("--enlargement-size","-s",
                      action = "store",
                      type = "int",
                      dest = "size",
                      default = 5000,
                      help = "The size of the region, which will be used for "+
                             "enlarging the genes at 5'end and 3'end. One region "+
                             "will be added to the 5' end and a second one to "+
                             "the 3' end. Default is '%default'.")

    parser.add_option("--genes","-e",
                      action = "store",
                      type = "string",
                      dest = "genes",
                      help = "If this is specified then only the genes from this "+
                             "list will be enlarged/covered. If this is not specified then "+
                             "all the genes will be enlarged (only where is possible, for "+
                             "example the genes which overlaps each other will be skipped "+
                             "from enlarging). If the file is empty "+
                             "then no gene will be enlarged!")


    parser.add_option("--gene-length","-l",
                      action = "store",
                      type = "int",
                      dest = "gene_length",
                      default = 1000,
                      help = "Genes having their lengths strictly less than this will have the enlargment done as one continuous region. "+
                             "For the rest of genes two regions will be added at the 5'end and 3'end of the gene. Default is '%default'.")


    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_directory",
                      help = "The output directory where the genes sequences "+
                             "are written. Default is '%default'.")


    (options,args) = parser.parse_args()

    # validate options
    if not (
            options.output_directory
            ):
        parser.print_help()
        parser.error("Missing argument(s)!")
        sys.exit(1)


    #
    #
    #

    size = options.size
    gene_length = options.gene_length

    database_filename = os.path.join(options.output_directory,"exons.txt")

    col = [
           'protein',
           'gene',
           'transcript',
           'exon',
           'exon_start',
           'exon_end',
           'exon_number',
           'gene_start',
           'gene_end',
           'transcript_start',
           'transcript_end',
           'strand',
           'chrom']
    col = dict([ (e,i) for i,e in enumerate(col)])

    database = [line.rstrip('\r\n').split('\t') for line in file(database_filename,'r').readlines() if line.rstrip("\r\n")]

    database = sorted(database, key = lambda x: (x[1],x[2],int(x[6])) )

    cr = col['chrom']
    start = col['gene_start']
    end = col['gene_end']
    gene = col['gene']
    strand = col['strand']
    chrom = sorted(set([line[cr] for line in database]))
    strands = sorted(set([line[strand] for line in database]))

    target = set()
    if options.genes:
        target = set([line.rstrip("\r\n") for line in file(options.genes,"r") if line.rstrip("\r\n")])
        print "Found %d genes to be enlarged!" % (len(target),)

    # take the a gene id and see how it starts
    head = database[0][1]
    m = len(head)
    if head.startswith("ENSG"):
        u = []
        for e in head:
            if e.isdigit():
                break
            else:
                u.append(e)
        head = ''.join(u)
        head_p7 = head[:-1]+"P07"
        head_p9 = head[:-1]+"P09"
        head_g7 = head[:-1]+"G07"
        head_g9 = head[:-1]+"G09"
        head_t7 = head[:-1]+"T07"
        head_t9 = head[:-1]+"T09"
        head_e7 = head[:-1]+"E07"
        head_e9 = head[:-1]+"E09"
        file(os.path.join(options.output_directory,"custom_genes_mark.txt"),"w").write(head_g9)
    else:
        print "Error: unknown Ensembl Id!"
        sys.exit(1)

    if options.genes and (not target):
        print "Exiting gracefully without doing any enlargement! This is ok!"
        sys.exit(0)


    #
    print "Mock-up enlargement of all genes..."
    pads_all = []
    for c in chrom:
        print " * Processing chromosome",c
#        for st in strands: # not a good idea
#            print " * processing strand",st
#            data = set([(line[gene],int(line[start]),int(line[end])) for line in database if line[cr] == c and line[strand] == st])
        if options.genes: # filter out the ENSG09
            data = set([(line[gene],int(line[start]),int(line[end])) for line in database if line[cr] == c and (not line[gene].startswith(head_g9))])
        else: # if not then keep all of them
            data = set([(line[gene],int(line[start]),int(line[end])) for line in database if line[cr] == c])
        data = sorted(data, key = lambda x: (x[1],x[2],x[0]) )
        n = len(data)
        pads = []
        if data:
            # first
            x = data[0][1] - size
            if x < 1:
                x = 1
            pads.append((data[0][0],x,data[0][1]))
            # last
            pads.append((data[-1][0],data[-1][2],data[-1][2]+size))

        for i in xrange(n-1):
            #print " *",data[i][0],data[i][1],data[i][2]
            for j in xrange(i+1,n):
                flag = True
                if data[i][2] > data[j][1]:
                    if data[i][2] < data[j][2]:
                        break
                    else:
                        continue
                for k in xrange(j-1,-1,-1):
                    if i == k:
                        continue
                    if ((data[i][2] < data[k][1] and data[k][1] < data[j][1]) or
                        (data[i][2] < data[k][2] and data[k][2] < data[j][1]) or
                        (data[i][2] > data[k][1] and data[k][2] > data[j][1])):
                        flag = False
                        break
                if flag:
                    # i and j genes are adjacent
                    a1 = data[i][2]
                    a2 = data[i][2] + size
                    b1 = data[j][1] - size
                    b2 = data[j][1]
                    if a2 > b1:
                        if data[i][2] > data[j][1]:
                            print "UPS!",data[i],data[j]
                            sys.exit(1)
                        else:
                            d = (data[i][2] + data[j][1])/2
                            a2 = d - 1
                            b1 = d + 1
                    if a2 > a1 and b1 < b2:
                        #print " *",data[j][0],data[j][1],data[j][2]
                        if a1 > data[j][1] or a2 > data[j][1]:
                            print "Error: Something wrong with",data[i],"new:",a1,a2
                        if b1 < data[i][2] or b2 < data[i][2]:
                            print "Error: Something wrong with",data[j],"new:",b1,b2
                        pads.append((data[i][0],a1,a2))
                        pads.append((data[j][0],b1,b2))
#                            if data[i][0] in ("ENSG00000187266","ENSG00000173928") or data[j][0] in ("ENSG00000187266","ENSG00000173928"):
#                                print "-------------------------------------------------------------------------"
#                                print " * original",data[i][0],data[i][1],data[i][2]
#                                print " * new",data[i][0],a1,a2
#                                print " * original",data[j][0],data[j][1],data[j][2]
#                                print " * new",data[j][0],b1,b2
                    break
        print "  *",len(pads),"enlargements for",len(pads)/2,"gene out of",len(data),"genes"
        if pads:
            pads_all.extend(pads)

    p = dict()
    for line in pads_all:
        if line[0] not in p:
            p[line[0]] = set()
        p[line[0]].add((line[1],line[2]))
    print "TOTAL:",len(p),"mock-up genes enlarged!"

#    print p["ENSG00000000003"]

    w = 0
    data = []
    z = 0
    nn = m - 2 - len(head)
    for ge in give_gene(database, gene):
        #print ge
        #raw_input()
        g = ge[0][gene]
        t = [1 for el in ge if el[0].startswith(head_p7) or el[0].startswith(head_p9)] # do not touch the ones which are already manually modified (which are marked by having a protein which starts with ENSG09 or ENSG07)
        r = False
        if target:
            if g in target:
                r = True
        else:
            r = True

#        if g == "ENSG00000259303" or g == "ENSG00000211976":
#            print "-------->",g,ge,p[g],t,r
        if (g in p) and (not t) and r:
            x = sorted(p[g])
            if len(x) == 2:
                y = [x[0][0], x[0][1], x[1][0], x[1][1],int(ge[0][start]),int(ge[0][end])]
                a = str(min(y))
                b = str(max(y))
            else:
                y = [x[0][0],x[0][1],int(ge[0][start]),int(ge[0][end])]
                a = str(min(y))
                b = str(max(y))
            nge = [(el[0],el[1],el[2],el[3],el[4],el[5],el[6],a,b,el[9],el[10],el[11],el[12]) for el in ge]

            w = w + 1
            if len(nge) == 1 or  abs(int(ge[0][start])-int(ge[0][end])) < gene_length :
                # here the entire gene is covered
                z = z + 1
                h = int2str(z,n=nn)

                nge.append([
                    head_p7+h,
                    g,
                    head_t7+h,
                    head_e7+h,
                    a,
                    b,
                    '1',
                    a,
                    b,
                    a,
                    b,
                    ge[0][strand],
                    ge[0][cr]])
            else:
                # here two paddings are added at the beginning and the end of the gene
                # here no covering of gene is done
                for v in x:
                    z = z + 1
                    h = int2str(z,n=nn)

                    nge.append([
                        head_p7+h,
                        g,
                        head_t7+h,
                        head_e7+h,
                        str(v[0]),
                        str(v[1]),
                        '1',
                        a,
                        b,
                        str(v[0]),
                        str(v[1]),
                        ge[0][strand],
                        ge[0][cr]])

            data.extend(nge)

        elif len(ge) != 1 and (g not in p) and (not t) and r:
            nge = ge[:]

            if abs(int(ge[0][start])-int(ge[0][end])) < gene_length :
                # here the entire gene is covered
                z = z + 1
                h = int2str(z,n=nn)
                w = w + 1

                a = ge[0][start]
                b = ge[0][end]

                nge.append([
                    head_p7+h,
                    g,
                    head_t7+h,
                    head_e7+h,
                    a,
                    b,
                    '1',
                    a,
                    b,
                    a,
                    b,
                    ge[0][strand],
                    ge[0][cr]])
            data.extend(nge)
        else:
            data.extend(ge)

    print "TOTAL:",w,"genes enlarged!"
    print "TOTAL:",z,"additions!"

    file(os.path.join(options.output_directory,"exons.txt"),"w").writelines(['\t'.join(line)+'\n' for line in data])
    x = set(['\t'.join((line[gene],line[end],line[start],line[strand],line[cr]))+'\n' for line in data])
    file(os.path.join(options.output_directory,"genes.txt"),"w").writelines(sorted(x))



    #
