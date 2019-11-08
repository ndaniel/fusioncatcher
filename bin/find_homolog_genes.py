#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It finds a list of genes that might homologous (there is a short read which maps on both genes).



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2019 Daniel Nicorici

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

"""
Input file looks like:


0001/2	ENST00000274276;ENSG00000145623
0001/2	ENST00000291547;ENSG00000160199
0001/2	ENST00000299163;ENSG00000166135
0001/2	ENST00000343325;ENSG00000173875
0001/2	ENST00000432907;ENSG00000160199
0001/2	ENST00000446165;ENSG00000173875
0001/2	ENST00000482342;ENSG00000177479
0001/2	ENST00000534447;ENSG00000166473
0001/2	ENST00000553458;ENSG00000119711
0001/2	ENST00000587343;ENSG00000267703
0001/2	ENST00000623232;ENSG00000256591
0001/2	ENST09000000012;ENSG09000000012
0002/1	ENST00000274276;ENSG00000145623
0002/1	ENST00000291547;ENSG00000160199
0002/1	ENST00000299163;ENSG00000166135
0002/1	ENST00000343325;ENSG00000173875
0002/1	ENST00000432907;ENSG00000160199
0002/1	ENST00000446165;ENSG00000173875
0002/1	ENST00000482342;ENSG00000177479
0002/1	ENST00000534447;ENSG00000166473
0002/1	ENST00000553458;ENSG00000119711
0002/1	ENST00000587343;ENSG00000267703
0002/1	ENST00000623232;ENSG00000256591
0002/1	ENST09000000012;ENSG09000000012
0003/1	ENST00000299163;ENSG00000166135
0003/1	ENST00000446165;ENSG00000173875
0003/1	ENST00000482342;ENSG00000177479
0003/1	ENST00000553458;ENSG00000119711
0003/1	ENST00000584148;ENSG00000264188
0003/1	ENST00000623232;ENSG00000256591


"""
import sys
import os
import optparse
import gc
import multiprocessing
import itertools

#########################
def line_from3(a_map_filename):
    # it gives chunks from a_map_filename which is assumed to be ordered by the name of transcripts (i.e. column 3)
    fin = None
    if a_map_filename == '-':
        fin = sys.stdin
    elif a_map_filename.lower().endswith('.gz'):
        fin = gzip.open(a_map_filename,'r')
    else:
        fin = open(a_map_filename,'r')

    while True:
        gc.disable()
        lines = fin.readlines(10**8)
        gc.enable()
        if not lines:
            break
        gc.disable()
        lines = [line.rstrip('\r\n').split('\t',2)[:2] for line in lines]
        gc.enable()
        for line in lines:
            if line:
                yield line
    fin.close()

#########################
def line_from4(a_map_filename):
    # it gives chunks from a_map_filename which is assumed to be ordered by the name of transcripts (i.e. column 3)
    fin = None
    if a_map_filename == '-':
        fin = sys.stdin
    elif a_map_filename.lower().endswith('.gz'):
        fin = gzip.open(a_map_filename,'r')
    else:
        fin = open(a_map_filename,'r')

    while True:
        gc.disable()
        lines = fin.readlines(10**8)
        gc.enable()
        if not lines:
            break
        gc.disable()
        lines = [line.rstrip('\r\n').split('\t',3)[:3] for line in lines]
        gc.enable()
        for line in lines:
            if line:
                yield line
    fin.close()

#########################
def read_from3(a_map_filename, database = None, filter_gene = None):
    last_r = ''
    chunk = set()
    last_g = '' # for speed purposes only
    for line in line_from3(a_map_filename):
        if not chunk:
            last_r = line[0]
        if last_r != line[0]: # line[2] is column no 3 in the BOWTIE MAP file which contains the reference sequence name
            if chunk and len(chunk) != 1:
                sc = sorted(chunk)
                if len(sc) > 50: # take only the first 100 genes
                    sc = sc[0:50]
                ex = []
                if database:
                    ex = [database[e] for e in sc]
                yield (last_r,sc,ex)
            last_r = line[0]
            chunk = set()
            last_g = ''
        #tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103
        #g=[el.split('ge=')[1] for el in line[2].split(';') if el.startswith('ge=')]
        #chunk.add(g[0])
#        g1 = line[2].find('ge=') + 3
#        g2 = line[2].find(';',g1+1)
#        g = line[2][g1:g2]
        #g = line[2].partition(';')[2]
        g = line[1]
        if filter_gene and g.startswith(filter_gene):
            continue
        if g != last_g:
            chunk.add(g)
            last_g = g
    if chunk and len(chunk) != 1:
        sc = sorted(chunk)
        ex = []
        if database:
            ex = [database[e] for e in sc]
        yield (last_r,sc,ex)

#########################
def read_from4(a_map_filename, database = None, filter_gene = None):
    last_r = ''
    chunk = set()
    m0 = False
    #m1 = set()
    m2 = set()
    last_g = '' # for speed purposes only
    for line in line_from4(a_map_filename):
        if not chunk:
            last_r = line[0]
        if last_r != line[0]: # line[2] is column no 3 in the BOWTIE MAP file which contains the reference sequence name
            if chunk and len(chunk) != 1:
                if m0 and m2:
                    chunk.difference_update(m2)
                if chunk and len(chunk) != 1:
                    sc = sorted(chunk)
                    if len(sc) > 50: # take only the first 100 genes
                        sc = sc[0:50]
                    ex = []
                    if database:
                        ex = [database[e] for e in sc]
                    yield (last_r,sc,ex)
            last_r = line[0]
            chunk = set()
            m0 = False
            m2 = set()
            last_g = ''
        #tr=ENST00000000233;ge=ENSG00000004059;pn=ENSP00000000233;chr=7;str=+;len=1103
        #g=[el.split('ge=')[1] for el in line[2].split(';') if el.startswith('ge=')]
        #chunk.add(g[0])
#        g1 = line[2].find('ge=') + 3
#        g2 = line[2].find(';',g1+1)
#        g = line[2][g1:g2]
        #g = line[2].partition(';')[2]
        g = line[1]
        if filter_gene and g.startswith(filter_gene):
            continue
        if g != last_g:
            chunk.add(g)
            last_g = g
            if line[2] == "0": # mismatches
                m0 = True
            elif line[2] == "2":
                    m2.add(g)

    if chunk and len(chunk) != 1:
        if m0 and m2:
            chunk.difference_update(m2)
        sc = sorted(chunk)
        ex = []
        if database:
            ex = [database[e] for e in sc]
        yield (last_r,sc,ex)


#########################
def is_overlapping(e1,e2,tolerance = 50):
    r = False
#    e1 = db[g1]
#    e2 = db[g2]
    if e1[2] == e2[2]: # same chromosomes
        if e1[1] - e2[0] > tolerance:
            if e2[1] - e1[0] > tolerance:
                r = True
    return r


#########################
#def compute_homology(ar,ge,ex):
def homology(stuff):
    # ar = read id
    # ge = list of gene ids
    # ex = list of exons coordinates

    ar = stuff[0]
    ge = stuff[1]
    ex = stuff[2]

    n = len(ge)

    flag = False
    hom = []
    # limit the amount of combinations
    if n > 100:
        n = 100
    for a in xrange(0,n-1):
        for b in xrange(a+1,n):
            k = '%s\t%s' % (ge[a],ge[b])
#            if g[a] > g[b]:
#                k = '%s\t%s' % (g[b],g[a])
            hom.append(k)
            if ex and (not is_overlapping(ex[a],ex[b])):
                flag = True
    return (hom, ar if flag else None )

#
# def shred(stuff):
#    return compute_homology(stuff[0],stuff[1],stuff[2])

#
#
#
if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It finds a list of genes that might homologous (there is a short read which maps on both genes)."""
    version="%prog 0.12 beta"

    parser=optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_map_filename",
                      help="""The input file in Bowtie MAP format (sorted by read name) containing the short reads mapped on the transcripts (can be also STDIN).""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output text tab-separated file containing the candidate homologous genes (the genes are sorted alphabetically on the each line).""")

    parser.add_option("--reads",
                      action="store",
                      type="int",
                      dest="reads",
                      default = 1,
                      help="""The minimum number of reads which map simultaneously on two genes in order to be considered as homolog genes. Default is %default.""")


#    parser.add_option("--output_offending_reads",
#                      action="store",
#                      type="string",
#                      dest="output_offending_reads_filename",
#                      help="""The output text file containing the reads names which mapp simultaneously on transcripts from at least two genes.""")

    parser.add_option("--output_offending_pair_reads",
                      action="store",
                      type="string",
                      dest="output_offending_pair_reads_filename",
                      help="""The output text file containing the reads names (and its mate) which mapp simultaneously on transcripts from at least two genes.""")

    parser.add_option("--input_exons",
                      action="store",
                      type="string",
                      dest="exons_filename",
                      help="""Database with exons position on chromosomes, e.g. 'more_exons_ensembl.txt'. This is used for filtering the UTRs extensions by removing any extension which overlaps with any exons from the database. This is optional.""")

    parser.add_option("--filter",
                      action="store",
                      type="string",
                      dest="filter_filename",
                      help="""Input file which contain a pattern for genes which should be ignored/skipped from the analysis.""")

    parser.add_option("--d0",
                      action="store_true",
                      dest="distance_mismatches_0",
                      default = False,
                      help="""If it set then only the alignments of a read are taken into consideration which are at maximum zero mismatches away. This expects that the input has 4 columns instead of 3, and the fourth column contains the mismatches from Bowtie.""")


    parser.add_option("--d1",
                      action="store_true",
                      dest="distance_mismatches_1",
                      default = False,
                      help="""If it set then only the alignments of a read are taken into consideration which are at maximum one mismatch away. This works only for maximum two mismatches. This expects that the input has 4 columns instead of 3, and the fourth column contains the mismatches from Bowtie.""")

    parser.add_option("--output_saved_from_pseudogenes",
                      action="store",
                      type="string",
                      dest="output_saved_from_pseudogenes_filename",
                      help="""A file containing paths to candidate fusion genes and transcripts together with the ids/names of supporting reads.""")

    parser.add_option("--input_pseudogenes",
                      action="store",
                      type="string",
                      dest="input_pseudogenes_filename",
                      help="""The input database with gene ids of the pseudogenes.""")



    parser.add_option("-p", "--processes",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 0,
                      help = """Number of parallel processes/CPUs to be used for computations. In case of value 0 then the program will use all the CPUs which are found. The default value is %default.""")



    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_map_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    #
    cpus = 0
    if options.processes:
        cpus = options.processes
    if cpus == 0:
        cpus = multiprocessing.cpu_count()
    print >>sys.stderr,"Using",cpus,"process(es)..."

    filterout = None
    if options.filter_filename:
        filterout = file(options.filter_filename,"r").readline().rstrip("\r\n")

    pseudogenes=set()
    if options.input_pseudogenes_filename:
        print "Reading...",options.input_pseudogenes_filename
        pseudogenes=set([line.rstrip('\r\n') for line in file(options.input_pseudogenes_filename,'r') if line.rstrip('\r\n')])


    database = {}
    if options.exons_filename:
        #print "Processing the exons database..."
        # ensembl_peptide_id             0
        # ensembl_gene_id                1
        # ensembl_transcript_id          2
        # ensembl_exon_id                3
        # exon_chrom_start               4
        # exon_chrom_end                 5
        # rank                           6
        # start_position                 7
        # end_position                   8
        # transcript_start               9
        # transcript_end                 10
        # strand                         11
        # chromosome_name                12
        # cds_start                      13
        # cds_end                        14
        # 5_utr_start                    15
        # 5_utr_end                      16
        # 3_utr_start                    17
        # 3_utr_end                      18
        exons = [line.rstrip('\r\n').split('\t') for line in file(options.exons_filename,'r').readlines() if line.rstrip('\r\n')]
        exons = [(line[1], # gene_id              0
                  line[7], # gene_start           1
                  line[8], # gene_end             2
                  line[12] # chromosome           3
                  ) for line in exons]

        for line in exons:
            gn = line[0]
            gs = int(line[1])
            ge = int(line[2])
            ch = line[3].upper()

            (gs,ge) = (ge,gs) if gs>ge else (gs,ge)
            if not database.has_key(gn):
                database[gn] = (gs,ge,ch)




    #print "Finding the homolog genes..."
    pool = multiprocessing.Pool(processes=cpus)
    homolog = dict()
    offenders = list()
    fo = None
    if options.output_offending_pair_reads_filename:
        fo = open(options.output_offending_pair_reads_filename,"w")
#    for (a_read,genes,exo) in read_from(options.input_map_filename, database = database, filter_gene = filterout):
#        (h,f) = compute_homology(a_read,genes,exo)
#
#        gc.disable()
#        for k in h:
#            homolog[k] = homolog.get(k,0) + 1
#        if f:
#            offenders.append(a_read)
#        gc.enable()

    my_iter = None
    if options.distance_mismatches_1:
        my_iter = read_from4(
            options.input_map_filename,
            database = database,
            filter_gene = filterout)
    else:
        my_iter = read_from3(
            options.input_map_filename,
            database = database,
            filter_gene = filterout)

    max_genes_per_read = 0
    max_genes_per_read_id = ''

    if fo:
        for w in pool.imap_unordered(
            homology,
            my_iter,
            chunksize = 100):

            h = w[0] # list of pairs of genes which are homologous to each other
            f = w[1] # read id if it the genes are NOT overlapping

            z = len(h)
            if z > max_genes_per_read:
                max_genes_per_read = z
                max_genes_per_read_id = f

            for k in h:
                gc.disable()
                homolog[k] = homolog.get(k,0) + 1
                gc.enable()
            if f:
                offenders.append(f)
                if len(offenders) > 100000:
                    d = list()
                    for e in offenders:
                        if e.endswith('/1'):
                            gc.disable()
                            d.append(e)
                            d.append(e[:-1]+'2')
                            gc.enable()
                        elif e.endswith('/2'):
                            gc.disable()
                            d.append(e[:-1]+'1')
                            d.append(e)
                            gc.enable()
                    fo.writelines([line+'\n' for line in d])
                    offenders = []

        if offenders:
            d = list()
            for e in offenders:
                if e.endswith('/1'):
                    gc.disable()
                    d.append(e)
                    d.append(e[:-1]+'2')
                    gc.enable()
                elif e.endswith('/2'):
                    gc.disable()
                    d.append(e[:-1]+'1')
                    d.append(e)
                    gc.enable()
            fo.writelines([line+'\n' for line in d])
        if fo:
            fo.close()

    else:

        for w in pool.imap_unordered(
            homology,
            my_iter,
            chunksize = 100):


#        for w in itertools.imap(
#            shred,
#            my_iter):



            h = w[0]
            #f = w[1]

            z = len(h)
            if z > max_genes_per_read:
                max_genes_per_read = z
                max_genes_per_read_id = w[1]

            for k in h:
                gc.disable()
                homolog[k] = homolog.get(k,0) + 1
                gc.enable()


#        g = genes
#        n = len(genes)
#        flag = False
#        for a in xrange(0,n-1):
#            for b in xrange(a+1,n):
#                k = '%s\t%s' % (g[a],g[b])
#                if g[a] > g[b]:
#                    k = '%s\t%s' % (g[b],g[a])
#                gc.disable()
#                homolog[k] = homolog.get(k,0) + 1
#                gc.enable()
#                if database and (not is_overlapping(g[a],g[b])):
#                    flag = True
#        if flag:
#            gc.disable()
#            offenders.add(a_read)
#            gc.enable()

    pool.close()
    pool.join()
    #print "Writing...",options.output_filename
    #homolog = sorted([k+'\t'+str(v)+'\n' for (k,v) in homolog.items() if v >= options.reads])
    homolog = [k+'\t'+str(v)+'\n' for (k,v) in homolog.items() if v >= options.reads]
    file(options.output_filename,'w').writelines(homolog)

    #print >> sys.stderr, "Read '%s' found mapping on %d genes!" % (max_genes_per_read_id,max_genes_per_read)

    #if options.output_offending_reads_filename:
        #print "Writing...",options.output_offending_reads_filename
        #gc.disable()
        #d = sorted(set(offenders))
        #gc.enable()
        #file(options.output_offending_reads_filename,'w').writelines([line+'\n' for line in offenders])


    #print "The end."
    #
