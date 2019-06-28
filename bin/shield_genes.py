#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Protects a given list of genes against their own pseudogenes and other genes or regions on genome with which share a very high similarity.



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
import sys
import os
import optparse
import gc

def shield(out,gids,readlen=60,step=1,logf=""):
    gids = dict([(k, ('','','','')) for k in gids])
    g = [line.rstrip("\r\n").split("\t") for line in file(os.path.join(out,"genes.txt"),"r") if line.rstrip("\r\n")]
    bed = []
    log = []
    # ensembl_gene_id
    # end_position
    # start_position
    # strand
    # chromosome_name
    for line in g:
        if line[0] in gids:
            log.append("%s\n" % (line[0],))
            if int(line[2]) < int(line[1]):
                gids[line[0]] = (line[4],int(line[2])-1,int(line[1]))
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[2])-1,line[1]))
            else:
                gids[line[0]] = (line[4],int(line[1])-1,int(line[2]))
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[1])-1,line[2]))
    file(os.path.join(out,"shielded_genes.bed"),"w").writelines(bed)
    
    if logf:
        file(os.path.join(out,logf),"w").writelines(log)
    
    bed = []
    bed_max = 10**8
    bedf = open(os.path.join(out,"shielded_reads_similarity.bed"),"w")
    for g in gids.keys():
        t = gids[g]
        c = t[0]
        s = t[1]
        e = t[2]

        k = 0
        gc.disable()
        x = [("%s\t%d\t%d\n" % (c,i,i+readlen)) for i in xrange(s,e,step) if i+readlen < e]

        bedf.writelines(x)
        gc.enable()
            
    bedf.close()
    
    
    
if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Protects a given list of genes against their own pseudogenes and other genes or regions on genome with which share a very high similarity.."""
    version = "%prog 0.17 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism","-g",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the output will be written. Default is '%default'.""")

    parser.add_option("--read-len","-r",
                      action = "store",
                      type = "int",
                      dest = "read_length",
                      default = 60,
                      help = "Read length used to computer the similarity between different regions/genes. "+
                             "Default is '%default'.")

    parser.add_option("--pseudo-genes-check","-x",
                      action = "store_true",
                      default = False,
                      dest = "check_pseudogenes",
                      help = "Skip the pseudogenes check. "+
                             "Default is '%default'.")

    parser.add_option("--use-synonyms","-z",
                      action = "store_true",
                      default = False,
                      dest = "use_synonyms",
                      help = "Use the synonyms symbols for genes. "+
                             "Default is '%default'.")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)

    log_file = os.path.join(options.output_directory,'shielded_genes.txt')
    file(log_file,"w").write("")

    gene_ids = []

    if options.organism.lower() == 'homo_sapiens':
        gene_ids = [
            "ENSG00000181163", # NPM1
            "ENSG00000260596", # DUX4
            "ENSG00000079432", # CIC
            "ENSG00000135486", # HNRNPA1
            "ENSG00000189403", # HMGB1
            "ENSG00000156508", # EEF1A1
            "ENSG00000198830", # HMGN2
            "ENSG00000154582", # TCEB1
            "ENSG00000029993", # HMGB3
            "ENSG00000115541", # HSPE1
            "ENSG00000067900", # ROCK1
            "ENSG00000189403", # HMGB1
            "ENSG00000198804", # MTCO1
            "ENSG00000263001", # GTF2I
            "ENSG00000161960", # EIF4A1
            "ENSG00000196712", # NF1
            "ENSG00000186716", # BCR

        ]
        shield(options.output_directory, gene_ids, options.read_length, logf=log_file)

    #
    # REPORT
    #
    # find candidate genes which might need protection from their corresponding pseudogenes
    #
    symbol = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(options.output_directory,'genes_symbols.txt'),"r").readlines() if line.rstrip("\r\n")] # ensg gene-id
    ensg = dict()
    genename = dict()
    for e,s in symbol:
        if s not in ensg:
            ensg[s] = set()
        ensg[s].add(e)
        genename[e]=s

    synonym = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(options.output_directory,'synonyms.txt'),"r").readlines() if line.rstrip("\r\n")] # ensg gene-id
#    sensg = dict()
#    for e,s in synonym.iteritems():
#        if s not in sensg:
#            sensg[s] = set()
#        sensg[s].add(e)
    
    f = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(options.output_directory,'known.txt'),"r").readlines() if line.rstrip("\r\n")]
    u = []
    for e in f:
        u.append(e[0])
        u.append(e[1])
    known = set(u)
    
    cosmic = set()
    if os.path.exists(os.path.join(options.output_directory,'cosmic.txt')):
        cosmic = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(options.output_directory,'cosmic.txt'),"r").readlines() if line.rstrip("\r\n")]
    u = []
    for e in cosmic:
        u.append(e[0])
        u.append(e[1])
    cosmic = set(u)


    pseudo = set()
    if options.check_pseudogenes:
        pseudo = set([line.rstrip('\r\n') for line in file(os.path.join(options.output_directory,'pseudogenes.txt'),"r").readlines() if line.rstrip("\r\n")])
        print "NOTE: Only pseudogenes are used! ("+str(len(pseudo))+" pseudogenes found)!"
    else:
        pseudo = set([e[0] for e in symbol]) # put here all the genes => simulate that no test is done
        print "NOTE: Pseudogenes info is not used!"


    ribosomal = set([line.rstrip('\r\n') for line in file(os.path.join(options.output_directory,'ribosomal_proteins.txt'),"r").readlines() if line.rstrip("\r\n")])

    #total = symbol + synonym
    total = symbol
    if options.use_synonyms:
        total = symbol + synonym
        print "NOTE: also the synonym symbols for genes are used!"
    r = []
    for w in ('P','L'):
        c = dict()
        en = dict()
        sym = dict()
        for line in total:
            lu1 = line[1].upper() # gene symbol
            if lu1 not in sym:
                sym[lu1] = line[0] # [gene symbol] =  ensembl id
                if lu1.startswith('MT-'):
                    xlu = "MT" + lu1[3:]
                    sym[xlu] = line[0] # [gene symbol] =  ensembl id
            g = lu1.rpartition(w) # gene symbol
            k = g[0] # base of gene symbol (before the P or L)
            if g[0] and g[1] and len(g[0])>2:
                if g[2]:
                    z = True
                    try:
                        y = int(g[2])
                    except:
                        z = False
                    if z and y < 0:
                        z = False
                    if z and (line[0] in pseudo):
                        if k not in c:
                            c[k] = set()
                            en[k] = set()
                        c[k].add(lu1)
                        en[k].add(line[0])
                elif (not g[2]) and (line[0] in pseudo) and g[0] != 'CEL': # add cases like CEL and CELP (but remove the CEL-CELP)
                    if k not in c:
                        c[k] = set()
                        en[k] = set()
                    c[k].add(lu1)
                    en[k].add(line[0])
            
        counts = []
        for k,v in c.iteritems():
            lv = len(v)
            ev = en[k]
            fl = True
            if lv < 2:
                for e in v:
                    if sym[e] not in pseudo:
                        fl = False
                        break
                if k not in sym:
                    fl = False
            if fl:
                counts.append((lv,k))
        counts = sorted(counts)
        counts.reverse()
        
        result = [(line[1] if sym.get(line[1],"") else sorted(c[line[1]])[0], sym.get(line[1],"") if sym.get(line[1],"") else sym.get(sorted(c[line[1]])[0]), line[0],w,','.join(sorted(c[line[1]])),','.join(sorted(en[line[1]]))) for line in counts]

        r.extend(result)
        
    r2 = sorted(r,key = lambda h:(-int(h[2]),h[0]))
    r2 = ["%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % (line[1],line[0],'known_fusion' if line[1] in known else "",'cosmic' if line[1] in cosmic else "",line[2],line[3],line[4],line[5]) for line in r2]
    file(os.path.join(options.output_directory,'shield__report.txt'),"w").writelines(r2)
    #   
    
    # get the pseudogenes of the gene_ids
    gis = set(gene_ids)
    pgs = set()
    pseudo_all = set()
    pseudo_ids = set()
    for line in r:
        if (line[1] in gis) and line[-1]:
            pgs.update(line[-1].split(','))
        if line[-1]:
            pseudo_all.update(line[-1].split(','))
            pseudo_ids.add(line[1])
    pgs.difference_update(gis)


    g = [line.rstrip("\r\n").split("\t") for line in file(os.path.join(options.output_directory,"genes_original.txt"),"r") if line.rstrip("\r\n")]
#    gremoved = [line.rstrip("\r\n").split("\t") for line in file(os.path.join(options.output_directory,"genes_removed.txt"),"r") if line.rstrip("\r\n")]
#    g.extend(gremoved)
    bed = []
    log = []
    # ensembl_gene_id
    # end_position
    # start_position
    # strand
    # chromosome_name

    found = set()
    for line in g:
        if line[0] in pgs:
            found.add(line[0])
            log.append("%s\n" % (line[0],))
            if int(line[2]) < int(line[1]):
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[2])-1,line[1]))
            else:
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[1])-1,line[2]))
    notfound = [e+"\n" for e in sorted(pgs.difference(found))]

    file(os.path.join(options.output_directory,"shield_erase-regions.bed"),"w").write("")
    # add manually the D4Z4 repeat region which contains DUX4 like genes
    if options.organism.lower() == "homo_sapiens":
        r = file(os.path.join(options.output_directory,"genome_information.txt"),'r').readline()
        if r and r.lower().find("grch38") != -1:
            #bed.append("4\t190064000\t190095000\n") # human genome version CRCH38
            eraser = ["4\t190064000\t190095000\n", # for DUX4
                      "10\t133660000\t133774000\n", # dor DUX4
                      "Y\t11303500\t11335000\n", # for DUX4
                      "18\t95000\t130000\n", # for DUX4 (but it covers gene ROCK1P1 which has some parts similar to DUX4)
                      "Y\t1187549\t1212750\n" # for CRLF2
#                      "15\t30355260\t30377930\n" # for gene CHRFAM7A (erase the part of CHRFAM7A which is CHRNA7 in order to detect CHRNA7-FAM7A fusion)
                ]
            file(os.path.join(options.output_directory,"shield_erase-regions.bed"),"w").writelines(eraser) # D4Z4 repeat region

        else:
            print >>sys.stderr,"ERROR: Wrong human genome version found"
            sys.exit(1)


   
    file(os.path.join(options.output_directory,"shield_against_pseudo-genes.bed"),"w").writelines(sorted(bed))
    file(os.path.join(options.output_directory,"shield_against_pseudo-genes.txt"),"w").writelines(sorted(log))
    file(os.path.join(options.output_directory,"shield_against_pseudo-genes_not-found.txt"),"w").writelines(sorted(notfound))

    bed = []
    pseudo_all.difference_update(pseudo_ids)
    ux = []
    for line in g:
        if line[4].lower() == 'mt' or line[0] in ribosomal: # skip the MT genes
            continue
        if line[0] in pseudo_all:
            ux.append(line[0]+"\t"+genename.get(line[0],'')+"\n")
            if int(line[2]) < int(line[1]):
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[2])-1,line[1]))
            else:
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[1])-1,line[2]))
    file(os.path.join(options.output_directory,"shield_pseudogenes-predicted-to-be-erased.bed"),"w").writelines(sorted(bed))
    file(os.path.join(options.output_directory,"shield_pseudogenes-predicted-to-be-erased.txt"),"w").writelines(sorted(ux))



    bed = []
    for line in g:
        if line[0] in pgs:
            if int(line[2]) < int(line[1]):
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[2])-1,line[1]))
            else:
                bed.append("%s\t%d\t%s\n" % (line[4],int(line[1])-1,line[2]))

    file(os.path.join(options.output_directory,"pseudogenes.bed"),"w").writelines(bed)



    #
