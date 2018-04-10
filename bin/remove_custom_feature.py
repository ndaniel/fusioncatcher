#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Remove human gene/transcripts/exons from Ensembl database which are wrongly annotated.



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
import symbols



def remove(outdir,
           gene_id = '*',
           transcript_id = '*',
           exon_id = '*',
           log=None):

    # Note:
    # - this initially was made only for transcript_id and therefore transcript_id has the highest priority here
    # - when exon_id is used then also transcript_id needs to be specified

    if log:
        log = file(log,'a')
        log.write("%s\t%s\t%s\n" % (gene_id,transcript_id,exon_id))

    headex = dict([(line.rstrip('\r\n'),i) for i,line in enumerate(file(os.path.join(outdir,'exons_header.txt'),'r').readlines()) if line.rstrip('\r\n')])
    exons = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(outdir,'exons.txt'),'r').readlines() if line.rstrip('\r\n')]
    eid = headex['ensembl_exon_id']
    tid = headex['ensembl_transcript_id']
    gid = headex['ensembl_gene_id']
    p1id = headex['start_position']
    p2id = headex['end_position']
    t1id = headex['transcript_start']
    t2id = headex['transcript_end']


    headge = dict([(line.rstrip('\r\n'),i) for i,line in enumerate(file(os.path.join(outdir,'genes_header.txt'),'r').readlines()) if line.rstrip('\r\n')])
    genes = [line.rstrip('\r\n').split('\t') for line in file(os.path.join(outdir,'genes.txt'),'r').readlines() if line.rstrip('\r\n')]
    gidg = headge['ensembl_gene_id']
    p1idg = headge['start_position']
    p2idg = headge['end_position']
    
    if gene_id != '*':
        # 
        # removing a gene
        #
        
        n = len(exons)
        exons = [line for line in exons if line[gid] != gene_id]
        m = len(exons)
        if n == m:
            print >>sys.stderr, "WARNING: %s gene not found in exons.txt!" % (gene_id,)
        else:
            print " - Found gene %s in 'exons' database!" % (gene_id,)
            file(os.path.join(outdir,'exons.txt'),'w').writelines(['\t'.join(line)+'\n' for line in exons])


        n = len(genes)
        genes = [line for line in genes if line[gidg] != gene_id]
        m = len(genes)
        if n == m:
            print >>sys.stderr, "WARNING: %s gene not found in genes.txt!" % (gene_id,)
        else:
            print " - Found gene %s in 'genes' database!" % (gene_id,)
            file(os.path.join(outdir,'genes.txt'),'w').writelines(['\t'.join(line)+'\n' for line in genes])


        # remove it from the GTF file also
        gtf = [line for line in file(os.path.join(outdir,'organism.gtf'),'r').readlines()]
        n = len(gtf)
        gtf = [line for line in gtf if line.find(gene_id) == -1]
        m = len(gtf)
        if n == m:
            print >>sys.stderr, "WARNING: %s gene not found in organism.gtf!" % (gene_id,)
        else:
            print " - Found gene %s in GTF file!" % (gene_id,)
            file(os.path.join(outdir,'organism.gtf'),'w').writelines(gtf)


    else:
        data = []
        fix = []
        target_gene_id = [line[gid] for line in exons if line[tid] == transcript_id]
        flag = False
        if len(set(target_gene_id)) == 1:

            target_gene_id = target_gene_id.pop()

            print "Processing transcript %s and exon %s (of gene %s) ..." % (transcript_id,target_gene_id,exon_id)

            for line in exons:
                if line[gid] == target_gene_id:
                    if line[tid] != transcript_id:
                        if exon_id == '*':
                            fix.append(line)
                        elif line[eid] != exon_id:
                            fix.append(line)
                else:
                    data.append(line)

            if fix:
                # find new start and end position for gene after the transcript has been removed
                t1 = fix[0][t1id]
                t2 = fix[0][t2id]
                for line in fix:
                    if line[t1id] < t1:
                        t1 = line[t1id]
                    if line[t2id] > t2:
                        t2 = line[t2id]
                for line in fix:
                    line[p1id] = t1
                    line[p2id] = t2
                # join everything
                data.extend(fix)
                # fix the genes.txt

                flag = False
                for line in genes:
                    if line[gidg] == target_gene_id:
                        line[p1idg] = t1
                        line[p2idg] = t2
                        flag = True
                if flag:
                    file(os.path.join(outdir,'genes.txt'),'w').writelines(['\t'.join(line)+'\n' for line in genes])
            else: # fix
                print " - Found gene %s with only one transcript %s. The gene and the transcript will be removed!" % (target_gene_id,transcript_id)
                # a gene with only one transcript which is removed => the gene is removed also
                file(os.path.join(outdir,'genes.txt'),'w').writelines(['\t'.join(line)+'\n' for line in genes if line[gidg] != target_gene_id])
                # remove it from the GTF file also
                gtf = [line for line in file(os.path.join(outdir,'organism.gtf'),'r').readlines() if line.find(target_gene_id) == -1]
                file(os.path.join(outdir,'organism.gtf'),'w').writelines(gtf)

            file(os.path.join(outdir,'exons.txt'),'w').writelines(['\t'.join(line)+'\n' for line in data])

            # remove it from the GTF file also
            if exon_id == '*':
                gtf = [line for line in file(os.path.join(outdir,'organism.gtf'),'r').readlines() if line.find(transcript_id) == -1]
            else:
                gtf = [line for line in file(os.path.join(outdir,'organism.gtf'),'r').readlines() if line.find(transcript_id) == -1 and line.find(exon_id) == -1]
            file(os.path.join(outdir,'organism.gtf'),'w').writelines(gtf)
        else:
            print >>sys.stderr, "WARNING: %s transcript not found!" % (transcript_id,)

    if log:
        log.close()


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Remove human gene/transcripts/exons from Ensembl database which are wrongly annotated."""
    version = "%prog 0.14 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--organism",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = "homo_sapiens",
                      help="""The name of the organism for which the list of allowed candidate fusion genes is generated, e.g. homo_sapiens, mus_musculus, etc. Default is '%default'.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the list of allowed candidate fusion genes is generated. Default is '%default'.""")

    parser.add_option("--transcripts",
                      action="store",
                      type="string",
                      dest="transcripts",
                      help="""List of transcripts which should be removed.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #

    log_file = os.path.join(options.output_directory,'removed_features.txt')
    
    print "Remove genes/transcripts/exons which are wrongly annotated in Ensembl database..."

    if options.transcripts:
        gl = [line.rstrip("\r\n") for line in file(options.transcripts,'r').readlines() if line.rstrip("\r\n")]
        for g in gl:
            remove(options.output_directory,g)

    if options.organism.lower() == 'homo_sapiens':
        remove(options.output_directory,transcript_id="ENST00000467125",log=log_file) # ENST00000467125 is listed as GOPC, but is actually the GOPC-ROS1
        remove(options.output_directory,transcript_id="ENST00000507166",log=log_file) # ENST00000507166 is listed as FIP1L1 but is really the FIP1L1-PDGFRA
#        remove(options.output_directory,transcript_id="ENST00000621209",log=log_file) # ENST00000621209 is listed as CEL but is really the CEL-CELP
        remove(options.output_directory,transcript_id="ENST00000562663",log=log_file) # ENST00000562663 is listed as RGL3, but most likely is EPOR (it overlaps EPOR)
        remove(options.output_directory,transcript_id="ENST00000563726",log=log_file) # ENST00000563726 is listed as RGL3, but most likely is EPOR (it overlaps EPOR)
        remove(options.output_directory,                                gene_id="ENSG00000129965",log=log_file) # remove gene INS-IGF2
        remove(options.output_directory,transcript_id="ENST00000628281",log=log_file) # ENST00000563726 is listed as PIGU, but most likely is PIGU-NCOA6 (it overlaps NCOA6)
        remove(options.output_directory,transcript_id="ENST00000603067",log=log_file) # ENST00000603067 is listed as TAF15, but it overlaps too many genes
#        remove(options.output_directory,transcript_id="ENST00000356424",log=log_file) # ENST00000356424 is listed as SERPINB3, but it overlaps also SERPINB4


        
        
    #
