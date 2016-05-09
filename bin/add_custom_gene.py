#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add the custom human genes which are missing from the Ensembl database.



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
import optparse
import symbols

global list_genes

def add(outdir,
        protein_id = '',
        gene_symbol = '',
        gene_id = '',
        transcript_id = '',
        exon_id = '',
        exon_number = '',
        start = '',
        end = '',
        chrom = '',
        strand = '',
        
        start_tr = '',
        end_tr ='',
        start_ex = '',
        end_ex = ''
        ):

    chrom_lens = dict([line.rstrip("\r\n").split("\t") for line in file(os.path.join(outdir,"chromosomes_lengths.txt"),"r") if line.rstrip('\r\n')])
    crlen = int(chrom_lens.get(chrom,'0'))
    if crlen < int(end) or crlen < int(start):
        print >>sys.stderr,"ERROR: The given coordinates of the gene are longer than the known length of the chromosome!"
        sys.exit(1)

    # get genes data
    g = [line.rstrip("\r\n").split("\t") for line in file(os.path.join(outdir,"genes.txt"),"r").readlines() if line.rstrip("\r\n")]
    gid = [(el[0], el[1], el[2], el[3], el[4]) for el in g if el[0] == gene_id]
#    gene id            end_pos     start_pos   strand  chrom 
#    ENSG00000000003	100639991	100627109	-1	X
#    ENSG00000000005	100599885	100584802	1	X
#    ENSG00000000419	50958555	50934867	-1	20
    

    start_tr = start_tr if start_tr else start
    end_tr = end_tr if end_tr else end
    
    start_ex = start_ex if start_ex else start_tr
    end_ex = end_ex if end_ex else end_tr
    #
    #
    #
    if gene_id and (not gid):
        #
        # It is completely new gene!
        #
        print "Gene %s (%s) not found in the database! Gene %s (%s), transcript %s, and exon %s added into database!" % (gene_id,gene_symbol,gene_id,gene_symbol,transcript_id,exon_id)
        file(os.path.join(outdir,'custom_genes.bed'),'a').write('%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,start,end,'%s-%s-%s-%s' % (gene_symbol,gene_id,transcript_id,exon_id),'0','+' if strand == '1' else '-'))
        file(os.path.join(outdir,'descriptions.txt'),'a').write('%s\t\n' % (gene_id,))
        file(os.path.join(outdir,'genes_symbols.txt'),'a').write('%s\t%s\n' % (gene_id,gene_symbol))
        file(os.path.join(outdir,'synonyms.txt'),'a').write('%s\t%s\n' % (gene_id,gene_symbol))
        file(os.path.join(outdir,'genes.txt'),'a').write('%s\t%s\t%s\t%s\t%s\n' % (gene_id,end,start,strand,chrom))
        file(os.path.join(outdir,'exons.txt'),'a').write(
            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                protein_id,
                gene_id,
                transcript_id,
                exon_id,
                start_ex,
                end_ex,
                exon_number,
                start,
                end,
                start_tr,
                end_tr,
                strand,
                chrom))
    elif gid:
        #
        # Gene ID found in the database => it is not new gene
        #
        print "Gene %s (%s) found already in the database!" % (gene_id,gene_symbol)
        if len(gid) == 1:
            gid = gid.pop(0)
        else:
            print "  * Error: Too many genes found in 'genes.txt'!",gid
            sys.exit(1)
        
        #
        # get exons' information for the gene
        #
        e = [line.rstrip("\r\n").split("\t") for line in file(os.path.join(outdir,"exons.txt"),"r").readlines() if line.rstrip("\r\n")]
# ENSP00000362111	ENSG00000000003	ENST00000373020	ENSE00001855382	100636608	100636806	1	100627109	100639991	100628670	100636806	-1	X
# ENSP00000362111	ENSG00000000003	ENST00000373020	ENSE00003662440	100635558	100635746	2	100627109	100639991	100628670	100636806	-1	X
# ENSP00000362111	ENSG00000000003	ENST00000373020	ENSE00003654571	100635178	100635252	3	100627109	100639991	100628670	100636806	-1	X
# ensembl_peptide_id
#                   ensembl_gene_id
#                                   ensembl_transcript_id
#                                                   ensembl_exon_id
#                                                                   exon_chrom_start
#                                                                               exon_chrom_end
#                                                                                           rank
#                                                                                               start_position
#                                                                                                           end_position
#                                                                                                                         transcript_start
#                                                                                                                                       transcript_end
#                                                                                                                                                 strand
#                                                                                                                                                      chromosome_name

        et = [line for line in e if line[1] == gene_id and line[2] == transcript_id]
        ete = [line for line in et if line[3] == exon_id]
        
        if ete:
            print "  * Gene %s (%s), transcript %s, exon %s are already in the database and they will NOT be added again!" % (gene_id,gene_symbol,transcript_id,exon_id)
        elif et: 
            print "  * New exon %s added for already present gene %s (%s) and transcript %s." % (exon_id,gene_id,gene_symbol,transcript_id)
            file(os.path.join(outdir,'custom_genes.bed'),'a').write('chr%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,start_ex,end_ex,'%s-%s-%s-%s-%s' % (gene_symbol,gene_id,transcript_id,exon_id,exon_number),'0','+' if strand == '1' else '-'))
            file(os.path.join(outdir,'exons.txt'),'a').write(
                '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                    protein_id,
                    gene_id,
                    transcript_id,
                    exon_id,
                    start_ex,
                    end_ex,
                    exon_number,
                    start,
                    end,
                    start_tr,
                    end_tr,
                    strand,
                    chrom))
        else: # exon and transcript not found in database
            # new transcript found (which is not in database)
            if int(end) > int(gid[1]) or int(start) < int(gid[2]) or strand != gid[3] or chrom != gid[4]: # testing for new transcript and new gene coordinates
                # new transcript found with different genes coordinates => just add it
                print "  * Gene %s (%s) requires changes in the entire database! New transcript %s and new exon %s added for this gene!" % (gene_id,gene_symbol,transcript_id,exon_id)
                # update the genes.txt
                g = [line for line in g if line[0] != gene_id]
                g.append([gene_id,end,start,strand,chrom])
                file(os.path.join(outdir,'genes.txt'),'w').writelines(['\t'.join(line)+'\n' for line in g])
                # update the exons.txt
                e_rest = [line for line in e if line[1] != gene_id]
                e_target = [line for line in e if line[1] == gene_id]
                e_target = [[line[0],line[1],line[2],line[3],line[4],line[5],line[6],start,end,line[9],line[10],strand,chrom] for line in e_target]
                e_target.append([
                        protein_id,
                        gene_id,
                        transcript_id,
                        exon_id,
                        start_ex,
                        end_ex,
                        exon_number,
                        start,
                        end,
                        start_tr,
                        end_tr,
                        strand,
                        chrom] )
                e = e_rest + e_target
                file(os.path.join(outdir,'exons.txt'),'w').writelines(['\t'.join(line)+'\n' for line in e])
                file(os.path.join(outdir,'custom_genes.bed'),'a').write('chr%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,start,end,'%s-%s-%s-%s' % (gene_symbol,gene_id,transcript_id,exon_id),'0','+' if strand == '1' else '-'))
                file(os.path.join(outdir,'custom_genes.bed'),'a').write('chr%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,start_ex,end_ex,'%s-%s-%s-%s-%s' % (gene_symbol,gene_id,transcript_id,exon_id,exon_number),'0','+' if strand == '1' else '-'))
            else:
                print "  * Gene %s (%s) does NOT require changes in the entire database! New transcript %s and new exon %s added for this gene!" % (gene_id,gene_symbol,transcript_id,exon_id)
                file(os.path.join(outdir,'exons.txt'),'a').write(
                    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        protein_id,
                        gene_id,
                        transcript_id,
                        exon_id,
                        start_ex,
                        end_ex,
                        exon_number,
                        start,
                        end,
                        start_tr,
                        end_tr,
                        strand,
                        chrom))
                file(os.path.join(outdir,'custom_genes.bed'),'a').write('chr%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,start_ex,end_ex,'%s-%s-%s-%s-%s' % (gene_symbol,gene_id,transcript_id,exon_id,exon_number),'0','+' if strand == '1' else '-'))
    print ""
    list_genes.append(gene_id)

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Add the custom human genes which are missing from the Ensembl database."""
    version = "%prog 0.16 beta"

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

    (options,args) = parser.parse_args()

    # validate options
    if not (options.output_directory
            ):
        parser.print_help()
        sys.exit(1)


    #
    #
    #
    list_genes = []

    print "Add/change/correct the genes which have errors or are missing from the Ensembl database..."

    file(os.path.join(options.output_directory,"custom_genes.txt"),"w").write('')
    file(os.path.join(options.output_directory,"custom_genes.bed"),"w").write('')
    file(os.path.join(options.output_directory,"custom_genes_mark.txt"),"w").write('')
    file(os.path.join(options.output_directory,"custom_transcripts_mark.txt"),"w").write('')
    
    database_filename = os.path.join(options.output_directory,"exons.txt")
    database = file(database_filename,'r').readline().rstrip('\r\n').split('\t')
    # take the a gene id and see how it starts
    head = database[1]
    m = len(head)
    if head.startswith("ENS"):
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
        file(os.path.join(options.output_directory,"custom_transcripts_mark.txt"),"w").write(head_t9+'\n'+head_t7)
    else:
        print "Error: unknown Ensembl Id!",head
        sys.exit(1)
        
        
        


    if options.organism.lower() == 'mus_musculus':
        pass
    elif options.organism.lower() == 'rattus_norvegicus':
        pass
    elif options.organism.lower() == 'canis_familiaris':
        pass
    elif options.organism.lower() == 'homo_sapiens':

        # find genome information
        d = [line for line in file(os.path.join(options.output_directory,'version.txt'),'r') if line.lower().startswith('genome version') ]
        if d:
            if d[0].lower().find('grch37') !=-1:
                print "Found version GRCh37 human genome version!"
                # coordinates valid only for GRCh37
                add(outdir = options.output_directory,
                    protein_id = '',
                    gene_symbol = 'C19MC',
                    gene_id = 'ENSG09000000001',
                    transcript_id = 'ENST09000000001',
                    exon_id = 'ENSE09000000001',
                    exon_number = '1',
                    start = '54146160',
                    end = '54280800',
                    chrom = '19',
                    strand = '1'
                )

                # coordinates valid only for GRCh37
                add(outdir = options.output_directory,
                    protein_id = '',
                    gene_symbol = 'MIR-371-CLUSTER',
                    gene_id = 'ENSG09000000002',
                    transcript_id = 'ENST09000000002',
                    exon_id = 'ENSE09000000002',
                    exon_number = '1',
                    start = '54281000',
                    end = '54295770',
                    chrom = '19',
                    strand = '1'
                )

#                coordinates valid only for GRCh37
#                add(outdir = options.output_directory,
#                    protein_id = '',
#                    gene_symbol = 'AL035685.1',
#                    gene_id = 'ENSG00000236127',
#                    transcript_id = 'ENST09000000003',
#                    exon_id = 'ENSE09000000003',
#                    exon_number = '1',
#                    start = '47933000',
#                    end = '47945900',
#                    chrom = '20',
#                    strand = '1'
#                )

                # coordinates valid only for GRCh37
                add(outdir = options.output_directory,
                    protein_id = '',
                    gene_symbol = 'DA750114',
                    gene_id = 'ENSG09000000004',
                    transcript_id = 'ENST09000000004',
                    exon_id = 'ENSE09000000004',
                    exon_number = '1',
                    start = '108541000',
                    end = '109040900',
                    chrom = '9',
                    strand = '1'
                )

                # coordinates valid only for GRCh37
                add(outdir = options.output_directory,
                    protein_id = '',
                    gene_symbol = 'AC008746.10',
                    gene_id = 'ENSG00000237955',
                    transcript_id = 'ENST09000000005',
                    exon_id = 'ENSE09000000005',
                    exon_number = '1',
                    start = '54883000', #54890500
                    end = '54926000', #54891700
                    chrom = '19',
                    strand = '1'
                )


            ####################################################################
            # human GRCh38/hg38
            ####################################################################
            elif d[0].lower().find('grch38') !=-1:
                print "Found version GRCh38 human genome version!"
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000001',
                    gene_symbol = 'C19MC',
                    gene_id = 'ENSG09000000001',
                    transcript_id = 'ENST09000000001',
                    exon_id = 'ENSE09000000001',
                    exon_number = '1',
                    start = '53641443',
                    end = '53780750',
                    chrom = '19',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000002',
                    gene_symbol = 'MIR-371-CLUSTER',
                    gene_id = 'ENSG09000000002',
                    transcript_id = 'ENST09000000002',
                    exon_id = 'ENSE09000000002',
                    exon_number = '1',
                    start = '53782000',
                    end = '53792600',
                    chrom = '19',
                    strand = '1'
                )

#                add(outdir = options.output_directory,
#                    protein_id = '',
#                    gene_symbol = 'AL035685.1',
#                    gene_id = 'ENSG00000236127',
#                    transcript_id = 'ENST09000000003',
#                    exon_id = 'ENSE09000000003',
#                    exon_number = '1',
#                    start = '49310000',
#                    end = '49329500',
#                    chrom = '20',
#                    strand = '1'
#                )

                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000004',
                    gene_symbol = 'DA750114',
                    gene_id = 'ENSG09000000004',
                    transcript_id = 'ENST09000000004',
                    exon_id = 'ENSE09000000004',
                    exon_number = '1',
                    start = '105817000',
                    end = '106278600',
                    chrom = '9',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000005',
                    gene_symbol = 'AC008746.10', # stjude
                    gene_id = 'ENSG00000237955',
                    transcript_id = 'ENST09000000005',
                    exon_id = 'ENSE09000000005',
                    exon_number = '1',
                    start = '54371500', # 54378500
                    end = '54414500', # 54380200
                    chrom = '19',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000006',
                    gene_symbol = 'CRLF2', # stjude
                    gene_id = 'ENSG00000205755',
                    transcript_id = 'ENST09000000006',
                    exon_id = 'ENSE09000000006',
                    exon_number = '1',
                    start = '1115000', # 54378500
                    end = '1250000', #'1220000', # 1267000 #1393000
                    chrom = 'X',
                    strand = '-1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000300006',
                    gene_symbol = 'CRLF2', # stjude
                    gene_id = 'ENSG00000205755',
                    transcript_id = 'ENST09000300006',
                    exon_id = 'ENSE09000300006',
                    exon_number = '1',
                    start = '1232040', 
                    end = '1232227', 
                    chrom = 'X',
                    strand = '-1'
                )

                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000007',
                    gene_symbol = 'CSF2RA', # stjude
                    gene_id = 'ENSG00000198223',
                    transcript_id = 'ENST09000000007',
                    exon_id = 'ENSE09000000007',
                    exon_number = '1',
                    start = '1250001', #'1220001', #'1213300'
                    end = '1322000', # 54380200
                    chrom = 'X',
                    strand = '1'
                )
                

                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000008',
                    gene_symbol = 'IL3RA', # stjude
                    gene_id = 'ENSG00000185291',
                    transcript_id = 'ENST09000000008',
                    exon_id = 'ENSE09000000008',
                    exon_number = '1',
                    start = '1322001', #'1213300'
                    end = '1383500', # 54380200
                    chrom = 'X',
                    strand = '1'
                )

                #
                # IGK locus -- split in several pieces
                #

                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000009',
                    gene_symbol = 'IGK_locus_(a)', # stjude
                    gene_id = 'ENSG09000000009',
                    transcript_id = 'ENST09000000009',
                    exon_id = 'ENSE09000000009',
                    exon_number = '1',
                    start = '88846000', #'1213300'
                    end =   '89154500', # 54380200
                    chrom = '2',
                    strand = '-1'
                )

                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000010',
                    gene_symbol = 'IGK_locus_(b)', # stjude
                    gene_id = 'ENSG09000000010',
                    transcript_id = 'ENSG09000000010',
                    exon_id = 'ENSE09000000010',
                    exon_number = '1',
                    start = '89154501', #'1213300'
                    end =   '89463000', # 54380200
                    chrom = '2',
                    strand = '-1'
                )
                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000011',
                    gene_symbol = 'IGK_locus_(c)', # stjude
                    gene_id = 'ENSG09000000011',
                    transcript_id = 'ENST09000000011',
                    exon_id = 'ENSE09000000011',
                    exon_number = '1',
                    start = '89521000', #'1213300'
                    end =   '89951250', # 54380200
                    chrom = '2',
                    strand = '-1'
                )
                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000012',
                    gene_symbol = 'IGK_locus_(d)', # stjude
                    gene_id = 'ENSG09000000012',
                    transcript_id = 'ENST09000000012',
                    exon_id = 'ENSE09000000012',
                    exon_number = '1',
                    start = '89951251', #'1213300'
                    end =   '90381500', # 54380200
                    chrom = '2',
                    strand = '-1'
                )
                
                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000009',
                    gene_symbol = 'IGK_locus_(e)', # stjude
                    gene_id = 'ENSG09000000009',
                    transcript_id = 'ENST09000000009',
                    exon_id = 'ENSE09000000009',
                    exon_number = '1',
                    start = '88846000', #'1213300'
                    end =   '89154500', # 54380200
                    chrom = '2',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001010',
                    gene_symbol = 'IGK_locus_(f)', # stjude
                    gene_id = 'ENSG09000001010',
                    transcript_id = 'ENSG09000001010',
                    exon_id = 'ENSE09000001010',
                    exon_number = '1',
                    start = '89154501', #'1213300'
                    end =   '89463000', # 54380200
                    chrom = '2',
                    strand = '1'
                )
                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001011',
                    gene_symbol = 'IGK_locus_(g)', # stjude
                    gene_id = 'ENSG09000001011',
                    transcript_id = 'ENST09000001011',
                    exon_id = 'ENSE09000001011',
                    exon_number = '1',
                    start = '89521000', #'1213300'
                    end =   '89951250', # 54380200
                    chrom = '2',
                    strand = '1'
                )
                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001012',
                    gene_symbol = 'IGK_locus_(h)', # stjude
                    gene_id = 'ENSG09000001012',
                    transcript_id = 'ENST09000001012',
                    exon_id = 'ENSE09000001012',
                    exon_number = '1',
                    start = '89951251', #'1213300'
                    end =   '90381500', # 54380200
                    chrom = '2',
                    strand = '1'
                )


                #
                # IGH locus -- split in several pieces
                #
                # IGH_locus: 14::+:chr14:105,556,000-106,883,700

                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000013',
                    gene_symbol = 'IGH_locus_(a)', # stjude
                    gene_id = 'ENSG09000000013',
                    transcript_id = 'ENST09000000013',
                    exon_id = 'ENSE09000000013',
                    exon_number = '1',
                    start = '105556000', #
                    end =   '105778000', #
                    chrom = '14',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000014',
                    gene_symbol = 'IGH_locus_(b)', # stjude
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000000014',
                    exon_id = 'ENSE09000000014',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000015',
                    gene_symbol = 'IGH_locus_(c)', # stjude
                    gene_id = 'ENSG09000000015',
                    transcript_id = 'ENST09000000015',
                    exon_id = 'ENSE09000000015',
                    exon_number = '1',
                    start = '106000001', #
                    end =   '106221250', #
                    chrom = '14',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000016',
                    gene_symbol = 'IGH_locus_(d)', # stjude
                    gene_id = 'ENSG09000000016',
                    transcript_id = 'ENST09000000016',
                    exon_id = 'ENSE09000000016',
                    exon_number = '1',
                    start = '106221251', #
                    end =   '106442500', #
                    chrom = '14',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000017',
                    gene_symbol = 'IGH_locus_(e)', # stjude
                    gene_id = 'ENSG09000000017',
                    transcript_id = 'ENST09000000017',
                    exon_id = 'ENSE09000000017',
                    exon_number = '1',
                    start = '106442501', #
                    end =   '106663100', #
                    chrom = '14',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000018',
                    gene_symbol = 'IGH_locus_(f)', # stjude
                    gene_id = 'ENSG09000000018',
                    transcript_id = 'ENST09000000018',
                    exon_id = 'ENSE09000000018',
                    exon_number = '1',
                    start = '106663101', #
                    end =   '106883700', #
                    chrom = '14',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001013',
                    gene_symbol = 'IGH_locus_(g)', # stjude
                    gene_id = 'ENSG09000001013',
                    transcript_id = 'ENST09000001013',
                    exon_id = 'ENSE09000001013',
                    exon_number = '1',
                    start = '105556000', #
                    end =   '105778000', #
                    chrom = '14',
                    strand = '-1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001014',
                    gene_symbol = 'IGH_locus_(h)', # stjude
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001014',
                    exon_id = 'ENSE09000001014',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001015',
                    gene_symbol = 'IGH_locus_(i)', # stjude
                    gene_id = 'ENSG09000001015',
                    transcript_id = 'ENST09000001015',
                    exon_id = 'ENSE09000001015',
                    exon_number = '1',
                    start = '106000001', #
                    end =   '106221250', #
                    chrom = '14',
                    strand = '-1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001016',
                    gene_symbol = 'IGH_locus_(j)', # stjude
                    gene_id = 'ENSG09000001016',
                    transcript_id = 'ENST09000001016',
                    exon_id = 'ENSE09000001016',
                    exon_number = '1',
                    start = '106221251', #
                    end =   '106442500', #
                    chrom = '14',
                    strand = '-1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001017',
                    gene_symbol = 'IGH_locus_(k)', # stjude
                    gene_id = 'ENSG09000001017',
                    transcript_id = 'ENST09000001017',
                    exon_id = 'ENSE09000001017',
                    exon_number = '1',
                    start = '106442501', #
                    end =   '106663100', #
                    chrom = '14',
                    strand = '-1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001018',
                    gene_symbol = 'IGH_locus_(l)', # stjude
                    gene_id = 'ENSG09000001018',
                    transcript_id = 'ENST09000001018',
                    exon_id = 'ENSE09000001018',
                    exon_number = '1',
                    start = '106663101', #
                    end =   '106883700', #
                    chrom = '14',
                    strand = '-1'
                )


                #
                #
                #

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000020',
                    gene_symbol = 'CRHR1-IT1_', # stjude # overlaps CRHR1-IT1 on opposite strand
                    gene_id = 'ENSG09000000020',
                    transcript_id = 'ENST09000000020',
                    exon_id = 'ENSE09000000020',
                    exon_number = '1',
                    start = '45614000', #
                    end =   '45651000', #
                    chrom = '17',
                    strand = '-1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000021',
                    gene_symbol = 'CRHR1-IT1', # stjude # overlaps CRHR1-IT1 on forward strand
                    gene_id = 'ENSG00000204650',
                    transcript_id = 'ENST09000000021',
                    exon_id = 'ENSE09000000021',
                    exon_number = '1',
                    start = '45614000', #
                    end =   '45651000', #
                    chrom = '17',
                    strand = '1'
                )


                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000030',
                    gene_symbol = 'WHSC1', #
                    gene_id = 'ENSG00000109685',
                    transcript_id = 'ENST09000000030',
                    exon_id = 'ENSE09000000030',
                    exon_number = '1',
                    start = '1865000', #1,871,424-1,982,207
                    end =   '1982500', #
                    chrom = '4',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000031',
                    gene_symbol = 'TMEM38B', #
                    gene_id = 'ENSG00000095209',
                    transcript_id = 'ENST09000000031',
                    exon_id = 'ENSE09000000031',
                    exon_number = '1',
                    start = '105670000', #chr9:105,670,062-105,816,113
                    end =   '105816100', #
                    chrom = '9',
                    strand = '1'
                )

                # big one
                # transcript 1 -- exon 1/1 - EPOR (overlaps entirely EPOR)
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000033',
                    gene_symbol = 'EPOR', #
                    gene_id = 'ENSG00000187266',
                    transcript_id = 'ENST09000000033',
                    exon_id = 'ENSE09000000033',
                    exon_number = '1',
                    start = '11377000', #
                    end =   '11394000', #
                    chrom = '19',
                    strand = '-1'
                )

                # transcript 2 - exon 1/2-- EPOR
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000035',
                    gene_symbol = 'EPOR', # stjude # overlaps EPOR on opposite strand
                    gene_id = 'ENSG00000187266',
                    transcript_id = 'ENST09000000036',
                    exon_id = 'ENSE09000000036',
                    exon_number = '1',
                    start = '11377000', #
                    end =   '11394000', #
                    chrom = '19',
                    strand = '-1',
                    
                    start_tr = '11384093',
                    end_tr = '11388928',
                    start_ex = '11388848', #
                    end_ex =   '11388928' #
                )

                # transcript 2 - exon 2/2-- EPOR
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000035',
                    gene_symbol = 'EPOR', # stjude 
                    gene_id = 'ENSG00000187266',
                    transcript_id = 'ENST09000000036',
                    exon_id = 'ENSE09000000037',
                    exon_number = '2',
                    start = '11377000', #
                    end =   '11394000', #
                    chrom = '19',
                    strand = '-1',
                    

                    start_tr = '11384093',
                    end_tr = '11388928',
                    start_ex = '11384093',
                    end_ex = '11384644'
                )

                # transcript 3 - exon 1/3-- EPOR
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000038',
                    gene_symbol = 'EPOR', # stjude 
                    gene_id = 'ENSG00000187266',
                    transcript_id = 'ENST09000000038',
                    exon_id = 'ENSE09000000037',
                    exon_number = '1',
                    start = '11377000', #
                    end =   '11394000', #
                    chrom = '19',
                    strand = '-1',

                    start_tr = '11384093',
                    end_tr = '11388928',
                    start_ex = '11388848', #
                    end_ex =   '11388928' #
                )

                # transcript 3 - exon 2/3 -- EPOR
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000038',
                    gene_symbol = 'EPOR', # stjude 
                    gene_id = 'ENSG00000187266',
                    transcript_id = 'ENST09000000038',
                    exon_id = 'ENSE09000000038',
                    exon_number = '2',
                    start = '11377000', #
                    end =   '11394000', #
                    chrom = '19',
                    strand = '-1',
                    
                    start_tr = '11384093',
                    end_tr = '11388928',
                    start_ex = '11384506', #
                    end_ex =   '11384644' #
                )

                # transcript 3 - exon 3/3 -- EPOR
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000038',
                    gene_symbol = 'EPOR', # stjude 
                    gene_id = 'ENSG00000187266',
                    transcript_id = 'ENST09000000038',
                    exon_id = 'ENSE09000000039',
                    exon_number = '3',
                    start = '11377000', #
                    end =   '11394000', #
                    chrom = '19',
                    strand = '-1',

                    start_tr = '11384093',
                    end_tr = '11388928',
                    start_ex = '11384093', #
                    end_ex =   '11384507' #
                )


                # transcript 1 - exon 1/1 -- RP11-167H9.4
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000050',
                    gene_symbol = 'RP11-167H9.4', # 
                    gene_id = 'ENSG00000243944',
                    transcript_id = 'ENST09000000050',
                    exon_id = 'ENSE09000000050',
                    exon_number = '1',
                    start = '150039214', #
                    end =   '150213726', #
                    chrom = '3',
                    strand = '1',
                )


                # transcript 1 - exon 1/1 -- RP11-167H9.4
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000055',
                    gene_symbol = 'OLFM1', # 
                    gene_id = 'ENSG00000130558',
                    transcript_id = 'ENST09000000055',
                    exon_id = 'ENSE09000000055',
                    exon_number = '1',
                    start = '135075243', #
                    end =   '135121179', #
                    chrom = '9',
                    strand = '1',
                )

                # transcript 1 - exon 1/1 -- EWSR1
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000060',
                    gene_symbol = 'EWSR1', # 
                    gene_id = 'ENSG00000182944',
                    transcript_id = 'ENST09000000060',
                    exon_id = 'ENSE09000000060',
                    exon_number = '1',
                    start = '29268009', #
                    end =   '29300525', #
                    chrom = '22',
                    strand = '1',
                )

                # TRA locus
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001070',
                    gene_symbol = 'TRA_locus_(a)', 
                    gene_id = 'ENSG09000001070',
                    transcript_id = 'ENST09000001070',
                    exon_id = 'ENSE09000001070',
                    exon_number = '1',
                    start = '21549840', #
                    end =   '22056909', #
                    chrom = '14',
                    strand = '1'
                )
                # TRA locus
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001071',
                    gene_symbol = 'TRA_locus_(b)', 
                    gene_id = 'ENSG09000001071',
                    transcript_id = 'ENST09000001071',
                    exon_id = 'ENSE09000001071',
                    exon_number = '1',
                    start = '22056910', #
                    end =   '22563979', #
                    chrom = '14',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001072',
                    gene_symbol = 'TRA_locus_(c)', 
                    gene_id = 'ENSG09000001072',
                    transcript_id = 'ENST09000001072',
                    exon_id = 'ENSE09000001072',
                    exon_number = '1',
                    start = '21549840', #
                    end =   '22056909', #
                    chrom = '14',
                    strand = '-1'
                )
                # TRA locus
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001073',
                    gene_symbol = 'TRA_locus_(d)', 
                    gene_id = 'ENSG09000001073',
                    transcript_id = 'ENST09000001073',
                    exon_id = 'ENSE09000001073',
                    exon_number = '1',
                    start = '22056910', #
                    end =   '22563979', #
                    chrom = '14',
                    strand = '-1'
                )


                # TRB locus
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001075',
                    gene_symbol = 'TRB_locus_(a)', 
                    gene_id = 'ENSG09000001075',
                    transcript_id = 'ENST09000001075',
                    exon_id = 'ENSE09000001075',
                    exon_number = '1',
                    start = '142274400', #
                    end =   '142821000', #
                    chrom = '7',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001076',
                    gene_symbol = 'TRB_locus_(b)', 
                    gene_id = 'ENSG09000001076',
                    transcript_id = 'ENST09000001076',
                    exon_id = 'ENSE09000001076',
                    exon_number = '1',
                    start = '142274400', #
                    end =   '142821000', #
                    chrom = '7',
                    strand = '-1'
                )


                # TRG locus
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001077',
                    gene_symbol = 'TRG_locus_(a)', 
                    gene_id = 'ENSG09000001077',
                    transcript_id = 'ENST09000001077',
                    exon_id = 'ENSE09000001077',
                    exon_number = '1',
                    start = '38232850', #
                    end =   '38381100', #
                    chrom = '7',
                    strand = '1'
                )
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001078',
                    gene_symbol = 'TRG_locus_(b)', 
                    gene_id = 'ENSG09000001078',
                    transcript_id = 'ENST09000001078',
                    exon_id = 'ENSE09000001078',
                    exon_number = '1',
                    start = '38232850', #
                    end =   '38381100', #
                    chrom = '7',
                    strand = '-1'
                )



                # transcript 1 - exon 2/2 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001100',
                    gene_symbol = 'IGH_locus_(h)', # stjude 
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001100',
                    exon_id = 'ENSE09000001101',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1',
                    
                    start_tr   = '105855906',
                    end_tr     = '105863400',
                    start_ex = '105863198', #
                    end_ex =   '105863400' #
                )

                # transcript 1 - exon 1/2-- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001100',
                    gene_symbol = 'IGH_locus_(h)', # stjude 
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001100',
                    exon_id = 'ENSE09000001100',
                    exon_number = '2',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1',

                    start_tr = '105855906',
                    end_tr   = '105863400',
                    start_ex = '105855906', #
                    end_ex =   '105856217' #
                )


                # transcript 1 - exon 1/2-- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001200',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001200',
                    exon_id = 'ENSE09000001200',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',

                    start_tr = '105855906',
                    end_tr =   '105863400',
                    start_ex = '105855906', #
                    end_ex =   '105856217' #
                )

                # transcript 1 - exon 2/2 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001200',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001200',
                    exon_id = 'ENSE09000001201',
                    exon_number = '2',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863400',
                    start_ex = '105863198', #
                    end_ex  = '105863400', #

                )


                # transcript 2 - exon 1/4-- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001202',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001202',
                    exon_id = 'ENSE09000001202',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',

                    start_tr = '105855906',
                    end_tr =   '105863400',
                    start_ex = '105855906', #
                    end_ex =   '105856217' #
                )

                # transcript 2 - exon 2/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001202',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001202',
                    exon_id = 'ENSE09000001203',
                    exon_number = '2',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863400',
                    start_ex = '105863198', #
                    end_ex  = '105863240', #
                )


                # transcript 2 - exon 3/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001202',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001202',
                    exon_id = 'ENSE09000001204',
                    exon_number = '3',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863400',
                    start_ex = '105863212', # rs74454466 --> -/ACC
                    end_ex  = '105863214', #
                )

                # transcript 2 - exon 4/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001202',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001202',
                    exon_id = 'ENSE09000001205',
                    exon_number = '4',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863400',
                    start_ex = '105863244', #
                    end_ex  = '105863400', #
                )




                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001210',
                    gene_symbol = 'RORC', #
                    gene_id = 'ENSG00000143365',
                    transcript_id = 'ENST09000001210',
                    exon_id = 'ENSE09000001210',
                    exon_number = '1',
                    start = '151806071',
                    end =   '151831872', #
                    chrom = '1',
                    strand = '-1'
                )


                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001250',
                    gene_symbol = 'ALK', #
                    gene_id = 'ENSG00000171094',
                    transcript_id = 'ENST09000001250',
                    exon_id = 'ENSE09000001250',
                    exon_number = '1',
                    start = '29192770',
                    end =   '30000000', #
                    chrom = '2',
                    strand = '-1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001251',
                    gene_symbol = 'NPM1', #
                    gene_id = 'ENSG00000181163',
                    transcript_id = 'ENST09000001251',
                    exon_id = 'ENSE09000001251',
                    exon_number = '1',
                    start = '171380000',
                    end =   '171415000', #
                    chrom = '5',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001252',
                    gene_symbol = 'FGFR1', #
                    gene_id = 'ENSG00000077782',
                    transcript_id = 'ENST09000001252',
                    exon_id = 'ENSE09000001252',
                    exon_number = '1',
                    start = '38411000',
                    end =   '38495000', #
                    chrom = '8',
                    strand = '-1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001253',
                    gene_symbol = 'FGFR2', #
                    gene_id = 'ENSG00000066468',
                    transcript_id = 'ENST09000001253',
                    exon_id = 'ENSE09000001253',
                    exon_number = '1',
                    start = '121451000',
                    end =   '121650000', #
                    chrom = '10',
                    strand = '-1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001254',
                    gene_symbol = 'FGFR3', #
                    gene_id = 'ENSG00000068078',
                    transcript_id = 'ENST09000001254',
                    exon_id = 'ENSE09000001254',
                    exon_number = '1',
                    start = '1790000',
                    end =   '1810000', #
                    chrom = '4',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001255',
                    gene_symbol = 'FGFR4', #
                    gene_id = 'ENSG00000160867',
                    transcript_id = 'ENST09000001255',
                    exon_id = 'ENSE09000001255',
                    exon_number = '1',
                    start = '177085000',
                    end =   '177100000', #
                    chrom = '5',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001256',
                    gene_symbol = 'RET', #
                    gene_id = 'ENSG00000165731',
                    transcript_id = 'ENST09000001256',
                    exon_id = 'ENSE09000001256',
                    exon_number = '1',
                    start = '43043000',
                    end =   '43135000', #
                    chrom = '10',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001257',
                    gene_symbol = 'MYC', #
                    gene_id = 'ENSG00000136997',
                    transcript_id = 'ENST09000001257',
                    exon_id = 'ENSE09000001257',
                    exon_number = '1',
                    start = '127735000',
                    end =   '127744000', #
                    chrom = '8',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001258',
                    gene_symbol = 'EML4', #
                    gene_id = 'ENSG00000143924',
                    transcript_id = 'ENST09000001258',
                    exon_id = 'ENSE09000001258',
                    exon_number = '1',
                    start = '42168000', # chr2:42169350-42332548
                    end =   '42333460', # chr2:42168000-42333460
                    chrom = '2',
                    strand = '1'
                )


                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001260',
                    gene_symbol = 'MALT1', # near MALT1
                    gene_id = 'ENSG00000172175',
                    transcript_id = 'ENST09000001260',
                    exon_id = 'ENSE09000001260',
                    exon_number = '1',
                    start = '58650000', #
                    end =   '58790000', #58754477
                    chrom = '18',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001261',
                    gene_symbol = 'NR4A3', #  NR4A3  Chromosome 9: 99,821,855-99,866,891 forward strand.
                    gene_id = 'ENSG00000119508',
                    transcript_id = 'ENST09000001261',
                    exon_id = 'ENSE09000001261',
                    exon_number = '1',
                    start = '99821200', #
                    end =   '99875400', #58754477
                    chrom = '9',
                    strand = '1'
                )

                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001262',
                    gene_symbol = 'TAF15', #  TAF15  Chromosome 17: 35,713,791-35,864,615 forward strand.
                    gene_id = 'ENSG00000270647',
                    transcript_id = 'ENST09000001262',
                    exon_id = 'ENSE09000001262',
                    exon_number = '1',
                    start = '35800000', #
                    end =   '35864650', #58754477
                    chrom = '17',
                    strand = '1'
                )



                # transcript 2 - exon 1/4-- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001300',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001300',
                    exon_id = 'ENSE09000001300',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',

                    start_tr = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105855906', #
                    end_ex =   '105856217' #
                )

                # transcript 2 - exon 2/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001300',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001300',
                    exon_id = 'ENSE09000001301',
                    exon_number = '2',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105863198', #
                    end_ex  = '105863240', #
                )


                # transcript 2 - exon 3/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001300',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001300',
                    exon_id = 'ENSE09000001302',
                    exon_number = '3',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105863296', # 
                    end_ex  = '105863298', #
                )

                # transcript 2 - exon 4/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001300',
                    gene_symbol = 'IGH_locus_(b)', # stjude 
                    gene_id = 'ENSG09000000014',
                    transcript_id = 'ENST09000001300',
                    exon_id = 'ENSE09000001304',
                    exon_number = '4',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '1',
                    
                    start_tr = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105863241', #
                    end_ex  = '105863248', #
                )


                # transcript 1 - exon 4/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001400',
                    gene_symbol = 'IGH_locus_(h)', # stjude 
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001400',
                    exon_id = 'ENSE09000001403',
                    exon_number = '4',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1',
                    
                    start_tr   = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105863241', #
                    end_ex  = '105863248', #
                )


                # transcript 1 - exon 3/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001400',
                    gene_symbol = 'IGH_locus_(h)', # stjude 
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001400',
                    exon_id = 'ENSE09000001402',
                    exon_number = '3',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1',
                    
                    start_tr   = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105863296', # 
                    end_ex  = '105863298', #
                )

                # transcript 1 - exon 2/4 -- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001400',
                    gene_symbol = 'IGH_locus_(h)', # stjude 
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001400',
                    exon_id = 'ENSE09000001401',
                    exon_number = '2',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1',
                    
                    start_tr   = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105863198', #
                    end_ex =   '105863400' #
                )

                # transcript 1 - exon 1/4-- IGHM
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000001400',
                    gene_symbol = 'IGH_locus_(h)', # stjude 
                    gene_id = 'ENSG09000001014',
                    transcript_id = 'ENST09000001400',
                    exon_id = 'ENSE09000001400',
                    exon_number = '1',
                    start = '105778001', #
                    end =   '106000000', #
                    chrom = '14',
                    strand = '-1',

                    start_tr = '105855906',
                    end_tr =   '105863248',
                    start_ex = '105855906', #
                    end_ex =   '105856217' #
                )



            else:
                print >>sys.stderr,"WARNING: Cannot identify correctly the human genome version!",d[0]

    file(os.path.join(options.output_directory,"custom_genes.txt"),"w").writelines([line+'\n'for line in list_genes])
    #
