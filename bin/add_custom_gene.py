#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add the custom human genes which are missing from the Ensembl database.



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
        strand = ''):
    #
    g = [line.rstrip("\r\n").split("\t") for line in file("genes.txt","r").readlines() if line.rstrip("\r\n")]
    gid = [(el[0], el[1], el[2], el[3], el[4]) for el in g if el[0] == gene_id]
    if gene_id and (not gid):
        print "Gene %s not found in the database! Gene %s and transcript %s added into database!" % (gene_id,gene_id,transcript_id)
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
                start,
                end,
                exon_number,
                start,
                end,
                start,
                end,
                strand,
                chrom))
    elif gid:
        print "Gene %s found already in the database!" % (gene_id,)
        if len(gid) == 1:
            gid = gid.pop(0)
        else:
            print "  Error: Too many genes found in 'genes.txt'!",gid
            sys.exit(1)

        e = [line.rstrip("\r\n").split("\t") for line in file("exons.txt","r").readlines() if line.rstrip("\r\n")]

        if end != gid[1] or start != gid[2] or strand != gid[3] or chrom != gid[4]:
            print "  Gene %s requires changes in the entire database! New transcript %s added for this gene!" % (gene_id,transcript_id)
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
                    start,
                    end,
                    exon_number,
                    start,
                    end,
                    start,
                    end,
                    strand,
                    chrom] )
            e = e_rest + e_target
            file(os.path.join(outdir,'exons.txt'),'w').writelines(['\t'.join(line)+'\n' for line in e])
        else:
            e_tr = [line for line in e if line[2] == transcript_id and line[1] == gene_id]
            if e_tr:
                print "  Transcript %s already in the database and it will not be added again!" % (transcript_id,)
            else:
                print "  New transcript %s added for already present gene %s!" % (transcript_id,gene_id)
                file(os.path.join(outdir,'exons.txt'),'a').write(
                    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        protein_id,
                        gene_id,
                        transcript_id,
                        exon_id,
                        start,
                        end,
                        exon_number,
                        start,
                        end,
                        start,
                        end,
                        strand,
                        chrom))
    list_genes.append(gene_id)

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Add the custom human genes which are missing from the Ensembl database."""
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

    print "Add/change the human genes which have mistakes or are missing from the Ensembl database..."

    file(os.path.join(options.output_directory,"custom_genes.txt"),"w").write('')

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
                    start = '105778720',
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
                    end = '1220000', # 1267000 #1393000
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
                    start = '1220001', #'1213300'
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

                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000009',
                    gene_symbol = 'IGK_locus_', # stjude
                    gene_id = 'ENSG09000000009',
                    transcript_id = 'ENST09000000009',
                    exon_id = 'ENSE09000000009',
                    exon_number = '1',
                    start = '88846000', #'1213300'
                    end =   '89154500', # 54380200
                    chrom = '2',
                    strand = '1'
                )

                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000010',
                    gene_symbol = 'IGK_locus__', # stjude
                    gene_id = 'ENSG09000000010',
                    transcript_id = 'ENSG09000000009',
                    exon_id = 'ENSE09000000010',
                    exon_number = '1',
                    start = '89154501', #'1213300'
                    end =   '89463000', # 54380200
                    chrom = '2',
                    strand = '1'
                )
                #coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000011',
                    gene_symbol = 'IGK_locus___', # stjude
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
                    gene_symbol = 'IGK_locus____', # stjude
                    gene_id = 'ENSG09000000012',
                    transcript_id = 'ENST09000000012',
                    exon_id = 'ENSE09000000012',
                    exon_number = '1',
                    start = '89951251', #'1213300'
                    end =   '90381500', # 54380200
                    chrom = '2',
                    strand = '-1'
                )

                #
                # IGH locus -- split in several pieces
                #
                # IGH_locus: 14::+:chr14:105,556,000-106,883,700

                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = 'ENSP09000000013',
                    gene_symbol = 'IGH_locus_', # stjude
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
                    gene_symbol = 'IGH_locus__', # stjude
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
                    gene_symbol = 'IGH_locus___', # stjude
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
                    gene_symbol = 'IGH_locus____', # stjude
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
                    gene_symbol = 'IGH_locus_____', # stjude
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
                    gene_symbol = 'IGH_locus______', # stjude
                    gene_id = 'ENSG09000000018',
                    transcript_id = 'ENST09000000018',
                    exon_id = 'ENSE09000000018',
                    exon_number = '1',
                    start = '106663101', #
                    end =   '106883700', #
                    chrom = '14',
                    strand = '1'
                )


#                add(outdir = options.output_directory,
#                    protein_id = 'ENSP09000000019',
#                    gene_symbol = 'SWSAP1_', # stjude # overlaps EPOR on opposite strand
#                    gene_id = 'ENSG09000000019',
#                    transcript_id = 'ENST09000000019',
#                    exon_id = 'ENSE09000000019',
#                    exon_number = '1',
#                    start = '11377000', #
#                    end =   '11394000', #
#                    chrom = '19',
#                    strand = '1'
#                )

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


            else:
                print >>sys.stderr,"WARNING: Cannot identify correctly the human genome version!",d[0]

    file(os.path.join(options.output_directory,"custom_genes.txt"),"w").writelines([line+'\n'for line in list_genes])
    #
