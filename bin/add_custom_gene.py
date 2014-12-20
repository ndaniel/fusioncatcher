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

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """Add the custom human genes which are missing from the Ensembl database."""
    version = "%prog 0.13 beta"

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

    print "Add the custom human genes which are missing from the Ensembl database..."

    if options.organism.lower() == 'homo_sapiens':

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

            elif d[0].lower().find('grch38') !=-1:
                print "Found version GRCh38 human genome version!"
                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = '',
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

                # coordinates valid only for GRCh38
                add(outdir = options.output_directory,
                    protein_id = '',
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

#                coordinates valid only for GRCh38
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
                    protein_id = '',
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
            else:
                print >>sys.stderr,"WARNING: Cannot identify correctly the human genome version!",d[0]

    #
