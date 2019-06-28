#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add the custom human genes which are missing from the Ensembl database.



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


def int2str(x,n=3):
    x = str(x)
    return '0' * int(n - len(x)) + x


def generate(ig_start, ig_end, ig_chrom, ig_interval, ig_overlap, ig_special, ig_id):

    ig_template = """
                add(outdir = options.output_directory,
                    protein_id = 'ENSP090%%S%%000%%COUNT%%',
                    gene_symbol = '%%ID%%_locus_(%%COUNT%%)',
                    gene_id = 'ENSG090%%S%%000%%COUNT%%',
                    transcript_id = 'ENST090%%S%%000%%COUNT%%',
                    exon_id = 'ENSE090%%S%%000%%COUNT%%',
                    exon_number = '1',
                    start = '%%START%%',
                    end =   '%%END%%',
                    chrom = '%%CHROM%%',
                    strand = '%%STRAND%%'
                )
    """


    intervals = range(ig_start,ig_end,ig_interval)

    last = intervals[-1]
    if ig_end - last < 5000:
        intervals = intervals[:-1]

    n = len(intervals)
    j = -1
    for i,interval in enumerate(intervals):
        start = interval
        end = interval + ig_interval + ig_overlap
        if i == n - 1:
            end = ig_end
        t = ig_template
        t = t.replace("%%S%%",str(ig_special)).replace("%%START%%",str(start)).replace("%%END%%",str(end)).replace("%%ID%%",ig_id).replace("%%CHROM%%",str(ig_chrom))
        j = j + 1
        print t.replace("%%COUNT%%",int2str(j,3)).replace("%%STRAND%%","1")
        j = j + 1
        print t.replace("%%COUNT%%",int2str(j,3)).replace("%%STRAND%%","-1")



################################################################################

ig_start = 88846000
ig_end = 90381500
ig_chrom = 2
ig_interval = 200000
ig_overlap = 300
ig_special = '11'
ig_id = 'IGK'
generate(ig_start, ig_end, ig_chrom, ig_interval, ig_overlap, ig_special, ig_id)

ig_start = 21549840
ig_end = 22563979
ig_chrom = 14
ig_interval = 200000
ig_overlap = 300
ig_special = '12'
ig_id = 'TRA'
generate(ig_start, ig_end, ig_chrom, ig_interval, ig_overlap, ig_special, ig_id)


ig_start = 142274400
ig_end = 142821000
ig_chrom = 7
ig_interval = 200000
ig_overlap = 300
ig_special = '13'
ig_id = 'TRB'
generate(ig_start, ig_end, ig_chrom, ig_interval, ig_overlap, ig_special, ig_id)

ig_start = 38232850
ig_end = 38381100
ig_chrom = 7
ig_interval = 200000
ig_overlap = 300
ig_special = '14'
ig_id = 'TRG'
generate(ig_start, ig_end, ig_chrom, ig_interval, ig_overlap, ig_special, ig_id)

ig_start = 105556000
ig_end = 106883700
ig_chrom = 14
ig_interval = 200000
ig_overlap = 300
ig_special = '15'
ig_id = 'IGH'
generate(ig_start, ig_end, ig_chrom, ig_interval, ig_overlap, ig_special, ig_id)


