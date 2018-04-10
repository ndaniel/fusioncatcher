#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It produces a very short summary of fusion genes and transcripts found.



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

################################################################################
################################################################################
################################################################################

if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It generates a very short summary of gene/transcript fusions found."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input report containg detailed information about fusion genes found.""")

    parser.add_option("--viruses","-u",
                      action = "store",
                      type = "string",
                      dest = "input_viruses_filename",
                      help = """The input report containg detailed information about viruses/bacteria/found found.""")

    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output summary of gene/transcript fusions found.""")




    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")
        sys.exit(1)


    top_virus = ''
    if options.input_viruses_filename and os.path.exists(options.input_viruses_filename):
        viruses = [line.rstrip("\r\n").strip().split(" ") for line in file(options.input_viruses_filename,'r') if line.rstrip("\r\n")]
        viruses.pop(0)
        if viruses:
            top_virus = [v[1] for v in viruses if v[1].find('virus') != -1 and v[0] != '1']
            tv = []
            for e in top_virus:
                ee = e.split(",_")[0]
                ee = ee.replace(",_complete_genome","").replace("_complete_genome","").replace("_complete_wild_type_genome","").replace("genome","").replace("complete","").replace("wild","").replace("wild_type","").replace("__","").replace("__","").replace("_"," ")
                tv.append(ee)
            if tv:
                # remove duplicates
                uniq = set()
                ntv = []
                for t in tv:
                    if t in uniq:
                        continue
                    else:
                        ntv.append(t)
                        uniq.add(t)
                top_virus = "', '".join(ntv[0:3]) # just top 3 viruses

    data = [line.upper().rstrip("\r\n").split("\t") for line in file(options.input_filename,'r').readlines() if line.rstrip("\r\n")]
    header = data.pop(0)



    if data:
        fusions_genes = set()
        for r in data:
            s = "%s--%s" % (r[0],r[1])
            if r[0] > r[1]:
                s = "%s--%s" % (r[1],r[0])
            fusions_genes.add(s)

        fusions_transcripts = set([ "%s---%s---%s---%s" % (r[0],r[1],r[8],r[9]) for r in data])

        # find reciprocal fusions
        rf = set(["%s--%s" % (r[0],r[1]) for r in data])
        reciprocal = set()
        for el in rf:
            rec = el.split("--")
            rev = "%s--%s" % (rec[1],rec[0])
            if rev in rf:
                reciprocal.add(rev)
                reciprocal.add(el)

        # known fusion
        known = set()
        for r in data:
            f1 = "%s--%s" % (r[1],r[0])
            f2 = "%s--%s" % (r[0],r[1])
            if r[2].lower().find('known') != -1 or r[2].lower().find('tcga') != -1 or r[2].lower().find('cosmic') != -1 or r[2].lower().find('cell_lines') != -1 or r[2].lower().find('prostates') != -1 or r[2].lower().find('pancreases') != -1 or r[2].lower().find('chimerdb') != -1:
                known.add(f1)
                known.add(f2)

        # questionable
        questionable = set()
        for r in data:
            f1 = "%s--%s" % (r[1],r[0])
            f2 = "%s--%s" % (r[0],r[1])
            if r[2].lower().find('healthy') != -1 or r[2].lower().find('conjoing') != -1 or r[2].lower().find('hpa') != -1 or r[2].lower().find('banned') != -1 or r[2].lower().find('paralogs') != -1 or r[2].lower().find('1000genomes') != -1 or r[2].lower().find('cortex') != -1:
                questionable.add(f1)
                questionable.add(f2)

        # readthrough
        readthrough = set()
        for r in data:
            f1 = "%s--%s" % (r[1],r[0])
            f2 = "%s--%s" % (r[0],r[1])
            if r[2].lower().find('readthrough') != -1:
                readthrough.add(f1)
                readthrough.add(f2)


        found = []
        bag = set()
        for r in data:
            s = "%s--%s" % (r[0],r[1])
            t = s[:]
            if r[0] > r[1]:
                s = "%s--%s" % (r[1],r[0])
            if s in bag:
                continue
            else:
                found.append(t)
                bag.add(s)

        results = []
        results.append("Very short summary of found candidate fusion genes\n")
        results.append("==================================================\n\n")
        if top_virus:
            results.append("The input sample contains sequencing reads mapping on: '%s'.\n\n" % (top_virus,))
        results.append("Found %d fusion gene(s), which are as follows:\n" % (len(fusions_genes),))
        for line in found:
            label = []
            if line in reciprocal:
                label.append('reciprocal fusion')
            if line in known:
                label.append('already known fusion')
            if line in questionable:
                label.append('probably false positive')
            if line in readthrough:
                label.append('readthrough')
            if label:
                label = "  ("+'; '.join(label)+")"
            else:
                label = ''
            results.append("  * %s%s\n" % (line,label))
        results.append("\nFound %d fusion transcript(s).\n\n" % (len(fusions_transcripts),))

        results.append("For more detailed information regarding these candidate fusions, see text file 'final-list_candidate-fusion-genes.txt'.\n")

        file(options.output_filename,"w").writelines(results)
    else:
        results = []
        results.append("Very short summary of found candidate fusion genes\n")
        results.append("==================================================\n\n")
        if top_virus:
            results.append("The input sample contains sequencing reads mapping on '%s'.\n\n" % (top_virus,))
        results.append("No fusion genes found.\n\n")
        results.append("No fusion transcripts found.\n")

        file(options.output_filename,"w").writelines(results)


#
