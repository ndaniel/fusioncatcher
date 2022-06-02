#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a list of chromosomal coordinates of fusion genes and it labels them if they use exon-exon borders for the fusion junction.

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2022 Daniel Nicorici

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
import gzip




#
#
#
def add_line(aline,adict):
    chrom = aline[0]
    pstart = int(aline[3])
    pend = int(aline[4])
    strand = aline[6]
    if pstart > pend:
        (pstart,pend) = (pend,pstart)
    ids = [l.replace('"','').replace("'","").strip().split(' ') for l in aline[8].split(";") if l]
    ids = dict([l for l in ids if len(l) == 2])
    g_id = ids["gene_id"]
    if g_id not in adict:
        adict[g_id] = set()
    adict[g_id].add((chrom,strand,pstart,pend))


#
#
#
def exonexon(gtf_file,
             input_file,
             output_file,
             verbose = True):


    if verbose:
        print >>sys.stderr,"Reading the list of fusion genes..."
    # read the list of fusion genes and their chromosomal positions
    data = [line.rstrip('\r\n').split('\t') for line in file(input_file,'r').readlines() if line.rstrip('')]
    header = data.pop(0)
    fusion1 = [tuple([line[10]] + line[8].split(':')) for line in data]
    fusion2 = [tuple([line[11]] + line[9].split(':')) for line in data]
    myg = set([line[10] for line in data]+[line[11] for line in data])

    if myg:

        # read and pre-process the GTF file
        if verbose:
            print >>sys.stderr,"Parsing the GTF file..."
        g = []
        if gtf_file == '-':
            g = [line.rstrip('\r\n').split("\t") for line in file(sys.stdin,"r").readlines() if (not line.startswith("#")) and line.rstrip()]
        elif gtf_file.endswith('.gz'):
            z = gzip.open(gtf_file,"r")
            g = [line.rstrip('\r\n').split("\t") for line in z.readlines() if (not line.startswith("#")) and line.rstrip()]
            z.close()
        else:
            g = [line.rstrip('\r\n').split("\t") for line in file(gtf_file,"r").readlines() if (not line.startswith("#")) and line.rstrip()]



        # get all exons per gene as a dictionary
        if verbose:
            print >>sys.stderr,"Building the database of exons..."
        exon = dict()
        q = 0
        for line in g:
            if line[2] == 'gene':
                q = q + 1
            flag = False
            x = line[8].partition(";")
            if x[0] and x[0].startswith("gene_id "):
                x = x[0][9:-1]
            if x not in myg:
                continue
            if line[2] == 'exon':
                add_line(line,exon)

        
        if verbose:
            print >>sys.stderr,"Checking exon borders..."
        for i in xrange(len(fusion1)):
            f1 = fusion1[i]
            f2 = fusion2[i]
            
            g1 = f1[0]
            c1 = f1[1]
            p1 = int(f1[2])
            s1 = f1[3]
            
            g2 = f2[0]
            c2 = f2[1]
            p2 = int(f2[2])
            s2 = f2[3]
            
            # check 5 prime partner
            e = exon.get(g1,None)
            ok5 = False
            if e:
                x = [1 for v in e if v[0] == c1 and v[1] == s1 and ((s1 == "+" and v[3] == p1) or (s1 == "-" and v[2] == p1))]
                if x:
                    ok5 = True
            # check 3 prime partner
            e = exon.get(g2,None)
            ok3 = False
            if e:
                x = [1 for v in e if v[0] == c2 and v[1] == s2 and ((s2 == "+" and v[2] == p2) or (s2 == "-" and v[3] == p2))]
                if x:
                    ok3 = True

            if ok5 and ok3:
                comma = ',' if data[i][2] else ''
                data[i][2] = data[i][2] + comma + 'exon-exon'
        
        data.insert(0,header)
        file(output_file,'w').writelines(['\t'.join(line)+'\n' for line in data])
    else:
        file(output_file,'w').write('\t'.join(header)+'\n')



################################################################################
################################################################################
################################################################################
def main():

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options]"

    description = """It takes as input a list of chromosomal coordinates of fusion genes and it labels them if they use exon-exon borders for the fusion junction."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2022 Daniel Nicorici

"""



    version = "%prog 0.04 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
                            )



    parser.add_option("-g","--gtf",
                      action = "store",
                      type = "string",
                      dest = "gtf_filename",
                      help = """The input GTF file containing the genome annotation.""")

    parser.add_option("-i","--input",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file containing the fusion (chromosomal) coordinates for each fusion genes.""")


    parser.add_option("-o","--output",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output file where the frame predictions are written. """)


    parser.add_option("-q", "--quiet",
                      action = "store_false",
                      dest = "verbose",
                      default = True,
                      help = "Do not print status messages to stdout.")

    parser.add_option("-a","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.gtf_filename and
            options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")


    #
    exonexon(options.gtf_filename,
             options.input_filename,
             options.output_filename,
             options.verbose)


if __name__ == '__main__':
    main()
    
    
#    
