#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes a list of fusion genes and their genomic coordinates and converts them
in another coordinates system using liftOver.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2017 Daniel Nicorici

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
import gc
import gzip
import tempfile

def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

def delete_file(some_file):
    if os.path.isfile(some_file) or os.path.islink(some_file):
        os.remove(some_file)


if __name__=='__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It takes a list of fusion genes and their genomic coordinates and converts them in another coordinates system using liftOver."""
    version="%prog 0.11 beta              Author: Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input list of fusion genes and their genome coordinates.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The input list of fusion genes and their genome coordinates.""")

    parser.add_option("--chain","-c",
                      action="store",
                      type="string",
                      dest="chain_filename",
                      help="""The chain files needed by liftOver to do the conversion.""")

    parser.add_option("--tmp_dir",
                      action="store",
                      type="string",
                      dest="tmp_dir",
                      default = None,
                      help="The directory which should be used as temporary directory. By default is the OS temporary directory.")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.chain_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("No inputs and outputs specified!")
        sys.exit(1)

    # read the fusion genes and their fusion coordinates
    coordinates = [line.rstrip("\r\n").split("\t") for line in file(options.input_filename,'r').readlines() if line.rstrip("\r\n")]
    head = coordinates.pop(0)
    if coordinates:
        coord = dict()
        bed = set()
        for line in coordinates:
            k = line[8]
            v = k.split(':')[0:2]
            if v[0].isdigit():
                v[0] = "chr"+v[0]
            elif v[0].lower() == 'mt':
                v[0] = "chrM"
            elif v[0].lower() == 'x':
                v[0] = "chrX"
            elif v[0].lower() == 'y':
                v[0] = "chrY"
            v = "%s\t%s\t%s\n" % (v[0],v[1],str(int(v[1])+1))
            bed.add(v)
            coord[k] = v
            coord[k.replace(":+",":-") if k.endswith("+") else k.replace(":+",":-")] = v
            k = line[9]
            v = k.split(':')[0:2]
            if v[0].isdigit():
                v[0] = "chr"+v[0]
            elif v[0].lower() == 'mt':
                v[0] = "chrM"
            elif v[0].lower() == 'x':
                v[0] = "chrX"
            elif v[0].lower() == 'y':
                v[0] = "chrY"
            v = "%s\t%s\t%s\n" % (v[0],v[1],str(int(v[1])+1))
            bed.add(v)
            coord[k] = v
            coord[k.replace(":+",":-") if k.endswith("+") else k.replace(":+",":-")] = v
        bed = sorted(bed)
        filetemp1 = give_me_temp_filename(options.tmp_dir) #"temp1.bed"
        filetemp2 = give_me_temp_filename(options.tmp_dir) #"temp2.bed"
        filetemp3 = give_me_temp_filename(options.tmp_dir) #"temp3.bed"
        file(filetemp1,'w').writelines(bed)
        cmd = ['liftOver',filetemp1,options.chain_filename,filetemp2,filetemp3]
        cmd = ' '.join(cmd)
        r = os.system(cmd)
        if r:
            print >>sys.stderr, "WARNING: Not able to run '%s'" % (cmd,)
            file(options.output_filename,'w').write('\t'.join(head)+'\n')
        else:
            unlifted = set(file(filetemp3,'r').readlines())
            lifted = [line for line in bed if line not in unlifted]
            newbed = file(filetemp2,'r').readlines()
            n = len(lifted)
            m = len(newbed)
            if n!=m:
                print >>sys.stderr, "WARNING: Lost into BED translation!"
                print >>sys.stderr, lifted
                print >>sys.stderr, unlifted
                print >>sys.stderr, newbed
                file(options.output_filename,'w').write('\t'.join(head)+'\n')
            else:
                translate = dict()
                for i in xrange(n):
                    v = newbed[i]
                    if v.startswith('chrM\t'):
                        v = "MT"+v[4:]
                    elif v.startswith("chr"):
                        v = v[3:]
                    v = ':'.join(v.split('\t')[0:2])+':'
                    translate[lifted[i]] = v
                for k in coord.keys():
                    x = coord[k]

                    y = translate.get(x,'not-converted')
                    if k.endswith(":+"):
                        y = y + "+"
                    elif k.endswith(":-"):
                        y = y + "-"
                    coord[k] = y
                # now to final conversion step using the fusion file
                fusions = []
                fusions.append(head)
                for line in coordinates:
                    line[8] = coord[line[8]]
                    line[9] = coord[line[9]]
                    fusions.append(line)
                file(options.output_filename,'w').writelines(['\t'.join(line)+'\n' for line in fusions])
        delete_file(filetemp1)
        delete_file(filetemp2)
        delete_file(filetemp3)
    else:
        file(options.output_filename,'w').write('\t'.join(head)+'\n')
    #
