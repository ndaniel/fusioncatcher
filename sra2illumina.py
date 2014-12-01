#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It convets the FASTQ file (PHRED-33 qualities and SRA read names) downloaded
from Short Read Archive (SRA) to Illumina FASTQ file (PHRED-64 Illumina v1.5
and Illumina read names).



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2014 Daniel Nicorici

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
import tempfile
import optparse
import phred
import shutil
import errno


def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

def int2str(x = 0, digits = 10):
    return '0' * int(digits - len(w)) + str(x)


def sra2illumina(input_file,
                 output_file,
                 tag_read = None,
                 tag='',
                 phred_conversion = False,
                 operation = 'change',
                 tmp_dir = None,
                 size_read_buffer = 10**8):
    """
    It converts the FASTQ file (PHRED-33 qualities and SRA read names) downloaded
    from Short Read Archive (SRA) to Illumina FASTQ file (PHRED-64 Illumina v1.5
    and Illumina read names).
    """
    temp_file = None
    if phred_conversion:
        temp_file = give_me_temp_filename(tmp_dir)
    else:
        temp_file = output_file

    read_name = file(input_file,'r').readline().rstrip('\r\n')
    sra = False
    e = read_name.partition(" ")[0]
    if read_name.startswith('@') and ( not(e.endswith('/1') or e.endswith('/2'))):
        sra = True

    if operation == 'change' or sra:
        fid = open(input_file,'r')
        fod = open(temp_file,'w')
        i = 0
        r = 0
        while True:
            lines = fid.readlines(size_read_buffer)
            if not lines:
                break
            n = len(lines)
            for j in xrange(n):
                r = r + 1
                i = i + 1
                if i == 1:
                    if tag_read:
                        lines[j] = '@%s%s%s\n' % (tag_read ,int2str(r,12) , tag)
                    else: # if there is no tag_read then the original SRA id is left
                        lines[j] = '%s%s\n' % (lines[j][:-1].partition(" ")[0], tag)
                    #lines[j] = lines[j].rstrip('\r\n').upper().split(' ')[1]+tag+'\n'
                elif i == 3:
                    lines[j] = "+\n"
                elif i == 4:
                    i = 0
            fod.writelines(lines)
        fid.close()
        fod.close()
        if phred_conversion == '64':
            phred.fq2fq(temp_file,'sanger',output_file,'illumina-1.5',tmp_dir = tmp_dir)
            os.remove(temp_file)
        elif phred_conversion == '33':
            phred.fq2fq(temp_file,'auto-detect',output_file,'sanger',tmp_dir = tmp_dir)
            os.remove(temp_file)
    else:
        print "No changes are done!"
        if os.path.isfile(output_file):
            os.remove(output_file)
        if operation == 'soft':
            if os.path.islink(input_file):
                linkto = os.readlink(input_file)
                os.symlink(linkto,ooutput_file)
            else:
                os.symlink(input_file,output_file)
        elif operation == 'hard':
            linkto = input_file
            if os.path.islink(input_file):
                linkto = os.readlink(input_file)
            try:
                os.link(linkto,output_file)
            except OSError as er:
                print >>sys.stderr,"WARNING: Cannot do hard links ('%s' and '%s')!" % (linkto,output_file)
                shutil.copyfile(linkto,output_file)
#                if er.errno == errno.EXDEV:
#                    # they are on different partitions
#                    # [Errno 18] Invalid cross-device link
#                    shutil.copyfile(linkto,output_file)
#                else:
#                    print >>sys.stderr,"ERROR: Cannot do hard links ('%s' and '%s')!" % (linkto,output_file)
#                    print >>sys.stderr,er
#                    sys.exit(1)

        elif operation == 'copy':
            shutil.copyfile(input_file, output_file)
        else:
            print >>sys.stderr, "ERROR: unknown operation of linking!", operation
            sys.exit(1)
#
#
#
if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It converts the FASTQ file (PHRED-33 qualities and SRA read names) downloaded from Short Read Archive (SRA) to Illumina FASTQ file (PHRED-64 Illumina v1.5 and Illumina read names).
Note: The command arguments input_1 and output_1 should be used in the same time. The same applies to input_2 and output_2. Also all four may be used in the same time."""
    version = "%prog 0.15 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input_1",
                      action="store",
                      type="string",
                      dest="input_1_filename",
                      help="""The input FASTQ file downloaded from SRA.""")

    parser.add_option("--input_2",
                      action="store",
                      type="string",
                      dest="input_2_filename",
                      help="""The input FASTQ file downloaded from SRA.""")

    parser.add_option("--output_1",
                      action="store",
                      type="string",
                      dest="output_1_filename",
                      help="""The output FASTQ file in Illumina format.""")

    parser.add_option("--output_2",
                      action="store",
                      type="string",
                      dest="output_2_filename",
                      help="""The output FASTQ file in Illumina format.""")

    parser.add_option("--tag_read_name",
                      action="store",
                      type="string",
                      dest="tag_read_name",
                      default = "",
                      help="""This tag is added to the beginning of each read name in order to make them unique in case that multiple FASTQ files need to be merged.""")

    parser.add_option("--phred64",
                      action="store_true",
                      dest="phred64_conversion",
                      default = False,
                      help="""If it is used then the PHRED-33 qualities from the input SRA are converted to PHRED-64 qualities.""")

    parser.add_option("--phred33",
                      action="store_true",
                      dest="phred33_conversion",
                      default = False,
                      help="""If it is used then the PHRED-64 qualities from the input SRA are converted to PHRED-33 qualities.""")

    parser.add_option("--no12",
                      action="store_true",
                      dest="no12",
                      default = False,
                      help="""By default to all reads /1 and /2 will be added to the reads ids. By using this no adding is done.""")

    parser.add_option("--tmp_dir",
                  action = "store",
                  type = "string",
                  dest = "tmp_dir",
                  default = None,
                  help = "The directory which should be used as temporary directory. By default is the OS temporary directory.")

    choices = ('soft','hard','copy','change')
    parser.add_option("--link",
                      action = "store",
                      choices = choices,
                      dest = "link",
                      default = 'change',
                      help = "If it is set to 'change' always the reads names from the "+
                             "input files are changed and no actual checking is done "+
                             "to see if it is indeed a FASTQ file containing a SRA-like reads names. "+
                             "For 'soft', 'hard', and 'copy' a checking is done to see if the "+
                             "reads names look like SRA names and if yes then their names will "
                             "be changed to new ones. In this case "+
                             "it creates a link from the output file to the input "+
                             "file of type ("+','.join(choices)+") in case that no "+
                             "operation is done on the input file. "+
                             "The choices are: "+','.join(choices)+ ". "+
                             "Default is '%default'.")


    (options,args)=parser.parse_args()

    # validate options
    if not ((options.input_1_filename and
            options.output_1_filename) or
            (options.input_2_filename and
            options.output_2_filename)
            ):
        parser.print_help()
        parser.error("Missing command arguments!")

    phred_conversion = None
    if options.phred33_conversion:
        phred_conversion = '33'
    elif options.phred64_conversion:
        phred_conversion = '64'


    if options.input_1_filename and options.output_1_filename:
        sra2illumina(options.input_1_filename,
                     options.output_1_filename,
                     options.tag_read_name,
                     '' if options.no12 else '/1',
                     phred_conversion,
                     options.link,
                     options.tmp_dir
                     )
    if options.input_2_filename and options.output_2_filename:
        sra2illumina(options.input_2_filename,
                     options.output_2_filename,
                     options.tag_read_name,
                     '' if options.no12 else '/2',
                     phred_conversion,
                     options.link,
                     options.tmp_dir
                     )
