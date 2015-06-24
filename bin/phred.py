#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It converts a FASTQ file from one PHRED quality score to another one.



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

#
"""
More info

http://en.wikipedia.org/wiki/FASTQ_format

Support PHRED quality scores for conversion are:
sanger -- refers to Sanger style FASTQ files which encode PHRED qualities using an ASCII offset of 33.
solexa -- refers to early Solexa/Illumina style FASTQ files which encode Solexa qualities using an ASCII offset of 64.
illumina -- refers to recent Solexa/Illumina style FASTQ files (from pipeline version 1.3+) which encode PHRED qualities using an ASCII offset of 64.
"""
import sys
import os
import shutil
import optparse
import Bio.SeqIO
import gzip
import errno

def give_me_temp_filename(tmp_dir = None):
    import tempfile
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

def quals_from_fastq(file_name, first = 10000000, size_read_buffer = 10**8):
    fid = open(file_name,'r')
    qual = ''
    idx = 0
    i = 0
    while True:
        lines = fid.readlines(size_read_buffer)
        if (not lines) or (i > first):
            break
        for a_line in lines:
            idx = idx + 1
            if idx == 4:
                qual = a_line.rstrip('\r\n')
                yield qual
                qual = ''
                idx = 0
                i = i + 1
    fid.close()

def detect_fastq_format(file_name):
    #
    t = 'sanger' # default format
    mm = 126
    for q in quals_from_fastq(file_name):
        m = min(map(ord,q))
        if m < mm:
            mm = m
            if mm < 59:
                break
    if mm < 59:
        t = 'sanger'
    elif mm < 64:
        t = 'solexa'
    elif mm <= 126:
        t = 'illumina'
    return t

def fq2fq(f_in,
          f_in_type,
          f_out,
          f_out_type,
          link = 'symbolic',
          tmp_dir = None):
    """
    It converts qualities in a FASTQ file.

    f_in_type     ('sanger','solexa','illumina','auto-detect')
    f_out_type    ('sanger','solexa','illumina','illumina-1.5')

    link          ('soft','hard','copy')

    tmp_dir       temporary directory
    """
    if f_in_type.lower() == 'auto-detect':
        # detect the input FASTQ format type
        f_in_type = detect_fastq_format(f_in)
        print  >>sys.stderr,"Auto-detect found "+f_in_type.upper()+" FASTQ format!"
    fit = 'fastq-'+f_in_type
    fot = ''
    if f_out_type == 'illumina-1.5':
        fot = 'fastq-illumina'
    else:
        fot = 'fastq-'+f_out_type
    if fit == fot and fit != '-':
        # input type is same as output type
        if os.path.isfile(f_out) or os.path.islink(f_out):
            os.remove(f_out)
        if link == 'soft':
            if os.path.islink(f_in):
                linkto = os.readlink(f_in)
                os.symlink(linkto,f_out)
            else:
                os.symlink(f_in,f_out)
        elif link == 'hard':
            linkto = f_in
            if os.path.islink(f_in):
                linkto = os.readlink(f_in)
            try:
                os.link(linkto,f_out)
            except OSError as er:
                print >>sys.stderr,"WARNING: Cannot do hard links ('%s' and '%s')!" % (linkto,f_out)
                shutil.copyfile(linkto,f_out)
#                if er.errno == errno.EXDEV:
#                    # they are on different partitions
#                    # [Errno 18] Invalid cross-device link
#                    shutil.copyfile(linkto,f_out)
#                else:
#                    print >>sys.stderr,"ERROR: Cannot do hard links ('%s' and '%s')!" % (linkto,f_out)
#                    print >>sys.stderr,er
#                    sys.exit(1)

        elif options.link == 'copy':
            shutil.copyfile(f_in, f_out)
        else:
            print >>sys.stderr, "ERROR: unknown operation of linking!", link
            sys.exit(1)

    else:
        if hasattr(Bio.SeqIO,'convert'):
        
            fin = f_in
            if f_in == "-":
                fin = sys.stdin
            elif f_in.lower().endswith('.gz'):
                fid = gzip.open(f_in,'r')
            else:
                fin = open(f_in,'r')
                

            fout = f_out
            if f_out == "-":
                fout = sys.stdout
            elif f_out.lower().endswith('.gz'):
                fout = gzip.open(f_out,'r')
            else:
                fout = open(f_out,'w')


            counts = Bio.SeqIO.convert(f_in, fit, f_out, fot)
            
            
            fin.close()
            fout.close()
            
        else:
            print  >>sys.stderr,"Bio.SeqIO.convert() not supported!"
            print  >>sys.stderr,"Trying to go around it!"
            if f_in == '-' or f_out == '-':
                print  >>sys.stderr,"ERROR: BioPython library from Python is tool old! Please, upgrade it!"
                sys.exit(1)
            input_handle = open(f_in, "rU")
            output_handle = open(f_out, "w")
            sequences = Bio.SeqIO.parse(input_handle, fit)
            counts = Bio.SeqIO.write(sequences,output_handle, fot)
            output_handle.close()
            input_handle.close()
        if f_out_type == 'illumina-1.5':
            ftemp = give_me_temp_filename(tmp_dir)
            fi = open(f_out,'rb')
            fo = open(ftemp,'wb')
            size_buffer = 10**8
            i = 0
            while True:
                lines = fi.readlines(size_buffer)
                if not lines:
                    break
                lines = [line.replace('@','B').replace('A','B') if ((i+j+1)%4 == 0) else line for (j,line) in enumerate(lines)]
                i = i + len(lines)
                fo.writelines(lines)
            fi.close()
            fo.close()
            shutil.move(ftemp,f_out)

        print  >>sys.stderr,"Converted %i records" % counts

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It converts a quality scores in a FASTQ file.

Support PHRED quality scores for conversion are:
sanger -- refers to Sanger style FASTQ files which encode PHRED qualities using an ASCII offset of 33.
solexa -- refers to early Solexa/Illumina style FASTQ files which encode Solexa qualities using an ASCII offset of 64.
illumina -- refers to recent Solexa/Illumina style FASTQ files (from pipeline version 1.3+) which encode PHRED qualities using an ASCII offset of 64.
illumina-1.5 -- refers to recent Solexa/Illumina style FASTQ files (from pipeline version 1.5+) which encode PHRED qualities using an ASCII offset of 64 (smallest quality is B).
    """
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in FASTQ format.""")

    choices = ('sanger','solexa','illumina','auto-detect')
    parser.add_option("--input_type",
                      action="store",
                      type="choice",
                      choices = choices,
                      dest="input_type",
                      default = 'auto-detect',
                      help="""Type quality encoding used in the FASTQ input file. The choices are: ["""+','.join(choices)+"""]""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTQ file.""")

    choices = ('sanger','solexa','illumina','illumina-1.5')
    parser.add_option("--output_type",
                      action="store",
                      type="choice",
                      choices = choices,
                      dest="output_type",
                      help="""Type quality encoding used in the FASTQ output file. The choices are: ["""+','.join(choices)+"""]""")

    parser.add_option("--tmp_dir",
                  action = "store",
                  type = "string",
                  dest = "tmp_dir",
                  default = None,
                  help = "The directory which should be used as temporary directory. By default is the OS temporary directory.")

    choices = ('soft','hard','copy')
    parser.add_option("--link",
                      action = "store",
                      choices = choices,
                      dest = "link",
                      default = 'soft',
                      help = """It creates a link from the output file to the input file of type ("""+','.join(choices)+""") in case that no operation is done on the input file. Default is '%default'.""")

    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.input_type and
            options.output_filename and
            options.output_type
            ):
        parser.print_help()
        parser.error("Missing input and output arguments!")

    # running
    print >>sys.stderr,"Converting..."
    fq2fq(options.input_filename,
          options.input_type,
          options.output_filename,
          options.output_type,
          options.link,
          options.tmp_dir)

    print >>sys.stderr,"The end."
