#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes a FASTQ file generated by Illumina Solexa pipeline and removes
(truncates them to length one) the reads marked as bad by Illumina.



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
'--skip_blat'. Fore more information regarding BLAT please see its license.

Please, note that FusionCatcher does not require BLAT in order to find
candidate fusion genes!

This file is not running/executing/using BLAT.
"""

#
"""
Example:
It changes:
@GQWE8:57:C00T6ABXX:2:1101:1233:2230 1:N:0:CTTGTA

to

@GQWE8:57:C00T6ABXX:2:1101:1233:2230/1
"""
import os
import sys
import optparse
import string
import shutil
import gc
import errno
import gzip

def reads_from_fastq(f_name,size_read_buffer=10**8):
    fid = None
    if f_name == '-':
        fid = sys.stdin
    elif f_name.lower().endswith('.gz'):
        fid = gzip.open(f_name,'r')
    else:
        fid = open(f_name,'r')
    piece = [None,None,None,None]
    i = 0
    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            i = i + 1
            piece[i-1] = a_line
            if i == 4:
                yield (piece[0],piece[1],piece[3])
                piece = [None,None,None,None]
                i = 0
    fid.close()



class lines_to_file:

    def __init__(self,file_name,size_buffer=10**8):
        self.file_name=file_name
        if file_name:
            self.file_handle = open(file_name,'w')
        self.size_buffer = size_buffer
        self.data = []
        self.size = 0

    def add_line(self,line):
        gc.disable()
        self.data.append(line)
        gc.enable()
        self.size = self.size + len(line)
        if self.size > self.size_buffer:
            self.__write_buffer()

    def __write_buffer(self):
        self.file_handle.writelines(self.data)
        self.size = 0
        self.data = []
    def is_filename_valid(self):
        if self.file_name:
            return True
        else:
            return False
    def close(self):
        if self.is_filename_valid():
            if self.data:
                self.__write_buffer()
            self.file_handle.close()
            self.file_name = None
    def __del__(self):
        self.close()

if __name__=='__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It takes a FASTQ file generated by Illumina Solexa pipeline and removes
(truncates them to length one) the reads (which look like this @GQWE8:57:C00T6ABXX:2:1101:1233:2230 1:N:0:CTTGTA) marked as bad by Illumina. Example: It changes  to @GQWE8:57:C00T6ABXX:2:1101:1233:2230/1"""
    version = "%prog 0.10 beta              Author: Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action = "store",
                      type = "string",
                      dest = "input_filename",
                      help = """The input file (in the newer Solexa FASTQ format, i.e. version 1.8 or newer) containing the short reads to be processed.""")

    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_filename",
                      help = """The output FASTQ file containing the short reads which are marked as good by Illumina.""")


    choices = ('soft','hard','copy')
    parser.add_option("--link",
                      action = "store",
                      choices = choices,
                      dest = "link",
                      default = 'soft',
                      help = """It creates a link from the output file to the input file of type ("""+','.join(choices)+""") in case that no operation is done on the input file. Default is '%default'.""")


    (options, args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("ERROR: No inputs and outputs specified!")


    # check if the FASTQ file is in format produced by CASAVA pipeline version 1.8
    t = file(options.input_filename,"r").readline()
    if not t:
        print >>sys.stderr,"ERROR: The input file '%s' is empty!" % (options.input_filename,)
        sys.exit(1)
    t = t.rstrip()

    if not t.startswith("@"):
        print >>sys.stderr,"ERROR: The input file '%s' is not in FASTQ file format!" % (options.input_filename,)
        print >>sys.stderr,"The read names in the input file look like '%s'" % (t,)
        sys.exit(1)


    flag_illumina = False
    tt = t.split(" ")
    if len(tt) >= 2:
        t0 = tt[0].split(":")
        t1 = tt[1].split(":")
        if len(t1) == 4 and (t1[0] == '1' or t1[0] == '2') and (t1[1] == 'N' or t1[1] == 'Y'):
            flag_illumina = True

    flag_length = False
    if t.lower().find(' length') != -1:
        flag_length = True

    #
    if not (flag_illumina or flag_length):
        print >>sys.stderr,"No changes are done (all reads are kept)!"
        if os.path.isfile(options.output_filename):
            os.remove(options.output_filename)
        if options.link == 'soft':
            if os.path.islink(options.input_filename):
                linkto = os.readlink(options.input_filename)
                os.symlink(linkto,options.output_filename)
            else:
                os.symlink(options.input_filename,options.output_filename)
        elif options.link == 'hard':
            linkto = options.input_filename
            if os.path.islink(options.input_filename):
                linkto = os.readlink(options.input_filename)
            try:
                os.link(linkto,options.output_filename)
            except OSError as er:
                print >>sys.stderr,"WARNING: Cannot do hard links ('%s' and '%s')!" % (linkto,options.output_filename)
                shutil.copyfile(linkto,options.output_filename)
#                if er.errno == errno.EXDEV:
#                    # they are on different partitions
#                    # [Errno 18] Invalid cross-device link
#                    shutil.copyfile(linkto,options.output_filename)
#                else:
#                    print >>sys.stderr,"ERROR: Cannot do hard links ('%s' and '%s')!" % (linkto,options.output_filename)
#                    print >>sys.stderr,er
#                    sys.exit(1)

        elif options.link == 'copy':
            shutil.copyfile(options.input_filename, options.output_filename)
        else:
            print >>sys.stderr, "ERROR: unknown operation of linking!", options.link
            sys.exit(1)
    else:
#        print "The input file is in FASTQ format compatible with Illumina pipeline CASAVA version 1.8!"
        data = lines_to_file(options.output_filename)

        i = 0
        j = 0
        for reads in reads_from_fastq(options.input_filename):
            i = i + 1
            r = reads[0].partition(" ")
            if r[2].startswith("1:Y:") or r[2].startswith("2:Y:"):
                data.add_line("%s\n%s\n+\n%s\n" % (r[0],reads[1][:1],reads[2][:1]))
                j = j + 1
            elif r[1]:
                data.add_line("%s\n%s+\n%s" % (r[0],reads[1],reads[2]))
            else:
                data.add_line("%s%s+\n%s" % (reads[0],reads[1],reads[2]))
        print >>sys.stderr,"%d (%.3f%%) reads out of %d total reads were removed due to being labeled as bad quality by Illumina!" % (j,(float(j)/i*100),i)
        data.close()
        #
