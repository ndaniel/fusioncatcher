#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It removes the reads from the input FASTQ file which do not form a pair (i.e. all the reads in the output FASTQ file are forming a pair).
It is assumed that reads are ending with /1 and /2.



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

import optparse
import os
import sys
import sort_ttdb
import tempfile
import gc
import multiprocessing

#
def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

#
def line_from_file(f_name, size_read_buffer = 10**8):
    fid = open(f_name,'r')
    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            p = a_line.rstrip('\r\n').split('\t')
            yield p
    fid.close()

#
def reads_from_fastq(f_name, size_read_buffer = 10**8):
    fid = sys.stdin
    if f_name != '-':
        fid = open(f_name,'r')
    piece = []
    i = 0
    while True:
        gc.disable()
        lines = fid.readlines(size_read_buffer)
        gc.enable()
        if not lines:
            break
        for a_line in lines:
            i = i + 1
            piece.append(a_line)
            if i == 4:
                yield piece
                i = 0
                piece = []

    if piece:
        yield piece
    fid.close()

def add_newlines(u):
    return ["%s\n" % (el,) for el in u]

###################################################################
###################################################################
###################################################################

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It removes the reads from the input FASTQ file which do not form a pair (i.e. all the reads in the output FASTQ file are forming a pair).
It is assumed that reads are ending with /1 and /2."""
    version="%prog 0.10 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""Input FASTQ file which contains the reads which will be processed for filtering.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTQ file containing the reads which form pairs.""")

    parser.add_option("--log",
                      action="store",
                      type="string",
                      dest="log_filename",
                      help="""The log file where the percentage of removed reads (as percentage) is written.""")

    parser.add_option("--interleaved","-x",
                      action = "store_true",
                      default = False,
                      dest="interleaved",
                      help="""If specified then the input FASTQ file is considered to have the reads interleaved (i.e. reads which form a pair to be one afte each other).""")

    parser.add_option("-p", "--processes",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 0,
                      help = """Number of parallel threads/processes/CPUs to be used for computations. In case of value 0 then the program will use all the CPUs which are found. The default value is %default.""")


    parser.add_option("--tmp_dir",
                  action="store",
                  type="string",
                  dest="tmp_dir",
                  default = None,
                  help="The directory which should be used as temporary directory. By default is the OS temporary directory.")

    ( options, args ) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        parser.error("One of the options has not been specified.")
        sys.exit(1)

    cpus = options.processes
    if cpus == 0:
        cpus = multiprocessing.cpu_count()

    if options.interleaved:
        # process the reads
        fod = open(options.output_filename,'w')
        last_read = None
        data = []
        limit = 10**5
        i = 0
        ii = 0
        t = 0
        for chunk in reads_from_fastq(options.input_filename):
            t = t + 1
            if not last_read:
                last_read = chunk[:]
                continue
            else:
                r1 = last_read[0].rstrip('\r\n')
                r2 = chunk[0].rstrip('\r\n')
                if r1[-2] == '/' and r1[-1] != r2[-1] and r1[:-1] == r2[:-1]:
                    gc.disable()
                    data.extend(last_read)
                    data.extend(chunk)
                    gc.enable()
                    last_read = None
                    i = i + 2
                    if i > limit:
                        fod.writelines(data)
                        data = []
                        ii = ii + i
                        i = 0
                else:
                    last_read = chunk[:]
        if data:
            fod.writelines(data)
            ii = ii + i
            data = []
        fod.close()
        print >>sys.stderr,"%d single short reads removed from a total of %d short reads." %(t-ii,t)

    else:
        # convert the fastq file into a text tab separated file
        buffer_size = 10**8
        temp_file_1 = give_me_temp_filename(tmp_dir = options.tmp_dir)
        fid = sys.stdin
        if options.input_filename != '-':
            fid = open(options.input_filename,'r')
        fod = open(temp_file_1,'w')
        i = 0
        while True:
            lines = fid.readlines(buffer_size)
            if not lines:
                break
            for j in xrange(len(lines)):
                i = i + 1
                if i != 4:
                    lines[j] = lines[j].rstrip('\r\n') + '\t'  # replace new line character with tab
                else:
                    i = 0
            fod.writelines(lines)
        fod.close()
        fid.close()

        # sort the temp_file_1 by first column which is the read name
        temp_file_2 = give_me_temp_filename(tmp_dir = options.tmp_dir)
        sort_ttdb.sort_columns(temp_file_1,
                               temp_file_2,
                               columns = '1', # short read name
                               header = False,
                               ignore_case = False,
                               tmp_dir = options.tmp_dir,
                               parallel = cpus
                               )
        os.remove(temp_file_1)

        # process the reads
        fod = open(options.output_filename,'w')
        last_read = None
        data = []
        limit = 10**5
        i = 0
        ii = 0
        t = 0
        for chunk in line_from_file(temp_file_2):
            t = t + 1
            if not last_read:
                last_read = chunk[:]
                continue
            else:
                r1 = last_read[0]
                r2 = chunk[0]
                if r1[-2] == '/' and r1[-1] != r2[-1] and r1[:-1] == r2[:-1]:
                    gc.disable()
                    data.extend(add_newlines(last_read))
                    data.extend(add_newlines(chunk))
                    gc.enable()
                    last_read = None
                    i = i + 2
                    if i > limit:
                        fod.writelines(data)
                        data = []
                        ii = ii + i
                        i = 0
                else:
                    last_read = chunk[:]
        if data:
            fod.writelines(data)
            data = []
            ii = ii + i
        fod.close()
        os.remove(temp_file_2)
        if t != 0:
            x = float(t-ii)/float(t)
            print >>sys.stderr,"%d (%.3f%%) single short reads removed from a total of %d short reads." %(t-ii,x*100,t)
            if options.log_filename:
                file(options.log_filename,'w').write("%f" %(x,))
        else:
            print >>sys.stderr,"WARNING (remove_single_reads.py): The input file is empty???!!!"

        #
