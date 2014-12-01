#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It trim the short reads from 5-end or 3-end.



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

if __name__ == '__main__':

    #command line parsing

    usage="%prog [options]"
    description="""It trim the short reads from 5-end or 3-end."""
    version="%prog 0.10 beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""Input FASTQ file which contanins the reads which will be trimmed.""")

    parser.add_option("--trim_size",
                      action="store",
                      type="int",
                      dest="trim_size",
                      help="""The size of part which is trimmed from each read. It should not be used together with --size!""")

    parser.add_option("--final_size",
                      action="store",
                      type="int",
                      dest="size",
                      help="""The size of read part which is left after trimming. If the read is shorter then no trimming is performed.""")

    choices = ( '5', '3' )
    parser.add_option("--trim_end",
                      action="store",
                      type="choice",
                      choices = choices,
                      dest="end",
                      help="""The end of the short read where the trimming will be done. The choices are ["""+','.join(choices)+"""].""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output FASTQ file containing the trimmed reads.""")



    ( options, args ) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename and
            options.end and
            (options.trim_size or options.size)
            ):
        parser.print_help()
        parser.error("One of the options has not been specified.")

    if options.trim_size and options.size:
        parser.error("Options --trim_size and --final_size are mutually exclusive!")

    fid = open(options.input_filename,'r')
    fod = open(options.output_filename,'w')
    size_buffer = 10**8
    m = -1
    i = 0
    if options.trim_size:
        trim = options.trim_size
        if options.end == '3':
            while True:
                lines = fid.readlines(size_buffer)
                if not lines:
                    break
                n = len(lines)
                for j in xrange(n):
                    i = i + 1
                    if i == 2:
                        lines[j] = lines[j].rstrip('\r\n')[:-trim].upper()+'\n'
                    elif i == 4:
                        lines[j] = lines[j].rstrip('\r\n')[:-trim].upper()+'\n'
                        i = 0
                    elif i == 3:
                        lines[j] = '+\n'
                fod.writelines(lines)
        elif options.end == '5':
            while True:
                lines = fid.readlines(size_buffer)
                if not lines:
                    break
                n = len(lines)
                for j in xrange(n):
                    i = i + 1
                    if i == 2:
                        lines[j] = lines[j][trim:].upper()
                    elif i == 4:
                        lines[j] = lines[j][trim:].upper()
                        i = 0
                    elif i == 3:
                        lines[j] = '+\n'
                fod.writelines(lines)
    elif options.size:
        size = options.size
        m = -1
        if options.end == '3':
            while True:
                lines = fid.readlines(size_buffer)
                if not lines:
                    break
                n = len(lines)

                for j in xrange(n):
                    i = i + 1
                    if i == 2:
                        lines[j] = lines[j].rstrip('\r\n')
                        m = len(lines[j])
                        if m > size:
                            lines[j] = lines[j][:size].upper()+'\n'
                        else:
                            lines[j] = lines[j].upper()+'\n'
                    elif i == 4:
                        if m > size:
                            lines[j] = lines[j].rstrip('\r\n')[:size]+'\n'
                        i = 0
                    elif i == 3:
                        lines[j] = '+\n'
                fod.writelines(lines)
        elif options.end == '5':
            size = options.size
            m = -1
            while True:
                lines = fid.readlines(size_buffer)
                if not lines:
                    break
                n = len(lines)
                for j in xrange(n):
                    i = i + 1
                    if i == 2:
                        lines[j] = lines[j].rstrip('\r\n')
                        m = len(lines[j])
                        if m > size:
                            lines[j] = lines[j][m-size:].upper()+'\n'
                        else:
                            lines[j] = lines[j].upper()+'\n'
                    elif i == 4:
                        if m > size:
                            lines[j] = lines[j].rstrip('\r\n')[m-size:]+'\n'
                        i = 0
                    elif i == 3:
                        lines[j] = '+\n'
                fod.writelines(lines)
    fod.close()
    fid.close()
