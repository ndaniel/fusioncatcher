#!/mnt/software/bin/python2.7-real
# -*- coding: utf-8 -*-
"""
It sorts the input file (text tab separated file) based on the specified columns.
It works like SELECT * SORT BY in SQL.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2019 Daniel Nicorici

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
sort -t $'\t' -k 1n,1 -k 2n,2 -k4rn,4 -k3,3 <my-file>

sort the column 1 (numeric), 2, 4 (numeric and reverse order), 3
sed 1d file | sort >output
date | tee -a directory_listing
cp header.orig header       #just to save your header for future use
cat sortedfile >> header

HEADER=` head -1 /applis/projets/COP_DEV/data/HMC/fichseq/ee1ESP0LOT16AUX_01_177_0_DS.NEWHEADER`
sed -i "1i$HEADER" /applis/projets/COP_DEV/data/HMC/fichseq/ee1ESP0LOT16AUX_01_177_0_DS

"""
import os
import sys
import optparse
import tempfile
import gc
import multiprocessing

######### Functions ############


def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

def delete_file(some_file):
    if os.path.isfile(some_file) or os.path.islink(some_file):
        os.remove(some_file)

def empty(a_file):
    f = True
    if (os.path.isfile(a_file) or os.path.islink(a_file)) and os.path.getsize(a_file) != 0:
        f = False
    return f

#
# sort
#
def sort_columns(input_filename = '-',
                 output_filename = '-',
                 columns=None,
                 header=False,
                 ignore_case=False,
                 unique = False,
                 tmp_dir=None,
                 buffer_size = '80%',
                 parallel = multiprocessing.cpu_count(),
                 compress_program = None):
    """
    It sorts the input file (text tab separated file) based on the specified columns.
    It works like SELECT * ORDER BY in SQL.
    """
    import locale

    locale.setlocale(locale.LC_ALL, 'C')

    # check options suppported by SORT command
    sh1 = give_me_temp_filename(tmp_dir)
    sh2 = give_me_temp_filename(tmp_dir)
    # check options suppported by SORT command
    sort_parallel = False
    r = os.system("sort --help | grep 'parallel' > '%s'" % (sh1,))
    if (not r) and (not empty(sh1)) and len(file(sh1,'r').readlines()) == 1:
        sort_parallel = True
    delete_file(sh1)
    # check options suppported by SORT command
    sort_buffer = False
    r = os.system("sort --help | grep 'buffer-size' > '%s'" % (sh1,))
    if (not r) and (not empty(sh1)) and len(file(sh1,'r').readlines()) == 1:
        sort_buffer = True
    delete_file(sh1)
    # check options suppported by SORT command
    sort_compress = False
    if compress_program:
        r = os.system("sort --help | grep 'compress-program' > '%s' ; %s --help 2>/dev/null | grep -i 'compress' > '%s'" % (sh1,compress_program,sh2))
        if (not r) and ((not empty(sh1)) and len(file(sh1,'r').readlines()) == 1 and
            (not empty(sh2)) and len(file(sh2,'r').readlines()) >= 1):
            sort_compress = True
        delete_file(sh1)
        delete_file(sh2)

    # treat the case when the input file is coming from the standard input
    fin = input_filename.strip('"').strip("'")
    if fin == '-':
        fin = give_me_temp_filename(tmp_dir)
        fod = open(fin,'w')
        fid = sys.stdin
        while True:
            lines = fid.readlines(10**8)
            if not lines:
                break
            fod.writelines(lines)
        fod.close()

    fon = output_filename.strip('"').strip("'")
    if fon == '-':
        fon = give_me_temp_filename(tmp_dir)

    if header:
        header_saved = file(fin,'r').readline()
        file(output_filename,'w').write(header_saved)
    else:
        file(output_filename,'w').write('')

    # process the type of the column, numeric, or string
    first_line = file(fin,'r').readline()
    if first_line:
        nc=len(file(fin,'r').readline().rstrip('\r\n').split('\t'))#read first line in order to find out the number of columns
        if columns:
            columns=columns.strip().lower()
            if columns=='d':
                columns=','.join([str(i+1)+'d' for i in range(nc)])
            elif columns=='n':
                columns=','.join([str(i+1)+'n' for i in range(nc)])
            elif columns=='nd' or columns=='dn':
                columns=','.join([str(i+1)+'nd' for i in range(nc)])
        else:
            columns=','.join([str(i+1) for i in range(nc)])

        # extra parameters
        extra = ""

        if sort_buffer and buffer_size and buffer_size != 'no' and buffer_size != 'none':
                extra = extra + ' --buffer-size=' + str(buffer_size) + ' '
        if sort_parallel and parallel and parallel > 1:
                extra = extra + ' --parallel=' + str(parallel) + ' '
        if sort_compress and compress_program and compress_program.lower() != 'no' and compress_program.lower() != 'none':
                extra = extra + ' --compress-program='+compress_program + ' '

        # processing the input columns
        columns = ['-k '+el+','+el.replace('n','').replace('r','') for el in columns.replace('d','r').split(',')]
        comd = "-s -t '\t' "+" ".join(columns)
        if ignore_case:
            comd = "-f "+comd
        if unique:
            comd = "-u "+comd
        if tmp_dir:
            comd = "-T '"+tmp_dir+"' "+comd
        if header:
            comd = "LC_ALL=C sed 1d '" + fin + "' | LC_ALL=C sort " + extra + comd + " >> '" + output_filename + "'"
        else:
            comd = "LC_ALL=C sort " + extra + comd + " '" + fin + "' >> '" + output_filename + "'"
        r = os.system(comd)
        if r != 0:
            print >>sys.stderr, "ERROR (sort_ttdb.py) while running:"
            print >>sys.stderr, comd
            sys.exit(1)

    if input_filename == '-':
        os.remove(fin)

    if output_filename == '-':
        fod = sys.stdout
        fid = open(fon,'r')
        while True:
            gc.disable()
            lines = fid.readlines(10**8)
            gc.enable()
            if not lines:
                break
            fod.writelines(lines)
        fid.close()
        os.remove(fon)

    #locale.setlocale(locale.LC_ALL, '')




#######################################################################
#######################################################################
#######################################################################



########### START ###########################




if __name__ == "__main__":
    #command line parsing
    usage="%prog [options] arg"
    description="It sorts the lines in a text tab separated file based on the given columns."
    version="%prog 0.11"

    parser=optparse.OptionParser(usage = usage, description = description, version = version)

    parser.add_option("--input",
                  action="store",
                  type="string",
                  dest="input_filename",
                  default='-',
                  help="The filename of the input text tab separated file used for sorting. For standard input '-' can be used. If missing it defaults to standard input.")

    parser.add_option("--input_header",
                  action="store_true",
                  dest="input_header",
                  default = False,
                  help="If it is used it means that the input file has a header which will be ignored during sorting.")

    parser.add_option("--input_columns",
                  action="store",
                  type="string",
                  dest="input_columns",
                  help="The index number of the column(s) from first input file to be used for sorting. The first column is numerotated as 1. The columns will be treatead as strings by default. If one wants a column to be treated as numeric should add append N to the column number, e.g. 1n. The columns will be ordered in ascending order by default in the order given here. If one wants a column to be sorted in descending order should append a D to the column number, e.g. 1d or 1nd. More than one column can be specified if comma is used as separator. If no columns are specified then all columns from left to right are sorted. By specifying only N then all columns will be treated as numeric. By specifying only D all columns will be sorted in descending order. By specifying only ND all columns will be treated as numeric and sorted in descending order.")

    parser.add_option("--output",
                  action="store",
                  type="string",
                  dest="output_filename",
                  default='-',
                  help="The output text tab separated file which contains the results of the sorting. For standard output '-' can be used. If missing it defaults to standard output.")

    parser.add_option("--tmp_dir",
                  action="store",
                  type="string",
                  dest="tmp_dir",
                  default = None,
                  help="The directory which should be used as temporary directory. By default is the OS temporary directory.")

    parser.add_option("--ignore_case",
                  action="store_true",
                  default = False,
                  dest="ignore_case",
                  help="If it is used it means that when the sorting is done the case will be ignore, i.e. 'BILL' is the same with 'bill'.")

    parser.add_option("--unique",
                  action="store_true",
                  default = False,
                  dest="unique",
                  help="Only those lines are written to the output which have the columns unique.")


    parser.add_option("--buffer-size",
                  action = "store",
                  type = "string",
                  default = "80%",
                  dest = "buffer_size",
                  help = "The main buffer size which is passed further to GNU sort command. For more see '--buffer-size' of GNU sort command. Default is '%default'.")

    parser.add_option("--parallel",
                  action = "store",
                  type = "int",
                  default = multiprocessing.cpu_count(),
                  dest = "parallel",
                  help = "The number of parallel processes to be used for sorting. This is passed further to GNU sort command. For more see '--parallel' of GNU sort command. Default is '%default'.")

    parser.add_option("--tmp-compress-program",
                  action = "store",
                  type = "string",
                  default = "lzop",
                  dest = "compress_program",
                  help = "The compress program to be used when reading/writing the temporary files. If no compression program should be used then 'none' or 'no' should be specified. For more see '--compress-program' of GNU sort command. Default is '%default'.")


    (options,args)=parser.parse_args()

    # validate command line arguments
    if not (
            options.input_columns
            ):
        parser.print_help()
        parser.error("Missing columns which should be used for sorting the input file!")


    # sort the columns
    sort_columns(
                 options.input_filename,
                 options.output_filename,
                 options.input_columns,
                 options.input_header,
                 options.ignore_case,
                 options.unique,
                 options.tmp_dir,
                 options.buffer_size,
                 options.parallel,
                 options.compress_program
                 )
