#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It runs BLAT in parallel by dividing the input file into equal parts.

Date: September 9, 2010.

Requirements:
- BLAT
- Python version >= 2.6
- BioPython version >= 1.50



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

This file is executing BLAT.

"""
import sys
import os
import multiprocessing
import subprocess
import time
import tempfile
import Bio.SeqIO
import math
import gc
import shutil

def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    delete_file(ft_name)
    return ft_name

def delete_file(some_file):
    if os.path.isfile(some_file) or os.path.islink(some_file):
        os.remove(some_file)
    elif os.path.isdir(some_file):
        shutil.rmtree(some_file)

def quote(txt):
    r = txt.strip()
    if txt.find(' ') != -1:
        r = '"%s"' % (txt,)
    return r

#
#
#
if __name__ == '__main__':

    cpus = multiprocessing.cpu_count()
    print >>sys.stderr,cpus,"processor(s) found!"

    cmd = sys.argv[1:]
    if not cmd:
        print >>sys.stderr,"It runs BLAT in parallel by dividing the input file into equal parts."
        print >>sys.stderr,"The temporary directory where the splitting is done can be specified using the option '--tmp_dir', e.g. '--tmp_dir=/some/temp/dir/'. If it is not specified the OS temporary directory is used."
        print >>sys.stderr,"The blat directory where the blat executable is placed, '--blat_dir', e.g. '--blat_dir=/some/blat/dir/'."
        print >>sys.stderr,"The option '--cpus' specifies the number of  CPUs to be used, e.g. '--cpus=10'. If it is not specified then all the CPUs found will be used."
        print >>sys.stderr,"The option '--filter-fusion' forces that all the lines in the PSL output are filtered according to finding gene fusions. If it is not specified then no filtering is done."


    blat_dir = [el for el in cmd if el.startswith('--blat_dir=')]
    if blat_dir:
        blat_dir = blat_dir[0][11:]
    else:
        blat_dir = None
    
    _BT_ = ""
    if blat_dir:
        _BT_ = blat_dir.rstrip("/")+"/"

    tmp_dir = [el for el in cmd if el.startswith('--tmp_dir=')]
    if tmp_dir:
        tmp_dir = tmp_dir[0][10:]
    else:
        tmp_dir = None

    cmd_cpus = [el for el in cmd if el.startswith('--cpus=')]
    if cmd_cpus:
        cmd_cpus = int(cmd_cpus[0][7:])
        if cpus > 0:
            cpus = cmd_cpus
            print >>sys.stderr,cpus,"CPU(s) will be used!"
    else:
        cmd_cpus = None

    filtered = [el for el in cmd if el.startswith('--filter-fusion')]
    if filtered:
        filtered = True
    else:
        filtered = False
    # remove the --tmp_dir and --cpus from the commands to be pass to BLAT
    cmd = [el for el in cmd if ((not el.startswith('--tmp_dir=')) and
                                (not el.startswith('--cpus=')) and
                                (not el.startswith('--blat_dir=')) and
                                (not el.startswith('--filter-fusion')))]

    if not tmp_dir:
        j = -1
        for i,el in enumerate(cmd):
            if el.startswith('--tmp_dir'):
                j = i
                break
        if j!=-1:
            tmp_dir = cmd[j+1]
            cmd.pop(j+1)
            cmd.pop(j)
        else:
            tmp_dir = None


    filenames = [el for el in cmd if not el.startswith('-')]
    cmd = [el for el in cmd if el.startswith('-')] # keep only the ones that starts with '-'

    if len(filenames)<3:
        print >>sys.stderr,"Error: not enough input and output files!"
        sys.exit(1)

    input_filename = filenames[1]
    database_filename = filenames[0]
    output_filename = filenames[2]


    print >>sys.stderr,"Counting the records in the input file..."
    input_handle = open(input_filename, "rU")
    input_seq_iterator = Bio.SeqIO.parse(input_handle, "fasta")
    count = len([1 for r in input_seq_iterator])
    input_handle.close()
    print >>sys.stderr," -",count,"records found in the input file!"
    if count > 0:
        print >>sys.stderr,"Splitting the input file into",cpus,"pieces..."
        list_input_temp_files = [give_me_temp_filename(tmp_dir) for i in range(cpus)]
        list_output_temp_files = [give_me_temp_filename(tmp_dir) for i in range(cpus)]
        pipes = [give_me_temp_filename(tmp_dir) for i in range(cpus)]

        input_handle = open(input_filename, "rU")
        input_seq_iterator = Bio.SeqIO.parse(input_handle, "fasta")
        empty_flag = []
        size_block = math.ceil(float(count) / float(cpus))
        if count <= 5*cpus:
            size_block = count
        i = -1
        j = -1
        seq = []
        for record in input_seq_iterator:
            i = i + 1
            seq.append(record)
            if ( (i+1) % size_block == 0) or (i + 1 == count):
                j = j + 1
                print >>sys.stderr," - writing part",j+1
                output_handle = open(list_input_temp_files[j], "w")
                cc = Bio.SeqIO.write(seq, output_handle, "fasta")
                output_handle.close()
                print >>sys.stderr,"       wrote",cc,"sequences"
                if seq:
                    empty_flag.append(False)
                else:
                    empty_flag.append(True)
                seq = []
        input_handle.close()
        if len(empty_flag) != cpus:
            empty_flag = empty_flag + [True] * (cpus - len(empty_flag))


        proc = []
        print >>sys.stderr,"Launching BLAT in parallel..."
        for i in xrange(j+1):
            if empty_flag[i]:
                continue
            parameters = [_BT_+'blat'] + cmd + [quote(database_filename), quote(list_input_temp_files[i]), quote(list_output_temp_files[i])]
            if filtered:
                parameters = [os.path.abspath(os.path.dirname(__file__))+'/blat-filter-fusion.sh',
                              _BT_ if _BT_ else '-',
                              quote(database_filename),
                              quote(list_input_temp_files[i]),
                              quote(pipes[i]),
                              quote(list_output_temp_files[i])] + cmd


            print >>sys.stderr,'-->JOB:'+str(i+1)+'-----------------------------------------------------------------'
            print >>sys.stderr,' '.join(parameters)

            # original
            #p = subprocess.Popen(parameters)
            # fix https://github.com/ndaniel/fusioncatcher/issues/62
            p = subprocess.Popen(parameters,close_fds=True)

            proc.append(p)
            time.sleep(1) # seconds
        print >>sys.stderr,'-------------------------------------------------------------------------'
        print >>sys.stderr,"Waiting for BLAT to finish running..."

        # original
        #for p in proc:
        #    p.communicate()
        # fix https://github.com/ndaniel/fusioncatcher/issues/62
        while proc:
            proc.pop().communicate()
            
        print >>sys.stderr,"BLAT finished running."
        print >>sys.stderr,'-------------------------------------------------------------------------'
        print >>sys.stderr,"Joining BLAT's output files..."
        fod = open(output_filename,'w')
        flag = False
        for i in xrange(j+1):
            if empty_flag[i] or (not os.path.exists(list_output_temp_files[i])):
                continue
            fid = open(list_output_temp_files[i],'r')
            while True:
                lines = fid.readlines(10**8)
                if not lines:
                    break
                if not lines[-1].endswith('\n'):
                    lines[-1] = lines[-1]+'\n'
                fod.writelines(lines)
                flag = True
            fid.close()
        if not flag:
            fod.write('')
        fod.close()

        for i in range(cpus):
            delete_file(list_input_temp_files[i])
            delete_file(list_output_temp_files[i])
            delete_file(pipes[i])
    print >>sys.stderr,"Done."


    #
