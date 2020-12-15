#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

WORFLOW
-------

It provides a framework for running workflows/pipelines in Python, i.e. running
external programs using command line.

================================================================================

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2020 Daniel Nicorici

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


================================================================================

The main goal is to compute checksums for:
- command line, which needs to be executed,
- the parameters which are passed as command line arguments,
- files which are given as inputs or outputs, and
- directories which are given as inputs or outputs (the checksum is computed for
  all files in the given directory, the subdirectoris and files contained in the
  subdirectories are ignored, the files in a directory are always considered in
  alphabetical order).

The checksums are saved by default in the file 'checksums.txt' in the current path.
This can be changed by changing for example the line

job = pipeline()

to

job = pipeline(log_filename='/other_path/log.txt',
               checksum_filename='/some_path/checksums.txt')

in the example 1 or 2 from below.

If a command line has the same overall checksum as a already computed checksum
it will skipped and not executed again. This is useful in situation when the
execution of a program takes very long time, e.g. hours or days, and has to be
executed several times even that the inputs and outputs have not changed.

If a task/job has only inputs and parameters always it will be run.



EXAMPLES:
========

Example 1:
----------

from workflow import pipeline # use this Python pipeline library
job = pipeline(log_filename='log.txt', checksums_filename='checksums.txt') # initialize the pipeline

job.add('notepad.exe',kind='program') # specify the main program
job.add('','aha.txt',kind='input') " specify the command line arguments
job.run() # run the "notepad.exe aha.txt"

This runs 'notepad.exe aha.txt' and it computes the checksum for file 'aha.txt'
and for string 'notepad.exe aha.txt'. The overall checksum is saved in
'checksums.txt' file and the log of the pipeline execution is saved in
'log.txt' file.


Example 2:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.add('dir',kind='program')
job.add('>','list.txt',kind='output')
job.run()

It runs "dir > list.txt".


Example 3:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.add('dir',kind='program')
job.add('|','',kind='parameter')
job.add('tee','list.txt',kind='output')
job.run()

It runs "dir | tee list.txt".


Example 4:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
list_files=['file1.txt','file2.txt','file3.txt']
job.add('zip',kind='program')
job.add('-9','',kind='parameter')
job.add('','archive.zip',kind='output')
job.add_list('',list_files,kind='input')
job.run()

it runs "zip -9 archive.zip file1.txt file2.txt file3.txt"


Example 5:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
list_files_input=['file1i.txt','file2i.txt','file3i.txt']
list_files_output=['file1o.txt','file2o.txt','file3o.txt','file4o.txt']
job.add('some_program',kind='program')
job.add('-m','1',kind='parameter')
job.add_list('--input',list_files_input,kind='input',list_separator=',')
job.add_list('--output',list_files_output,kind='output',list_separator=',')
job.run()

it runs "some_program -m 1 --input file1i.txt,file2i.txt,file3i.txt --output file1o.txt,file2o.txt,file3o.txt,file4o.txt"

Example 6:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.add('zip',kind='program')
job.add('-9','',kind='parameter')
job.add('','some_path/archive.zip',kind='output',dest='zip_archive')
job.add('','*',kind='input')
job.run()
job.add('unzip',kind='program')
job.add('',job.zip_archive,kind='input')
job.run()

it runs:
zip -9 some_path/archive.zip *
unzip some_path/archive.zip

and later in the python code the path "some_path/archive.zip" can be accessed
thru 'job.zip_archive'.

Example 7:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.add('',kind='program')
if job.run():
    # some python code
    for i in xrange(1,10):
        print i

it runs the python code:
for i in xrange(1,10):
    print i
as a workflow step (i.e. it has a step count, it can be skipped later from
execution with 'start_step').


Example 8:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.link('foo_1.txt','foo_2.txt')

it creates a soft link from 'foo_1.txt' to 'foo_2.txt'.


Example 9:
----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.link('foo_1.txt','foo_2.txt',temp_path='yes')

it moves 'foo_1.txt' to 'foo_2.txt'.

Example 10:
-----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.clean('foo.txt',temp_path='yes')

it deletes 'foo.txt'.

Example 11:
-----------

from workflow import pipeline # use this Python pipeline library
job=pipeline() # initialize the pipeline
job.write("This is a test message!")

it outputs to the standard output the text "This is a test message!" and also
adds it to the log_filename.


Example 12:
-----------

from workflow import pipeline # use this Python pipeline library
job.add('some_program',kind='program')
job.add('-m','some_parameter.txt',kind='parameter',from_file='yes') # 'some_parameter.txt' contains this text 47
job.add('--output','foo.txt',kind='output')
job.run()

it runs the following:

some_progra -m 47 --output foo.txt

This is useful in cases that one wants to skip this step later from execution.

Example 13:
-----------

from workflow import pipeline # use this Python pipeline library
job = pipeline() # initialize the pipeline

# step 1
job.add('ls',kind='program')
job.add('-alh',kind='parameter')
job.add('>','list_files.txt',kind='output')
job.run()

# step 2 (mock step)
if job.run():
    # step 3
    job.add('cp',kind='program')
    job.add('','list_files.txt',kind='input')
    job.add('','copy_list_files.txt',kind='output')
    job.run()


it executes the step 1 and in case the the step 1 is not skipped from execution
(e.g. when start_step >= 4) then step 2 (which does nothing) and step 3 are
executed. This is useful in cases when there are operations with files and
parameters which are passed to the pipeline are not fully supported by the
pipeline/workflow.

Example 14:
-----------

from workflow import pipeline # use this Python pipeline library
from random import randint # random numbers

job = pipeline() # initialize the pipeline

# step 1
job.add('ls',kind='program')
job.add('-alh',kind='parameter')
job.add('>','list_files.txt',kind='output')
job.run()

i = randint(0,10) # generate a random number between 0 and 10 (including 0 and 10)

# step 2 (IF with memory)
if job.iff(i == 2, id = "random-test"):
    # step 3
    job.add('cp',kind='program')
    job.add('','list_files.txt',kind='input')
    job.add('','copy_list_files.txt',kind='output')
    job.run()


When running first time the pipeline lets assume that randint(0,10) gives 2.
In this case the step 1, step 2 (here the job.iff is TRUE becasue i == 2), and
step 3 are executed.
In case of a re-run from step 1 (or step 2) lets assume that randint(0,10) gives
5 but the IF proceeds as the last run of the pipeline (when
'if job.iff(i == 2, id = "random-test"):' was TRUE because previously randint(0,10)
value 2 and not 5), that is job.iff gives again the value TRUE because is a re-run).
Therefore the step 1, step 2, and step 3 are executed even though 'i' has value 5.
This is useful in cases when there are operations with files which do not exist
at the moment when there is a re-run.

"""



import os
import sys
import datetime
import hashlib
import time
import shutil
import math
import errno
import subprocess
import tempfile

#import multiprocessing


def _expand(*p):
    """
    Expands a path.
    """
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))


def _islink(alink = None):
    """
    Wrapper for: os.path.islink()
    """
    f = False
    if alink:
        alink = alink[:-1] if alink.endswith(os.sep) else alink
        if os.path.islink(alink):
            f = True
    return f


#############################
class _crc32:
    """
    Wraps up zlib.crc32 to make it suitable for use as a faster but less
    accurate alternative to the hashlib.* classes.
    """
    def __init__(self, initial = None):
        self.crc = 0
        if initial is not None:
            self.update(initial)

    def update(self, block):
        import zlib
        self.crc = zlib.crc32(block, self.crc)

    def hexdigest(self):
        return "%X" % self.crc

    def digest(self):
        # ...
        return self.crc

#############################
class _adler32:
    """
    Wraps up zlib.adler32 to make it suitable for use as a faster but less
    accurate alternative to the hashlib.* classes.
    """
    def __init__(self, initial = None):
        self.crc = 0
        if initial is not None:
            self.update(initial)

    def update(self, block):
        import zlib
        self.crc = zlib.adler32(block, self.crc)

    def hexdigest(self):
        return "%X" % self.crc

    def digest(self):
        # ...
        return self.crc


#############################
class pipeline:
    """
    Pipeline class
    """
    __version__ = '0.98.3 beta'
    __author__  = 'Daniel Nicorici'
    __copyright__ = "Copyright (c) 2009-2020 Daniel Nicorici"
    __credits__ = ["Henrikki Almusa"]
    __license__ = "GPL v3.0"
    __maintainer__ = "Daniel Nicorici"
    __email__ = "daniel.nicorici@gmail.com"
    __status__ = "Beta"


    blank = ' ' # the character used for concatenating the command line arguments for running the job
    count = 0   # counts of pipelines
    __escape = set([" ","\t","(",")"]) # escape characters for which quotes are added to the paths
    __devs = set(['/dev/null','/dev/stdin','/dev/stdout','/dev/stderr','-'])
    #
    #
    #
    def __init__(self,
                 log_filename = 'log_pipeline.txt',
                 checksums_filename = 'checksums.txt',
                 hash_library = 'crc32',
                 threads = 1, # not in use yet
                 start_step = 1 # the number of the starting step (in case that one wants to execute again some specific part of the workflow
                 ):
        """
        Initialization.

        log_filename       - name of the file where the log of running the pipeline is written.
        checksums_filename - name of the file where the checksums for the given pipeline are saved.
        hash_library       - type of hash library used for computing the checksums, e.g. sha512,
                             sha256, md5, crc32, adler32. If it is set to '' or 'no' then no
                             checksums are used and everything is executed.
        threads            - number of threads to be used for runnind the job.
        start_step         - the count of the step from where the execution of workflow should
                             start, default is 1 (in case that one wants to execute again some
                             specific part of the workflow). If it is set to 0 then workflow
                             is restarted automatically from the last step where it was previously
                             stopped.
        """

        self.task = []
        self.workflow = []
        self.log_filename = log_filename
        self.ifs_steps = dict()
        self.ifs_ids = dict()
        self.iffs = set()
        self.threads = 1
        self.task_count = 0
        self.checksums_filename = checksums_filename
        self.start_time = datetime.datetime.now()
        self.hash_library = hash_library
        self.exit_flag = True # True for "no error" reported; False for error found!
        self.closed = False
        self.temp_paths = set()
        self.start_step = start_step
        self.header = dict()
        self.screen_length = 80
        self.protected_paths = set()

        self.hash_type = 'smart' # it can be 'smart' or 'all'
        # - 'smart' - if a file has had its checksum computed before it will
        #             not be computed again. This will work only and only if
        #             the files are not changed "in between" the workflow
        #             tasks/steps using Python code
        # - 'all'   - checksum of all files are computed even if that the given
        #             files have had their checksum computed before. This works
        #             even if the files are changed "in between" the workflow
        #             tasks/steps using some Python code
        self.hash_files = dict() # in case the self.hash_type == 'smart'
                                 # then here are saved the files and their checksums


        if os.path.isfile(log_filename):
            previous_log = [line.rstrip('\r\n') for line in file(log_filename,'r').readlines()]
            if previous_log:
                # build database of IFs' results
                k = "running: step = "
                q = "result of if ="
                qt = "result of if = true"
                qf = "result of if = false"
                t = [line.lower().strip() for line in previous_log if line.lower().find(k) != -1 or line.lower().find(q) != -1]
                self.ifs_steps = dict()
                self.ifs_ids = dict()
                for i in xrange(1,len(t)):
                    if t[i].find(q) == -1:
                        continue
                    s1 = int(t[i-1].partition(k)[2].partition(' ')[0])
                    #
                    mark = t[i].partition('[unique id = ')[2].partition(']')[0].upper()
                    if mark in self.ifs_ids:
                        print >>sys.stderr,"ERROR: the IDs for IFs should be unique! '%s' it is not unique id!" % (mark,)
                        sys.exit(1)
                    self.ifs_ids[mark] = s1
                    #
                    s2 = ''
                    if t[i].find(qt) != -1:
                        s2 = 'yes'
                    elif t[i].find(qf) != -1:
                        s2 = 'no'
                    else:
                        continue
                    self.ifs_steps[s1] = s2

                # automatic restart
                if start_step == 0:
                    self.write("Trying to restart automatically...",log = False)
                    # find the last step where it stopped
                    k = 'running: step = '
                    q = 'execution time:'
                    t = [line for line in previous_log if line.lower().find(k) != -1 or line.lower().find(q) != -1]
                    self.start_step = 1
                    if t and len(t) > 1:
                        if t[-1].lower().find(q) != -1:
                            t = t[-2].lower().rstrip('\r\n').strip() # get the line before the last one
                            t = int(t.split(k,1)[1].split(' ')[0])
                            t = t + 1
                        else:
                            t = t[-1].lower().rstrip('\r\n').strip() # get the last line
                            t = int(t.split(k,1)[1].split(' ')[0])
                        self.start_step = int(t)


            new_log_filename = log_filename
            for i in xrange(1,1000):
                t = "%.3f" % (float(i)/float(1000))
                t = new_log_filename+'.'+t[2:]+'.bak'
                if not os.path.isfile(t):
                    new_log_filename = t
                    break
            if new_log_filename != log_filename:
                shutil.move(log_filename, new_log_filename)
        elif start_step == 0:
            print >> sys.stderr, "WARNING: Cannot restart automatically because the previous log file '%s' cannot be found!" % (log_filename,)
            print >> sys.stderr, "The workflow will be restarted from the beginning with step 1!"
            self.start_step = 1
        temp = ["Log of the pipeline:\n"+
                "-"*self.screen_length]
        self.write(temp)
        self.write("Starting execution with step %d." % (self.start_step,))

    #
    #
    #
    def new(self):
        """
        It prepares the ground for a new task. It is not necessary to be used.
        """
        self.task = []

    #
    #
    #
    def add(self,
            identifier = '',
            value = '',
            kind = 'parameter',
            space = 'yes',
            io = 'input',
            from_file = 'no',
            command_line = 'yes',
            checksum = 'yes',
            dest = None,
            dest_list = None,
            temp_path = 'no'):
        """
        IDENTIFIER   - identifier for command line argument, e.g. '--time'.
                       Always it should be specified, e.g. ''.
        VALUE        - value given to the IDENTIFIER, e.g. '10s'
                       Always it should be specified, e.g. ''.
        KIND         - type of the VALUE, which can be only: 'program', 'parameter',
                       'path', 'input', 'output'. Here 'input'
                       is equivalent to (KIND='path' and IO='input'), 'output'
                       is equivalent to (KIND='path' and IO='output'). Default is 'parameter'.
                       If it is set to 'program' the parameters IO, COMMAND_LINE,
                       and CHECKSUM do not need to be specified.
                       If it is set to 'parameter' the parameter IO does not need
                       to be specified.
                       If it is set to 'path' the VALUE should contain the name of
                       a file or a directory. If it is a directory then the path
                       should end with '*' or with a '/'.
        SPACE        - join the parameters IDENTIFIER and VALUE using a blank space
                       if it is set to 'yes', otherwise is set to 'no'. Default
                       is 'yes'.
        IO           - Input/output type of the VALUE when KIND='path'. It can be
                       only: 'input' or 'output'. Default is 'input'.
        FROM_FILE    - It specifies if the VALUE is passed directly (i.e. 'no')
                       or if it should be read from the file (i.e. 'yes' and the
                       VALUE contains the file name; only the first line of the
                       file is passed. It is useful when a step needs to be skipped
                       but there it is not known its VALUE. Default is 'no'.
        COMMAND_LINE - specifies if the parameters IDENTIFIER and VALUE are passed
                       to the command line for execution. It can be only 'yes' or
                       'no'. Default is 'yes'.
                       If it is set to 'no' the checksum will still be computed for
                       parameters IDENTIFIER and VALUE.
        CHECKSUM     - specifies if the checksum is computed for the parameters
                       IDENTIFIER and VALUE. It can be only 'yes' or 'no'.
                       Default is 'yes'.
        DEST         - destination variable where the VALUE is stored for later use,
                       for example: job.DEST.
        DEST_LIST    - destination variable which is a list and where the VALUE is
                       appended for later use, for example: job.DEST_LIST.
        TEMP_PATH    - in case that KIND='path' (e.g. file or directory), it specifies
                       if the file or directory specified in VALUE should be deleted
                       as soon as possible. It can be only 'yes' or 'no'.
                       Default is 'no', i.e. it is not a temporary path or file
                       and nothing is deleted.

        """
        identifier = str(identifier)
        if not (type(identifier).__name__ == 'str' and
                   (type(value).__name__ == 'str' or
                    type(value).__name__ == 'list' or
                    type(value).__name__ == 'int' or
                    type(value).__name__ == 'float') and
                type(kind).__name__ == 'str' and
                type(space).__name__ == 'str' and
                type(io).__name__ == 'str' and
                type(from_file).__name__ == 'str' and
                type(command_line).__name__ == 'str' and
                type(checksum).__name__ == 'str' and
                type(temp_path).__name__ == 'str'
                ):
            print >> sys.stderr, "ERROR: One of these: IDENTIFIER, VALUE, KIND, SPACE, IO, COMMAND_LINE, CHECKSUM, TEMP_PATH is not of type STR (or LIST for VALUE)!"
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)
        if type(value).__name__ == 'list':
            value = map(str, value)
        elif type(value).__name__ == 'int' or type(value).__name__ == 'float':
            value = str(value)

        kind = str(kind).lower()
        if kind not in ('program','parameter','path','input','output'):
            print >> sys.stderr, "ERROR: ",kind, "is an invalid choice for 'kind'!"
            print >> sys.stderr, "The valid choices for 'kind' are: 'program','parameter', or 'path'."
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)
        elif kind == 'input':
            kind = 'path'
            io = 'input'
        elif kind == 'output':
            kind = 'path'
            io = 'output'
        if kind == 'path' and not value:
            print >> sys.stderr, "ERROR: VALUE cannot be empty in JOB.ADD(...)!"
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        if (kind == 'path' and
            io == 'output' and
            temp_path == 'yes' and
            value and
            (not value.endswith('/')) and
            (not value.endswith('/*')) and
            (not value.endswith('\\*')) and
            (not value.endswith('\*')) and
            (not value.endswith('\\')) and
            (not value.endswith('/')) and
            (not os.path.isdir(value)) and
            (value not in self.__devs) and
            ((not self.hash_library) or self.hash_library == 'no') and
            from_file == 'no'
            ):
            kind = 'parameter'
            io = 'input'
            value = '/dev/null'
            temp_path = 'no'


        space = str(space).lower()
        if space not in ('yes','no'):
            print >> sys.stderr, "ERROR: ",space, "is an invalid choice for 'space'!"
            print >> sys.stderr, "The valid choices for 'space' are: 'yes' or 'no'."
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        io = str(io).lower()
        if io not in ('input','output'):
            print >> sys.stderr, "ERROR: ",io, "is an invalid choice for 'io'!"
            print >> sys.stderr, "The valid choices for 'io' are: 'input' or 'output'."
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        command_line = str(command_line).lower()
        if command_line not in ('yes','no'):
            print >> sys.stderr, "ERROR: ",command_line, "is an invalid choice for 'command_line'!"
            print >> sys.stderr, "The valid choices for 'command_line' are: 'yes' or 'no'."
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        checksum = str(checksum).lower()
        if checksum not in ('yes','no'):
            print >> sys.stderr, "ERROR: ",checksum, "is an invalid choice for 'checksum'!"
            print >> sys.stderr, "The valid choices for 'checksum' are: 'yes' or 'no'."
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        from_file = str(from_file).lower()
        if from_file not in ('yes','no'):
            print >> sys.stderr, "ERROR: ",checksum, "is an invalid choice for 'checksum'!"
            print >> sys.stderr, "The valid choices for 'from_file' are: 'yes' or 'no'."
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        if (not identifier) and (not value) and (kind != 'program'):
            print >> sys.stderr, "ERROR: IDENTIFIER and VALUE cannot be empty in JOB.ADD(...)!"
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        temp_path = str(temp_path).lower()
        if temp_path == 'yes' and kind != 'path' and from_file == 'no':
            print >> sys.stderr, "ERROR: kind = '",kind,"' and temp_path = 'yes' are not allowed!"
            print >> sys.stderr, "The valid choice is kind = 'path' and temp_path = 'yes'!"
            print >> sys.stderr, identifier,value,kind,space,io,from_file,command_line,checksum,temp_path
            self.exit_flag = False
            sys.exit(1)

        if dest is not None:
            setattr(self, dest, value)

        if dest_list is not None:
            if hasattr(self, dest_list):
                gl = getattr(self, dest_list)
                setattr(self, dest_list, gl + [value])
            else:
                setattr(self, dest_list, [value])

        if temp_path == 'yes':
            if from_file == 'yes':
                if kind == 'parameter' and value not in self.__devs:
                    self.temp_paths.add(_expand(value))
                elif kind == 'path':
                    x = [_expand(line.rstrip('\r\n')) for line in file(value,'r')]
                    x.append(_expand(value))
                    self.temp_paths.update(set(x))
                    for e in self.temp_paths:
                        if e in self.__devs:
                            self.temp_paths.remove(e)
            elif kind == 'path' and value not in self.__devs:
                self.temp_paths.add(_expand(value))


        self.task.append({'identifier': identifier,
                          'value': value,
                          'kind': kind,
                          'space': space,
                          'io': io,
                          'from_file': from_file,
                          'command_line': command_line,
                          'checksum': checksum,
                          'temp_path': temp_path})

    #
    #
    #
    def add_list(self,
            identifier = '',
            value = '',
            kind = 'program',
            space = 'yes',
            io = 'input',
            from_file = 'no',
            command_line = 'yes',
            checksum = 'yes',
            list_separator = ' ',
            dest = None,
            dest_list = None,
            temp_path = 'no'):
        """
        It is a wrapper for previous method ADD in case when the VALUE is a list and IDENTIFIER is not a list.
        It is useful if for example one wants to pass a list of files as input or outputs.

        LIST_SEPARATOR - The separator used for concatenating the in a string the list
                         of VALUES. Default is a blank space.
        """
        if type(value).__name__ == 'list' or type(value).__name__ == 'tuple':
            list_separator = str(list_separator)
            value = map(str,value)
            if dest is not None:
                setattr(self,dest,value)
            if dest_list is not None:
                if hasattr(self, dest_list):
                    gl = getattr(self, dest_list)
                    setattr(self, dest_list, gl + [value])
                else:
                    setattr(self, dest_list, [value])
            if kind in ('path','input','output'):
                # add quotes for each file
                v = []
                for val in value:
                    v.append(self.__quote(val,force=True))
                    #v.append(self.__quote(val))
                self.add(identifier = identifier,
                         value = list_separator.join(v),
                         kind = kind,
                         space = space,
                         io = io,
                         from_file = from_file,
                         command_line = command_line,
                         checksum = 'no') # fake
            else:
                # if it is not a path=> no quotes added
                self.add(identifier = identifier,
                         value = list_separator.join(value),
                         kind = kind,
                         space = space,
                         io = io,
                         from_file = from_file,
                         command_line = command_line,
                         checksum = 'no') # fake
            if identifier:
                self.add(identifier = identifier,
                         value = '',
                         kind = 'parameter',
                         space = space,
                         io = 'input',
                         from_file = from_file,
                         command_line = 'no',
                         checksum = checksum) # real - part 1
            for element_value in value:
                self.add(identifier = '',
                         value = element_value,
                         kind = kind,
                         space = space,
                         io = io,
                         from_file = from_file,
                         command_line = 'no',
                         checksum = checksum,
                         temp_path = temp_path) # real - part 2
        else:
            print >> sys.stderr, "ERROR: VALUE is not a list in JOB.ADD_LIST(...)!"
            print >> sys.stderr, identifier,value,kind,space,io,command_line,checksum,temp_path
            self.exit_flag=False
            sys.exit(1)

    ###
    ### RUN
    ###
    def run(self,
            comment = 'no',
            error_message = '',
            successful_exit_status = (0,0),
            exit_code = 0):
        """
        It runs what has been added using method ADD.

        comment       - It can be 'yes' or 'no'. If it is no it means that the entire job
                        has been commented out and it will not be executed. It is a fast
                        way to comment out a job.

        error_message - An additional error message to be displayed if the job/task
                        fails to run.

        It returns:

        True   - if the task has been executed succesfully
        False  - if the task has been skipped from execution

        """
        elk = [1 for elk in self.task if type(elk['value']).__name__ == 'list']
        if not elk: # only if there is only one task to be run

            self.task_count = self.task_count + 1


            # check if there are more or less than one program to run
            x = [1 for element in self.task if element['kind'] == 'program']
            if len(x) != 1:
                self.add(identifier = '',
                         value = '',
                         kind = 'program')
                # if not then added it!
                #print >> sys.stderr, "ERROR: It should be only at least one 'identifier' of type 'program'. If there are more then one program then the rest can be passed as 'parameter'."
                #self.exit_flag = False
                #sys.exit(1)
            # check for  kind = 'program' and value = ''
            empty_program = False
            x = [1 for element in self.task if (element['kind'] == 'program' and
                                                (not element['value']) and
                                                (not element['identifier']))]
            if x:
                empty_program = True

            # build the command line and get the program
            captured_error_message = [] # try to find the '2>&1', '2>', '&>'
            hit_redirect = False
            cmd_line = []
            for element in self.task:
                if element['kind'] == 'program':
                    cmd_line.insert(0,self.__quote(element['identifier']))
                    continue
                if element['command_line'] == 'no':
                    continue
                temp = ''
                element_value = ''
                if element['from_file'] == 'yes': # "from_file"
                    if comment == 'no' and self.task_count >= self.start_step:
                        # read the first line from the text file and remove the new line character
                        element_value = file(element['value'],'rt').readline().rstrip('\r\n')
                    else:
                        element_value = "[ARGUMENT from file '%s']" % (element['value'],)
                else:
                    element_value = self.__quote(element['value']) if element['kind'] in ('path','input','output') else element['value']
                if element['space'] == 'no':
                    temp = element['identifier'] + element_value
                else:
                    temp = element['identifier'] + ' ' + element_value
                if element['value'] and element['value'] not in self.__devs:
                    if element['identifier'] == '2>&1 >' or element['identifier'] == '2>&1>' or element['identifier'] == '2>'  or element['identifier'] == '&>':
                        captured_error_message.append(element['value'])
                        hit_redirect = True
                    if element['identifier'] == '>' or element['identifier'] == 'tee' or element['identifier'] == '| tee':
                        captured_error_message.append(element['value'])
                if element['identifier'] == '2>&1':
                    hit_redirect = True
                cmd_line.append(temp)


            # print the program and command line arguments
            self.__show_step_header_start()

            temp = ' \\\n'.join(cmd_line)
            self.write(temp)
            temp = "-" * self.screen_length
            self.write(temp)
            executed = True

            #execute the program with the given command line arguments
            if comment == 'no': # it is not commented out
                if self.task_count >= self.start_step:
                    if self.__run_again():

                        # EXECUTE IT!
                        proc = 0
                        if not empty_program:
                           #proc=subprocess.call(cmd_line,shell=True)
                            self.write('+-->EXECUTING...')
                            proc = os.system(' '.join(cmd_line))
                        else:
                            self.write('+-->MOCK EXECUTION (i.e. code executed outside of workflow)...')
                        #print "---------------------->",proc,max(successful_exit_status),min(successful_exit_status)
                        exit_code = float(proc)/float(256)
                        if exit_code > max(successful_exit_status) or exit_code < min(successful_exit_status):
                            temp = "\n\nERROR: Workflow execution failed at step %d while executing:\n----------------\n   %s\n----------------\n" % (self.task_count,' \\\n   '.join(cmd_line),)
                            self.write(temp, stderr = True)
                            if error_message:
                                self.write(error_message, stderr = True)
                            # print the input and output file sizes
                            for elem in self.task:
                                if element['command_line'] == 'yes' and ((elem['kind'] == 'path' and elem['io'] in ('input','output')) or elem['from_file'] == 'yes'):
                                    ap = elem['value']
                                    temp = "  * Size '%s' = %d bytes" % (ap,self.__path_size(ap))
                                    self.write(temp, stderr = True)

                            # print the captured error message from '2>&1>' or '2>'
                            if hit_redirect:
                                for il in captured_error_message:
                                    if os.path.isfile(il) or _islink(il):
                                        temp = []
                                        try:
                                            temp = file(il,'r').readlines()
                                        except:
                                            pass
                                        self.write(temp, stderr = True)
                            else:
                                #pass
                                temp = "\n\nExecuting second time the same step/command in order to capture error messages (i.e. STDERR)...\n\n-------------------------------------------"
                                self.write(temp, stderr = True)
                                # no redirection was found so then try to execute again the command and capture the STDERR
                                #p = subprocess.Popen(cmd_line, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = False)
                                #if p.returncode != 0:
                                #    self.write(p.communicate()[0].splitlines(), stderr = True)
                                temp = []
                                temp_file = self.__give_me_temp_filename()
                                procx = os.system(' '.join(cmd_line+['2>',temp_file]))
                                newprocx = float(procx)/float(256)
                                if newprocx > max(successful_exit_status) or newprocx < min(successful_exit_status):
                                    if os.path.isfile(temp_file):
                                        # read error message
                                        try:
                                            temp = file(temp_file,'r').readlines()
                                        except:
                                            pass
                                        if temp:
                                            temp.append("")
                                            temp.append("")
                                        self.write(temp, stderr = True)
                                        os.remove(temp_file)
                                else:
                                    temp = "\n\nWARNING: First execution ended with error but second execution did not! Therefore cannot capture the STDERR!\n\n"
                                    self.write(temp, stderr = True)
                            self.exit_flag = False
                            sys.exit(1)
                        elif self.hash_library and self.hash_library != 'no': # DON'T EXECUTE IT
                            temp = "==> Saving checksum..."
                            self.write(temp)
                            self.__save_checksum()
                    else:
                        temp="|==> SKIPPED because it has not changed since last run."
                        self.write(temp)
                        executed = False
                else:
                    temp="|==> SKIPPED because of the step count."
                    self.write(temp)
                    executed = False
            elif comment == 'yes':
                temp="|==> SKIPPED because it has been commented out."
                self.write(temp)
                executed = False

            # erase the 'temp_path'
            if self.task_count >= self.start_step:
                for task in self.task:
                    if task['temp_path'] == 'yes':
                        v = task['value']
                        if task['from_file'] == 'yes':
                            if task['kind'] == 'parameter':
                                self.__delete_path(v)
                            elif task['kind'] == 'path':
                                x = [line.rstrip('\r\n') for line in file(v,'r')]
                                x.append(v)
                                self.__delete_path(x)
                        elif task['kind'] == 'path':
                            self.__delete_path(v)

            # time difference
            self.__show_step_header_end()
        else:
            # hmm ... there is a list of tasks to be run
            self.write("ERROR: Not implemented this yet!", stderr = True)
            sys.exit(1)
        # cleaning
        #temp = '-'*self.screen_length
        #self.write(temp)
        self.task = []
        self.__class__.count = self.__class__.count + 1
        self.workflow.append(self.task)

        return executed # return if the task has been executed or skipped

    ###
    ###  __RUN_AGAIN
    ###
    def __run_again(self):
        """
        It tests if the the task should be run again or not based on checking if the output and input files have changed since last run.
        """
        flag_run = True
        if self.hash_library and self.hash_library != 'no':
            checksum_now = self.__compute_checksum()
            checksums_old = set()
            if os.path.isfile(self.checksums_filename):
                checksums_old = set([line.strip() for line in file(self.checksums_filename,'rt').readlines() if line.strip()])
            count_outputs = len([0 for elem in self.task if elem['io'] == 'output'])
            if (checksum_now in checksums_old) and count_outputs > 0:
                flag_run = False
        if flag_run:
            if (self.hash_library and self.hash_library != 'no'):
                self.__build_paths(directories = True, files = True)
            else:
                self.__build_paths(directories = True, files = False)
        return flag_run

    ###
    ### __COMPUTE_CHECKSUM
    ###
    def __compute_checksum(self):
        """
        It computes the checksum of the input and output files and command line string.
        """
        cmd_line = ''
        list_files_to_check = []
        output_files = set()
        for elem in self.task:
            ident = elem['identifier']
            value = elem['value']
            kind = elem['kind']
            io = elem['io']
            checksum = elem['checksum']
            if checksum == 'no':
                continue
            # building list of files for the checksum
            if kind == 'path':
                deposit = [value]
                for a_value in deposit:
                    if a_value.endswith('\\*') or a_value.endswith('/*'):
                        a_value = a_value[:-1]
                    if os.path.exists(a_value):
                        if os.path.isfile(a_value) or _islink(a_value):
                            list_files_to_check.append(a_value[:])
                            if io == 'output':
                                output_files.add(a_value[:])
                        elif os.path.isdir(a_value):
                            temp = []
                            for el in sorted(os.listdir(a_value)):
                                el = os.path.join(a_value,el)
                                if os.path.isfile(el):
                                    temp.append(el)
                            temp.sort() # sort the list of files
                            list_files_to_check.extend(temp)
                            if io == 'output':
                                for el in temp:
                                    output_files.add(el)
                        else:
                            print >>sys.stderr, "\nERROR: '%s' is not a file nor a directory." % (value,)
                            self.exit_flag = False
                            sys.exit(1)
            # bulding the command line for hash
            if kind in ('path','program'):
                value = os.path.basename(value) # get the filename
            cmd_line = cmd_line+'___'.join([ident,
                                          value,
                                          kind,
                                          io
                                          ])
        # compute the checksum
        checksum = ''

        fingerprint = self.__hashlib_init()
        fingerprint.update(cmd_line)
        checksum = checksum + fingerprint.hexdigest()

        length_piece = 2**26 # it is ~262144
        for a_file in list_files_to_check:
            if a_file in self.__devs:
                continue
            a = _expand(a_file)
            b = _expand(os.readlink(a)) if _islink(a) else None
            dig = None
            if (self.hash_type == 'smart' and
                (a not in output_files) and
                (b not in output_files) and
                ((a in self.hash_files) or (b in self.hash_files))):
                temp = "===> Checksum has been computed previously for: '%s'" % (a_file,)
                self.write(temp)
                dig = self.hash_files[a]
                if b:
                    self.hash_files[b] = self.hash_files[a]
            else:
                temp = "===> Computing checksum for: '%s'" % (a_file,)
                self.write(temp)
                file_fingerprint = self.__hashlib_init()
                ff = open(a_file,'rb')
                while True:
                    dd = ff.read(length_piece)  # originally was dd=ff.read(8096)
                    if not dd:
                        break
                    file_fingerprint.update(dd)
                ff.close()
                dig = file_fingerprint.hexdigest()
                self.hash_files[a] = dig
                if b:
                    self.hash_files[b] = dig
            checksum = checksum + dig


        return checksum # checksum

    ###
    ### __hashlib
    ###
    def __hashlib_init(self):
        """
        Gives the finger print to the corresponding hash library.
        """
        finger = None
        if self.hash_library.lower() == 'sha512':
            finger = hashlib.sha512()
        elif self.hash_library.lower() == 'md5':
            finger = hashlib.md5()
        elif self.hash_library.lower() == 'sha256':
            finger = hashlib.sha256()
        elif self.hash_library.lower() == 'crc32':
            finger = _crc32()
        elif self.hash_library.lower() == 'adler32':
            finger = _adler32()
        else:
            print >> sys.stderr,"ERROR: Unknown hash library!"
            self.exit_flag = False
            sys.exit(1)
        return finger

    ###
    ###  __BUILD_PATHS
    ###
    def __build_paths(self,
                      directories = True,
                      files = True,
                      links = True):
        """
        It creates empty output files and the output directories when they do not exist

        directories - If True then directories are created else they are not.
        files       - If True then files are created else they are not.
        links       - If True then links are deleted.
        """
        for elem in self.task:
            ident = elem['identifier']
            value = elem['value']
            kind = elem['kind']
            io = elem['io']
            checksum = elem['checksum']
            #
            if (io == 'output' and
                kind == 'path' and
                checksum == 'yes'):
                deposit = [value]
                for a_value in deposit:
                    if a_value in self.__devs:
                        continue
                    if a_value.endswith('*'):
                        a_value = a_value[:-1]
                    if not os.path.exists(a_value):
                        if (a_value.endswith('\\') or
                            a_value.endswith('/') or
                            os.path.isdir(a_value)):
                            # it is a directory
                            # NOTE: Assumption is that the directory ends with a path separator and a file does not!!!
                            if directories:
                                if (a_value.endswith('\\') or a_value.endswith('/')) and os.path.exists(a_value[:-1]):
                                    # test the there is no file with exactly the same name
                                    self.__delete_path(a_value[:-1])
                                os.makedirs(a_value)
                        else: # it is a file
                            (dir_name,file_name) = os.path.split(a_value)
                            if (dir_name and
                                (not os.path.exists(dir_name)) and
                                directories):

                                os.makedirs(dir_name)
                            if links and _islink(a_value):
                                self.__delete_path(a_value)
                            if files:
                                file(a_value,'w').write('')

    ###
    ### __SAVE_CHECKSUM
    ###
    def __save_checksum(self):
        """
        It saves the current checksum.
        """
        checksum = self.__compute_checksum()
        checksums_old = set([])
        if os.path.isfile(self.checksums_filename):
            checksums_old = set([line.rstrip("\r\n") for line in file(self.checksums_filename,'rt').readlines() if line.rstrip("\r\n")])
        if not (checksum in checksums_old):
            no_tries = 10
            for ix in xrange(no_tries + 1): # try 10 time if it cannot write
                try:
                    file(self.checksums_filename,'a').write(checksum+'\n')
                except IOError:
                    print >> sys.stderr,"==> Warning: Cannot write to file: ",self.checksums_filename
                    print >> sys.stderr,"==> Waiting 5 seconds and trying one more time..."
                    if no_tries == ix:
                        print >> sys.stderr,"==> ERROR: Cannot write to file: ",self.checksums_filename
                        self.exit_flag = False
                        sys.exit(1)
                    else:
                        time.sleep(5)
                else:
                    break


    ###
    ### __SHOW_STEP_HEADER
    ###
    def __show_step_header_start(self, id = None):
        if not id:
            id = self.task_count
        time_now = datetime.datetime.today()
        time_start = datetime.datetime.now()
        self.header[id] = time_start
        str_time_now = time_now.strftime("%H:%M")
        str_date_now = time_now.strftime("%Y-%m-%d")

        time_difference = time_start - self.start_time
        minutes = math.ceil(float(time_difference.seconds) / float(60))
        hours, minutes = divmod(minutes, 60)

        temp = ["/"*self.screen_length,
                "  Running: step = %d   Time: %s   Date: %s (elapsed time: %dd:%dh:%dm)\n" % (self.task_count,str_time_now,str_date_now,time_difference.days,hours,minutes),
                "\\"*self.screen_length,
                "==> Current working directory: '%s'\n" % (os.getcwd(),)]
        self.write(temp)



    ###
    ### __SHOW_STEP_HEADER
    ###
    def __show_step_header_end(self, id = None):
        if not id:
            id = self.task_count
        time_stop = datetime.datetime.now()
        time_start = self.header[id]
        time_difference = time_stop - time_start
        minutes, seconds = divmod(time_difference.seconds, 60)
        hours, minutes = divmod(minutes, 60)
        temp = ['-'*self.screen_length,
                "==> Execution time: %d day(s), %d hour(s), %d minute(s), and %d second(s)" % (time_difference.days,hours,minutes,seconds)]
        self.write(temp)


    ###
    ### remove files and directories
    ###
    def __delete_path(self, paths):
        if type(paths).__name__ == 'set':
            paths = list(paths)
        if type(paths).__name__ != 'list':
            paths = [paths]
        paths = list(set(paths))
        for a_path in paths:
            if not a_path:
                continue
            if a_path.endswith('\\*') or a_path.endswith('/*'):
                a_path = a_path[:-1]
            if self.__isprotected(a_path):
                temp = "==> Protected path (will not be erased/moved): '%s' (size: %d bytes)" % (a_path,self.__path_size(a_path))
                self.write(temp)
                if a_path in self.temp_paths:
                    self.temp_paths.remove(a_path)
                continue
            if os.path.exists(a_path) or _islink(a_path):
                if os.path.isdir(a_path):
                    temp = "==> Erasing directory (and all directories and files within): '%s' (size: %d bytes)" %(a_path,self.__path_size(a_path))
                    self.write(temp)
                    shutil.rmtree(a_path)
                    if a_path in self.temp_paths:
                        self.temp_paths.remove(a_path)
                elif os.path.isfile(a_path) or _islink(a_path):
                    temp = "==> Erasing file: '%s' (size: %d bytes)"%(a_path,self.__path_size(a_path))
                    self.write(temp)
                    os.remove(a_path)
                    if a_path in self.temp_paths:
                        self.temp_paths.remove(a_path)
                elif a_path in self.__devs:
                    temp = "==> Pretending to erase the file: '%s' (size: %d bytes)"%(a_path,self.__path_size(a_path))
                    self.write(temp)
                    if a_path in self.temp_paths:
                        self.temp_paths.remove(a_path)
                else:
                    temp = "'%s' is not a directory, not file, and not a link. It will not be deleted!" % (a_path,)
                    self.write(temp)


    ###
    ### __DEL__
    ###
    def close(self):
        """
        Close the pipeline.
        """
        if not self.closed:
            # delete the paths and files marked as temporary
            if self.exit_flag: # if there is no error
                self.__delete_path(self.temp_paths)

            # time counting
            end_time = datetime.datetime.now()
            time_difference = end_time - self.start_time
            minutes, seconds = divmod(time_difference.seconds, 60)
            hours, minutes = divmod(minutes, 60)
            temp = ["#"*self.screen_length,
                    "#"*self.screen_length,
                    "TOTAL RUNNING TIME: %d day(s), %d hour(s), %d minute(s), and %d second(s) \n" % (time_difference.days,hours,minutes,seconds),
                    "#"*self.screen_length,
                    "#"*self.screen_length,
                    ]
            self.write(temp)
            self.closed = True



    ###
    ### __DEL__
    ###
    def __del__(self):
        """
        Destructor
        """
        if not self.closed: # if it has not been closed yet
            self.close()


    ###
    ### PRINT
    ###
    def write(self,
              text = '',
              stdout = True,
              stderr = False,
              log = True):
        """
        Prints/writes the text on screen and adds it also the to log file.

        text   - the text to print
        stdout - the text is printed to the stdout if True. If True then stderr is excluded.
        stderr - the text is printed to the stderr if True. If True then stderr is excluded.
        log    - the text is added also to the log file.

        """
        if type(text).__name__ == 'str':
            text = [text]

        t = [e.rstrip('\r\n') for e in text]
        if stdout:
            for e in t:
                print >>sys.stdout, e
        elif stderr:
            for e in t:
                print >> sys.stderr, e
        if log:
            t = [e+'\n' for e in t]
            file(self.log_filename,'a').writelines(t)


    ###
    ### LINK
    ###
    def link(self,
             fin,
             fout,
             temp_path = 'no',
             kind = 'soft',
             dest = None,
             dest_list = None
             ):
        """
        It creates a 1) soft link, 2) hard link, or 3) it copies the input file
        or path to the output. It also deletes the output file/path before
        the process.

        fin        - input path/file name
        fout       - output path/file name
        temp_path  - it specifies if the input file or directory should be
                     deleted as soon as possible. It can be only 'yes' or 'no'.
                     Default is 'no', i.e. it is not a temporary path or file
                     and nothing is deleted.
        kind       - it can take the values ['soft', 'hard', 'copy','move'] where:
                     'soft' will create a soft link from the input to the output,
                     'hard' will create a hard link from the input to the output,
                     'copy' will copy the input to the output,
                     'move' will move the input to the output.
        dest         - destination variable where the VALUE is stored for later use,
                       for example: job.dest.
        dest_list    - destination variable which is a list and where the FOUT is
                       appended for later use, for example: job.dest_list.


        It returns:

        True   - if the task has been executed succesfully
        False  - if the task has been skipped from execution

       """
        self.task_count = self.task_count + 1
        self.__show_step_header_start()
        executed = True

        if dest is not None:
            setattr(self, dest, fout)

        if dest_list is not None:
            if hasattr(self, dest_list):
                gl = getattr(self, dest_list)
                setattr(self, dest_list, gl + [fout])
            else:
                setattr(self, dest_list, [fout])

        if self.task_count >= self.start_step:
            if os.path.exists(fout) or _islink(fout):
                self.__delete_path(fout)
            linkfrom = fin
            if _islink(fin):
                linkfrom = _expand(os.readlink(fin))
            if (not self.__isprotected(fin)) and ((temp_path == 'yes' and kind == 'soft') or kind == 'move'):
                try:
                    shutil.move(fin, fout)
                except OSError as er:
                    self.write("ERROR: Cannot move ('%s' and '%s')!  [Error number: %s][%s]" % (fin,fout,str(er.errno),str(er)))
                    sys.exit(1)
                self.write("Moving from:\n'%s'\nto:\n'%s'" % (fin,fout) )
            elif kind == 'copy':
                if os.path.isdir(linkfrom):
                    try:
                        shutil.copytree(linkfrom,fout)
                    except OSError as er:
                        self.write("ERROR: Cannot copy ('%s' and '%s')!  [Error number: %s][%s]" % (linkfrom,fout,str(er.errno),str(er)))
                        sys.exit(1)
                else:
                    try:
                        shutil.copyfile(linkfrom,fout)
                    except OSError as er:
                        self.write("ERROR: Cannot copy ('%s' and '%s')!  [Error number: %s][%s]" % (linkfrom,fout,str(er.errno),str(er)))
                        sys.exit(1)
                self.write("Copying from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
            elif kind == 'soft':
                linkfrom = _expand(linkfrom)
                fout = _expand(fout)
                fox = False
                try:
                    os.symlink(linkfrom,fout)
                except OSError as er:
                    self.write("WARNING: Cannot do symbolic links ('%s' and '%s')! Copying will be used instead of symbolic links and this will be very very slow! Please, fix this by testing that it is possible to create symbolic links on this drive/path. [Error number: %s][%s]" % (linkfrom,fout,str(er.errno),str(er)))
                    shutil.copyfile(linkfrom,fout)
                    fox = True
                if fox:
                    self.write("Copying from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
                else:
                    self.write("Soft linking from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
            elif kind == 'hard':
                if os.path.isdir(linkfrom):
                    try:
                        shutil.copytree(linkfrom,fout)
                    except OSError as er:
                        self.write("ERROR: Cannot copy ('%s' and '%s')!  [Error number: %s][%s]" % (linkfrom,fout,str(er.errno),str(er)))
                        sys.exit(1)
                    self.write("Copying from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
                else:
                    linkfrom = _expand(linkfrom)
                    fout = _expand(fout)
                    try:
                        os.link(linkfrom,fout)
                        self.write("Hard linking from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
                    except OSError as er:
                        self.write("WARNING: Cannot do hard links ('%s' and '%s')! [Error number: %s][%s]" % (linkfrom,fout,str(er.errno),str(er)))
                        shutil.copyfile(linkfrom,fout)
#                        if er.errno == errno.EXDEV:
#                            # they are on different partitions
#                            # [Errno 18] Invalid cross-device link
#                            shutil.copyfile(linkfrom,fout)
#                            self.write("Copying from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
#                        else:
#                            print >>sys.stderr,"ERROR: Cannot do hard links ('%s' and '%s')!" % (linkfrom,fout)
#                            print >>sys.stderr,er
#                            sys.exit(1)

        else:
            linkfrom = fin
            if os.path.exists(fin) and _islink(fin):
                linkfrom = _expand(os.readlink(fin))
            if (not self.__isprotected(fin)) and ((temp_path == 'yes' and kind == 'soft') or kind == 'move'):
                self.write("Skipping the moving from:\n'%s'\nto:\n'%s'" % (fin,fout) )
            elif kind == 'copy':
                self.write("Skipping the copying from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
            elif kind == 'soft':
                if os.path.exists(linkfrom):
                    linkfrom = _expand(linkfrom)
                if os.path.exists(fout):
                    fout = _expand(fout)
                self.write("Skipping the soft linking from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
            elif kind == 'hard':
                if os.path.exists(linkfrom):
                    linkfrom = _expand(linkfrom)
                if os.path.exists(fout):
                    fout = _expand(fout)
                self.write("Skipping the hard linking from:\n'%s'\nto:\n'%s'" % (linkfrom,fout) )
            executed = False
        self.__show_step_header_end()
        return executed # return if the task has been executed or skipped

    ###
    ### SINK
    ###
    def sink(self,
             variable,
             fout):
        """
        It writes in an output file the content of variable. It also deletes
        the output file/path before the process.

        variable   - a string/integer/float or a (python) list
        fout       - output path/file name

        It returns:

        True   - if the task has been executed succesfully
        False  - if the task has been skipped from execution

       """
        self.task_count = self.task_count + 1
        self.__show_step_header_start()
        executed = True
        if self.task_count >= self.start_step:
            if os.path.exists(fout) or _islink(fout):
                self.__delete_path(fout)
            if variable:
                data = []
                if type(variable).__name__ == 'list':
                    data = map(str,variable)
                else:
                    data.append(str(variable))
                data = [line+'\n' for line in data]
                file(fout,'w').writelines(data)
            else:
                file(fout,'w').write('')
            self.write("Writing internal variable(s)/parameter(s) into:\n'%s'" % (fout,) )
        else:
            self.write("Skipping from writing internal variable(s)/parameter(s) into:\n'%s'" % (fout,) )
            executed = False
        self.__show_step_header_end()
        return executed # return if the task has been executed or skipped



    ###
    ### IF for WORKFLOW (it remembers what was the last value of it)
    ###
    def iff(self,
           boolean = False,
           id = ''):
        """
        This should be used in Python in IF constructions.
        In case this IF has not been executed in a previous workflow run then
        FALSE is returned. In case this IF has been executed in a previous
        workflow run then the previous IF result is returned (that is TRUE or
        FALSE). This is executed as a mock step in the workflow.

        It returns:

        True   - if the IF task has been executed in a previous workflow run and
                 if the result of the the entire IF was true
        False  - if the task has been skipped from execution

       """
        self.task_count = self.task_count + 1
        self.__show_step_header_start()

        if type(boolean).__name__ != 'bool':
            self.write("ERROR: IF accepts for 'boolean' only True/False values!",stderr=True,stdout=False)
            sys.exit(1)
        if type(id).__name__ != 'str':
            self.write("ERROR: IF accepts for 'id' only strings!",stderr=True,stdout=False)
            sys.exit(1)
        if id.find(']') != -1 or id.find('[') != -1 or id.find(' ') != -1 or id.find('"') != -1 or id.find("'") != -1 or id.find("\\") != -1 or id.find("/'") != -1 or id.find("*") != -1 or id.find(",") != -1 or id.find(";") != -1:
            self.write("ERROR: IF accepts for 'id' only strings which do not contain the following characters '[] \"\''\\/*!,;",stderr=True,stdout=False)
            sys.exit(1)
        if id.upper() in self.iffs:
            self.write("ERROR: IF accepts for 'id' which is unique! '%s' is not unique in this workflow!" % (id.upper(),),stderr=True,stdout=False)
            sys.exit(1)

        self.iffs.add(id.upper())

        result = boolean
        if self.task_count >= self.start_step:
            self.write("+-->EXECUTING workflow's IF...")
            if boolean:
                self.write("Result of IF = TRUE [Unique Id = %s]" % (id.upper(),))
                result = True
            else:
                self.write("Result of IF = FALSE [Unique Id = %s]" % (id.upper(),))
                result = False
        else:
            self.write("|==> SKIPPED workflow's IF because of the step count (value returned in the previous run is used).")
            if self.ifs_steps.has_key(self.task_count):
                r = self.ifs_steps.get(self.task_count)
                k = self.ifs_ids.get(id.upper(),-1)
                if k != self.task_count:
                    self.write("WARNING: The step count and unique id of this IF do not match the ones from the previous run!")
                if r == 'yes':
                    self.write("Result of IF = TRUE [Unique Id = %s]" % (id.upper(),))
                    result = True
                elif r == 'no':
                    self.write("Result of IF = FALSE [Unique Id = %s]" % (id.upper(),))
                    result = False
                else:
                    self.write("Result of IF = cannot be parsed out!!! [Unique Id = %s]" % (id.upper(),))
            else:
                self.write("Result of IF = not found in the previous log run!!! [Unique Id = %s]" % (id.upper(),))
        self.__show_step_header_end()
        return result

    ###
    ### clean
    ###
    def clean(self, a_path, list_file = 'no', temp_path = 'yes'):
        """
        If the path 'a_path' exists and 'temp_path' is set to 'yes' then the
        path 'a_path' will be deleted.

        temp_path = 'yes' or 'no'

        It returns:

        True   - if the task has been executed succesfully
        False  - if the task has been skipped from execution
        """
        self.task_count = self.task_count + 1
        self.__show_step_header_start()

        if type(a_path).__name__ == 'str':
            if temp_path == 'yes':
                self.write("Erasing '%s' (size: %d bytes)" %(a_path, self.__path_size(a_path)))
            elif temp_path == 'no':
                self.write("Skipping from erasing '%s' (size: %d bytes)" % (a_path, self.__path_size(a_path)))
            else:
              print >> sys.stderr,"Error: unknown temp_path value!"
              sys.exit(1)
              
            executed = True
            if self.task_count >= self.start_step:
                if (os.path.exists(a_path) or _islink(a_path)) and temp_path == 'yes':
                    if list_file == 'yes':
                        for ap in file(a_path,'r').readlines():
                            apr = ap.rstrip('\r\n')
                            if apr and (os.path.exists(apr) or _islink(apr)):
                                self.__delete_path(apr)
                    self.__delete_path(a_path)
            else:
                executed = False
        ########################################################################
        elif type(a_path).__name__ == 'list':
            for ap in a_path:
                if temp_path == 'yes':
                    self.write("Erasing '%s' (size: %d bytes)" %(ap, self.__path_size(ap)))
                elif temp_path == 'no':
                    self.write("Skipping from erasing '%s' (size: %d bytes)" % (ap, self.__path_size(ap)))
                else:
                  print >> sys.stderr,"Error: unknown temp_path value!"
                  sys.exit(1)
                  
            executed = True
            if self.task_count >= self.start_step:
                for ap in a_path:
                    if (os.path.exists(ap) or _islink(ap)) and temp_path == 'yes':
                        if list_file == 'yes':
                            for apl in file(ap,'r').readlines():
                                apr = apl.rstrip('\r\n')
                                if apr and (os.path.exists(apr) or _islink(apr)):
                                    self.__delete_path(apr)
                        self.__delete_path(ap)
            else:
                executed = False
        

        self.__show_step_header_end()
        return executed # return if the task has been executed or skipped


    ###
    ### PROTECTS paths from deletion
    ###
    def protect(self, paths = None):
        """
        It adds a path to a list of protected paths which can be erased even
        if they are marked later as temporary paths (i.e. temp_path = 'yes').
        For example, this can be used to protect the input files from deletion.

        paths    - can be one path or a list of paths. If it is empty then it
                   returns the protected paths.
        """
        if paths:
            if type(paths).__name__ != 'list':
                paths = [paths]
            for apath in paths:
                if apath in self.__devs:
                    continue
                a = _expand(apath)
                if a.endswith('\\') or a.endswith('/'):
                    a = a [:-1]
                print "Protecting from future deletion attempts '%s'" %(a,)
                self.protected_paths.add(a)
        else:
            # return the list of protected paths
            return self.protected_paths


    ###
    ### Removes paths from list of PROTECTED paths
    ###
    def unprotect(self, paths = None):
        """
        It removes a path from a list of protected paths which can be erased even
        if they are marked later as temporary paths (i.e. temp_path = 'yes').

        paths    - can be one path or a list of paths. If it is empty then all
                   protected paths are removed from the list protected paths
                   (but not deleted).
        """
        if paths:
            if type(paths).__name__ != 'list':
                paths = [paths]
            for apath in paths:
                if apath in self.__devs:
                    continue
                a = _expand(apath)
                if a.endswith('\\') or a.endswith('/'):
                    a = a [:-1]
                if a in self.protected_paths:
                    self.protected_paths.remove(a)
        else:
            # reset the list of protected paths
            self.protected_paths = set()

    ###
    ### Tests if a given path is in the list of PROTECTED paths
    ###
    def __isprotected(self, path):
        """
        It if a given path is the list of protected paths which can be erased even
        if they are marked later as temporary paths (i.e. temp_path = 'yes').

        path     - one path.
        """
        a = _expand(path)
        if a.endswith('*'):
            a = a[:-1]
        if a.endswith('\\') or a.endswith('/'):
            a = a[:-1]
        r = True
        if a in self.protected_paths:
            r = True
        else:
            r = False
            for p in self.protected_paths:
                if a.startswith(p):
                    r = True
                    break
        return r

    ###
    ### QUOTE
    ###
    def __quote(self, a_path, force = False):
        """
        It adds quotes to a path if it contains blank characters.

        force    - It adds always the quotes if it is True.

        """
        flag = False
        new_path = a_path
        if ((new_path.startswith('"') and new_path.endswith('"') ) or
            (new_path.startswith("'") and new_path.endswith("'") )):
            pass
        else:
            for c in self.__escape:
                if a_path.find(c)!=-1:
                    flag = True
                    break
            if flag or force:
                new_path = '"%s"' % (a_path,)
        return new_path

    ###
    ### QUOTE
    ###
    def __path_size(self, apath = '.'):
        total_size = 0
        if os.path.isfile(apath):
            try:
                total_size = os.path.getsize(apath)
            except:
                total_size = 0
        elif os.path.isdir(apath):
            for dirpath, dirnames, filenames in os.walk(apath):
                for f in filenames:
                    fp = os.path.join(dirpath, f)
                    if os.path.exists(fp) and (not _islink(fp)):
                        t = 0
                        try:
                            t = os.path.getsize(fp)
                        except:
                            continue
                        total_size = total_size + t
        return total_size

    ###
    ### temporary file
    ###
    def __give_me_temp_filename(self, tmp_dir = None):
        if tmp_dir:
            if not (os.path.isdir(tmp_dir) or _islink(tmp_dir)):
                os.makedirs(tmp_dir)
        (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
        os.close(ft)
        return ft_name

###
### MAIN
###
if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.dirname(sys.argv[0])))
    os.system('python -c "import workflow; help(workflow)"')
    print "=================================================================================="
    print "Testing..."
    job = pipeline(hash_library='')
    job.add('dir', kind = 'program')
    job.add('>','list.txt', kind = 'output')
    job.run()
#
