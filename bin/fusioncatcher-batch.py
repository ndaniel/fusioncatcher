#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It is a wrapper for 'fusioncatcher.py' which allows to run 'fusioncatcher.py' in
batch mode where the input is a directory which contains the directories for
which FusionCatcher should be run for each of them separately.

The input file (given using -i or --input) may be:
- a text file with two columns (that are tab separated), where the first column
contains the a list of paths (one line = one path = one sample) such that
a path is a SRA file, or a directory containing Fastq files, and the second
column contains the corresponding sample name which will be used to construct
the output directory (also URLs may be given here)
- a text file which contains the a list of paths (one line = one path = one sample)
such that a path is a SRA file, or a directory containing Fastq files, and the
second column contains the corresponding sample name which will be used to
construct the output directory (also URLs may be given here)
- a directory which contains subdirectories and each subdirectory contains Fastq
file corresponding to one sample (one subdirectory = one sample); all
found subdirectories will be analyzed one by one by FusionCatcher



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
import shutil

if __name__ == "__main__":
    #command line parsing
    cmds = sys.argv
    head = cmds.pop(0)

    # hijack input
    # --input OR -i
    xin = None
    found = [el.replace('--input=','').replace('-i=','') for el in cmds if el.startswith('--input=') or el.startswith('-i=')]
    if not found:
        found = [cmds[i] for i in xrange(1,len(cmds)) if cmds[i-1] == '--input' or cmds[i-1] == '-i']
    if found:
        xin = found[0]

    # hijack output
    # --output OR -o
    xou = None
    found = [el.replace('--output=','').replace('-o=','') for el in cmds if el.startswith('--output=') or el.startswith('-o=')]
    if not found:
        found = [cmds[i] for i in xrange(1,len(cmds)) if cmds[i-1] == '--output' or cmds[i-1] == '-o']
    if found:
        xou = found[0]

    # hijack normal matched
    # --normal
    xno = None
    found = [el.replace('--normal=','') for el in cmds if el.startswith('--normal=')]
    if not found:
        found = [cmds[i] for i in xrange(1,len(cmds)) if cmds[i-1] == '--normal']
    if found:
        xno = found[0]

    flag_reverse = False
    if "--reverse" in set(cmds):
        cmds = [el for el in cmds if el != '--reverse']
        flag_reverse = True

    newcmds = cmds[:]
    if xin and xou:
        # new command line options
        exclusion = ('--output','-o','--input','-i','--normal')
        newcmds = [newcmds[i] for i in xrange(0,len(newcmds)) if (i == 0) or (newcmds[i-1] not in exclusion)]
        newcmds = [el for el in newcmds if ((not el.startswith('-i=')) and
                                            (not el.startswith('-o=')) and
                                            (not el.startswith('--input=')) and
                                            (not el.startswith('--output=')) and
                                            (not el.startswith('--normal=')) and
                                            el!='-i' and
                                            el!='--input' and
                                            el!='--output' and
                                            el!='--normal' and
                                            el!='-o')]
        if head.endswith('fusioncatcher-batch.py'):
            newcmds.insert(0,head.replace('fusioncatcher-batch.py','fusioncatcher.py --keep-preliminary'))
        else:
            newcmds.insert(0,head.replace('fusioncatcher-batch','fusioncatcher.py --keep-preliminary'))
        if xin and os.path.isdir(xin):
            # input is a directory and contains subdirectories, one subdirectory is one sample
            dirs = sorted([el for el in os.listdir(xin) if os.path.isdir(os.path.join(xin,el)) and not el.startswith('.')],reverse=flag_reverse)
            nos = None
            if xno and os.path.isdir(xno):
                nos =  sorted([el for el in os.listdir(xno) if os.path.isdir(os.path.join(xno,el)) and not el.startswith('.')])
            if nos:
                if xou and not os.path.isdir(xou):
                    os.makedirs(xou)
                normaltemp = os.path.join(xou,'preliminary-candidate-fusions-found-in-normal.log')
                file(normaltemp,"w").write('')
                partialnormaltemp = os.path.join(xou,'partial-preliminary-candidate-fusions-found-in-normal.log')
                file(partialnormaltemp,"w").write('')
                first = True
                for d in nos:
                    t = newcmds[:]
                    t.append('--input')
                    t.append(os.path.join(xin,d))
                    t.append('--output')
                    t.append(os.path.join(xou,d))
                    if first:
                        first = False
                    else:
                        t.append('--label-title')
                        t.append('partial-matched-normal,matched-normal')
                        t.append('--label-file')
                        t.append(partialnormaltemp+','+normaltemp)
                        t.append('--label-threshold')
                        t.append('2,0')
                    t = ' '.join(t)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)
                    t = 'sed "1 d" "%s" | cut -f 11,12 | uniq >> "%s"' % (os.path.join(xou,d,'final-list_candidate-fusion-genes.txt'),normaltemp)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)
                    t = 'sed "1 d" "%s" | cut -f 1-3 | uniq >> "%s"' % (os.path.join(xou,d,'preliminary-list_candidate-fusion-genes.txt'),partialnormaltemp)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)

                for d in dirs:
                    t = newcmds[:]
                    t.append('--input')
                    t.append(os.path.join(xin,d))
                    t.append('--output')
                    t.append(os.path.join(xou,d))
                    t.append('--label-title')
                    t.append('partial-matched-normal,matched-normal')
                    t.append('--label-file')
                    t.append(partialnormaltemp+','+normaltemp)
                    t.append('--label-threshold')
                    t.append('2,0')
                    t = ' '.join(t)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)

                os.remove(normaltemp)

            else:
                for d in dirs:
                    t = newcmds[:]
                    t.append('--input')
                    t.append(os.path.join(xin,d))
                    t.append('--output')
                    t.append(os.path.join(xou,d))
                    t = ' '.join(t)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)

        elif xin and os.path.isfile(xin):
            # the input is text file containing the files/directories which should be given as input to FusionCatcher
            # ignore the lines which are empty of start with #
            # a second column (separated by tab) can contain sample name which will be used for name of output directory
            txt = [line.rstrip('\r\n').split('\t') for line in file(xin,'r').readlines() if line.rstrip("\r\n") and not line.strip().startswith("#")]
            nos = None
            if xno:
                nos = [line.rstrip('\r\n').split('\t') for line in file(xno,'r').readlines() if line.rstrip("\r\n") and not line.strip().startswith("#")]

            if nos:
                if xou and not os.path.isdir(xou):
                    os.makedirs(xou)
                normaltemp = os.path.join(xou,'preliminary-candidate-fusions-found-in-normal.log')
                file(normaltemp,"w").write('')
                partialnormaltemp = os.path.join(xou,'partial-preliminary-candidate-fusions-found-in-normal.log')
                file(partialnormaltemp,"w").write('')

                first = True
                ix = 0
                for line in nos:
                    ix = ix + 1
                    t = newcmds[:]
                    pin = None
                    pou = None
                    pin = os.path.join(line[0])
                    if len(line) >= 2 and line[1]:
                        pou = os.path.join(xou,line[1])
                    else:
                        #pou = os.path.join(xou,os.path.basename(line[0].rstrip(os.sep)))
                        vu = line[0]
                        if vu.find(',') != -1:
                            vu = str(ix) + "_" + os.path.basename(vu.split(",")[0].rstrip(os.sep))
                        else:
                            vu = os.path.basename(vu.rstrip(os.sep))
                        pou = os.path.join(xou,vu)
                    t.append('--input')
                    t.append(pin)
                    t.append('--output')
                    t.append(pou)
                    if first:
                        first = False
                    else:
                        t.append('--label-title')
                        t.append('partial-matched-normal,matched-normal')
                        t.append('--label-file')
                        t.append(partialnormaltemp+','+normaltemp)
                        t.append('--label-threshold')
                        t.append('2,0')
                    t = ' '.join(t)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)
                    t = 'sed "1 d" "%s" | cut -f 11,12 | uniq >> "%s"' % (os.path.join(pou,'final-list_candidate-fusion-genes.txt'),normaltemp)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)
                    t = 'sed "1 d" "%s" | cut -f 1-3 | uniq >> "%s"' % (os.path.join(pou,'preliminary-list_candidate-fusion-genes.txt'),partialnormaltemp)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)

                for line in txt:
                    ix = ix + 1
                    t = newcmds[:]
                    pin = None
                    pou = None
                    pin = os.path.join(line[0])
                    if len(line) >= 2 and line[1]:
                        pou = os.path.join(xou,line[1])
                    else:
                        #pou = os.path.join(xou,os.path.basename(line[0].rstrip(os.sep)))
                        vu = line[0]
                        if vu.find(',') != -1:
                            vu = str(ix) + "_" + os.path.basename(vu.split(",")[0].rstrip(os.sep))
                        else:
                            vu = os.path.basename(vu.rstrip(os.sep))
                        pou = os.path.join(xou,vu)
                    t.append('--input')
                    t.append(pin)
                    t.append('--output')
                    t.append(pou)
                    t.append('--label-title')
                    t.append('partial-matched-normal,matched-normal')
                    t.append('--label-file')
                    t.append(partialnormaltemp+','+normaltemp)
                    t.append('--label-threshold')
                    t.append('2,0')
                    t = ' '.join(t)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)
                os.remove(normaltemp)
                os.remove(partialnormaltemp)
            else:
                ix = 0
                for line in txt:
                    ix = ix + 1
                    t = newcmds[:]
                    #print t
                    pin = None
                    pou = None
                    pin = os.path.join(line[0])
                    if len(line) >= 2 and line[1]:
                        pou = os.path.join(xou,line[1])
                    else:
                        #pou = os.path.join(xou,os.path.basename(line[0].rstrip(os.sep)))
                        vu = line[0]
                        if vu.find(',') != -1:
                            vu = str(ix) + "_" + os.path.basename(vu.split(",")[0].rstrip(os.sep))
                        else:
                            vu = os.path.basename(vu.rstrip(os.sep))
                        pou = os.path.join(xou,vu)
                    t.append('--input')
                    t.append(pin)
                    t.append('--output')
                    t.append(pou)
                    t = ' '.join(t)
                    print "------------------------------------------"
                    print t
                    print "------------------------------------------"
                    r = os.system(t)


    else:
        #newcmds.insert(0,head.replace('fusioncatcher-dir.py','fusioncatcher.py'))
        #newcmds = ' '.join(newcmds)
        #os.system(newcmds)
        print """
FUSIONCATCHER-BATCH.PY version 0.12

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2017 Daniel Nicorici

It is a wrapper for 'fusioncatcher.py' which allows to run 'fusioncatcher.py' in
batch mode where the input is a directory which contains the directories for
which FusionCatcher should be run for each of them separately.

The input file (given using -i or --input) may be:
 - a text file with two columns (that are tab separated), where the first column
contains a list of paths (one line = one path = one sample) such that
a path is a SRA file, or a directory containing one or more Fastq files, and the
second column (this is optional and may be missing) contains the corresponding
sample name which will be used to construct the output directory (also URLs may
be given here)
 - a directory which contains subdirectories and each subdirectory contains one
or more Fastq files corresponding to one sample (one subdirectory = one sample);
all found subdirectories will be analyzed one by one by FusionCatcher.

"""
#
