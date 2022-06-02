#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script just updates the 'fusioncatcher/etc/configuration.cfg' file with
the paths that it finds when running in the PATH variable.
It only needs to have pre-installed:
- Python version >=2.6.0 and < 3.0.



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
'--skip-blat'. Fore more information regarding BLAT please see its license.

Please, note that FusionCatcher does not require BLAT in order to find
candidate fusion genes!

This file is NOT running/executing/using BLAT.
"""



import os
import sys
import optparse
import shutil
import subprocess
import time
import tempfile




################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

#############################################
# expand path
#############################################
def expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))

#############################################
# Find the path
#############################################
def findpath(exe,title=None,d=""):
    if not title:
        title = exe
    print >>sys.stderr,"Finding",title,"..."
    p = which(exe)
    if p:
        p = os.path.dirname(expand(p))
        print >>sys.stderr,"  * Ok! '%s' found!" % (p,)
    else:
        p = d
        print >>sys.stderr,"  * WARNING: '%s' NOT found!" % (exe,)
    if p:
        p = p.rstrip("/")+"/"
    return p


#############################################
# test Python modules
#############################################
def findmodule(module, title = "",d="", verbose = True):
    """ Test is a given module is installed
        Example:
            module = 'Bio'
            title = 'BioPython'
    """
    if not title:
        title = module
    if verbose:
        print >>sys.stderr,"Checking if the Python library named '%s' is installed..." % (title,)
    flag = True
    try:
        __import__(module)
    except:
        flag = False
        if verbose:
            print >>sys.stderr, "  * WARNING: The Python library '%s' is not installed!\n" % (title,)
    module_path = None
    try:
        module_path = getattr(__import__(module),'__path__')[0]
    except:
        if verbose:
            print >>sys.stderr, "  * WARNING: Cannot find the path of the Python library '%s'!" % (title,)

    if verbose:
        if flag:
            if module_path:
                print "  * Ok! Python library '%s' found at '%s'!" % (title,module_path)
            else:
                print "  * Ok! Python library '%s' found!" % (title,)
        else:
            print >>sys.stderr,"  * WARNING! Python library '%s' not found!" % (title,)
    if flag == False or (not module_path):
        module_path = d
    if module_path:
        module_path = module_path.rstrip("/")+"/"
    return module_path



#############################################
# simulates which
#############################################
def which(program, cwd = True):
    """
    Simulates which from Linux
    
    Usage example:
    
    p = which(exe, cwd = False)
    if p:
        p = expand(p)
        
    """
    if os.path.dirname(program):
        if os.access(program,os.X_OK) and os.path.isfile(program):
            return program
    else:
        paths = os.environ["PATH"].split(os.pathsep)
        if cwd:
            paths.append(os.getcwd())
        for path in paths:
            if path:
                p = os.path.join(path.strip('"'),program)
                if os.access(p,os.X_OK) and os.path.isfile(p):
                    return p
    return None



################################################################################
# MAIN
################################################################################
if __name__ == '__main__':




    #command line parsing

    usage = "%prog [options]"
    description = ("This script just updates the 'fusioncatcher/etc/configuration.cfg'\n"+
                   "file with the paths that it finds when running in the PATH variable.")
    version = "%prog 0.99.7c beta"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("-w","--write-changes",
                      action = "store_true",
                      default = False,
                      dest = "write_changes",
                      help = """If specified than the updates/changes will be written to '%s'.""" % ( os.path.abspath(os.path.join(os.path.dirname(expand(__file__)),"..","etc","configuration.cfg")),))

    (options, args) = parser.parse_args()


################################################################################
################################################################################
################################################################################

    os.system("set +e") # make sure that the shell scripts are still executed if there are errors

    PATH = dict()

    ############################################################################
    # Absolute path to the Python executable
    ############################################################################
    print >>sys.stderr,"Obtaining the absolute path of the Python executable..."
    PATH["python"] = expand(sys.executable)
    print >>sys.stderr,"  * Ok! '%s' found!" % (PATH["python"],)


    ############################################################################
    # Python version
    ############################################################################
    print >>sys.stderr,"Checking Python version..."
    version = sys.version_info
    if version >= (2,6) and version < (3,0):
        print >>sys.stderr,"  * Ok! Found Python version: %s.%s" % (version[0],version[1])
    else:
        print >>sys.stderr, "  * ERROR: Found Python version: %s.%s !\n" % (version[0],version[1])
        print >>sys.stderr, "           The Python version should be >=2.6.0 and < 3.0 . If there is another"
        print >>sys.stderr, "           Python version installed you could run again this script using that"
        sys.exit(1)

    ############################################################################
    # Test 64-bit environment
    ############################################################################
    print >>sys.stderr,"Checking if this environment is a 64-bit environment..."
    import struct
    if struct.calcsize("P") * 8 >= 64:
        print >>sys.stderr,"  * Ok! 64-bit environment found."
    else:
        print >>sys.stderr, "  * WARNING: Not a 64-bit environment! 64-bit environment is needed!"

    ############################################################################
    # FUSIONCATCHER
    ############################################################################
    print >>sys.stderr,"Finding FusionCatcher's path..."
    PATH["scripts"] = os.path.dirname(expand(__file__))
    print >>sys.stderr,"  * Ok! '%s' found!" % (PATH["scripts"],)
    PATH["data"] = os.path.abspath(os.path.join(os.path.dirname(expand(__file__)),"..","data","current"))


    ############################################################################
    # BIOPYTHON
    ############################################################################
    PATH["biopython"] = findmodule("Bio","BioPython")

    ############################################################################
    # Python module: XLRD
    ############################################################################
    PATH["xlrd"] = findmodule("xlrd","Xlrd")

    ############################################################################
    # Python module: OPENPYXL
    ############################################################################
    PATH["openpyxl"] = findmodule("openpyxl","OpenPyXL")

    ############################################################################
    # BOWTIE
    ############################################################################
    PATH["bowtie"] = findpath("bowtie")

    ############################################################################
    # BOWTIE2
    ############################################################################
    PATH["bowtie2"] = findpath("bowtie2")


    ############################################################################
    # SRATOOLKIT
    ############################################################################
    PATH["sra"] = findpath("fastq-dump")

    ############################################################################
    # LILTFOVER
    ############################################################################
    PATH["liftover"] = findpath("liftOver")

    ############################################################################
    # BLAT
    ############################################################################
    PATH["blat"] = findpath("blat")
    
    ############################################################################
    # FATOTWOBIT
    ############################################################################
    PATH["fatotwobit"] = findpath("faToTwoBit")

    ############################################################################
    # SEQTK
    ############################################################################
    PATH["seqtk"] = findpath("seqtk")

    ############################################################################
    # STAR
    ############################################################################
    PATH["star"] = findpath("STAR")

    ############################################################################
    # PIGZ
    ############################################################################
    PATH["pigz"] = findpath("pigz")

    ############################################################################
    # EXTRA (SORT and LZOP)
    ############################################################################
    PATH["bwa"] = findpath("bwa")

    ############################################################################
    # SAMTOOLS
    ############################################################################
    PATH["samtools"] = findpath("samtools")

    ############################################################################
    # VELVET
    ############################################################################
    PATH["velvet"] = findpath("velveth")

    ############################################################################
    # PARALLEL
    ############################################################################
    PATH["parallel"] = findpath("parallel")

    ############################################################################
    # JAVA
    ############################################################################
    PATH["java"] = findpath("java")

#    ############################################################################
#    # PICARD
#    ############################################################################
#    PICARD_PATH = findpath("picard")


#    ############################################################################
#    # LZO
#    ############################################################################
#    LZO_PATH = r
#    LZOP_PATH = r

#    ############################################################################
#    # COREUTILS
#    ############################################################################
#    COREUTILS_PATH = r


#    ############################################################################
#    # PXZ
#    ############################################################################
#    PXZ_PATH = r


    if options.write_changes:
        c = os.path.abspath(os.path.join(os.path.dirname(expand(__file__)),"..","etc","configuration.cfg"))
        n = os.path.abspath(os.path.join(os.path.dirname(expand(__file__)),"..","etc","configuration.cfg.bak"))
        print >>sys.stderr, "\n\nWARNING: Writting updates/changes to the configuration file '%s'!\n\n" %(c,)
        d = [line for line in file(c,"r").readlines()]
        # save the original into BAK file
        file(n,"w").writelines(d)
        r = []
        for line in d:
            t = line
            if line and line.rstrip("\r\n") and (not line.startswith("#")) and (not line.startswith("[")):
                #
                x = line.split("=")
                k = x[0].strip()
                v = PATH.get(k,None)
                if v:
                    t = "%s = %s\n" % (k,v)
                    print >>sys.stderr,"  * Changed: %s" % (t.strip(),)
            r.append(t)
                
        # write the changes
        file(c,"w").writelines(r)
    else:
        print >>sys.stderr, "\n\nWARNING: No changes have been done to the configuration file!\n\n"



