#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A bootstrap script to automatically install FusionCatcher <http://github.com/ndaniel/fusioncatcher>.
It only needs to have pre-installed:
- Python version >=2.6.0 and < 3.0
- NumPy <http://pypi.python.org/pypi/numpy>



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

This file is running/executing/using BLAT.
"""



import os
import sys
import optparse
import shutil
import subprocess
import time
import tempfile
import ftplib

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = "\033[1m"
UNDERLINE = '\033[4m'
HIGHLIGHT = '\033[33;7m'


################################################################################
################################################################################
################################################################################

def PATHS(exe = None, prefix = None, installdir = None, internet = True):
    global DEFAULT_ALIGNERS
    global PYTHON_EXE
    global FUSIONCATCHER_PREFIX
    global FUSIONCATCHER_PATH
    global FUSIONCATCHER_BIN
    global FUSIONCATCHER_URL
    global FUSIONCATCHER_DATA
    global FUSIONCATCHER_CURRENT
    global FUSIONCATCHER_ORGANISM
    global FUSIONCATCHER_TOOLS
    global FUSIONCATCHER_CONFIGURATION
    global FUSIONCATCHER_VERSION
    global FUSIONCATCHER_THREADS
    global NUMPY_PATH
    global NUMPY_URL
    global BIOPYTHON_PATH
    global BIOPYTHON_URL
    global XLRD_PATH
    global XLRD_URL
    global OPENPYXL_PATH
    global OPENPYXL_URL
    global BWA_PATH
    global BWA_URL
    global BWA_VERSION
    global BBMAP_PATH
    global BBMAP_URL
    global BBMAP_VERSION
    global SAMTOOLS_PATH
    global SAMTOOLS_URL
    global SAMTOOLS_VERSION
    global SETUPTOOLS_PATH
    global SETUPTOOLS_URL
    global BOWTIE_OLD_PATH
    global BOWTIE_OLD_URL
    global BOWTIE_OLD_VERSION
    global BOWTIE_PATH
    global BOWTIE_URL
    global BOWTIE_VERSION
    global BOWTIE2_PATH
    global BOWTIE2_URL
    global BOWTIE2_VERSION
    global BLAT_PATH
    global BLAT_URL
    global BLAT_VERSION
    global STAR_PATH
    global STAR_URL
    global STAR_VERSION
    global FATOTWOBIT_PATH
    global FATOTWOBIT_URL
    global FATOTWOBIT_VERSION
    global SRATOOLKIT_PATH
    global SRATOOLKIT_URL
    global SRATOOLKIT_VERSION
    global VELVET_PATH
    global VELVET_URL
    global VELVET_VERSION
    global LZO_PATH
    global LZO_URL
    global LZO_VERSION
    global LZOP_PATH
    global LZOP_URL
    global LZOP_VERSION
    global COREUTILS_PATH
    global COREUTILS_URL
    global COREUTILS_VERSION
    global PIGZ_PATH
    global PIGZ_URL
    global PIGZ_VERSION
    global PXZ_PATH
    global PXZ_URL
    global PXZ_VERSION
    global SEQTK_PATH
    global SEQTK_URL
    global SEQTK_VERSION
    global PARALLEL_PATH
    global PARALLEL_URL
    global PARALLEL_VERSION
    global PICARD_PATH
    global PICARD_URL
    global PICARD_VERSION
    global LIFTOVER_PATH
    global LIFTOVER_URL
    global LIFTOVER_VERSION
    global JAVA_PATH
    global ENSEMBL_VERSION



    # python
    PYTHON_EXE = ''
    if exe:
        PYTHON_EXE = expand(exe)
    else:
        PYTHON_EXE = '/usr/bin/python'

    # ALIGNERS
    DEFAULT_ALIGNERS = "blat,star"


    # FUSIONCATCHER
    if prefix:
        FUSIONCATCHER_PREFIX = expand(prefix)
    else:
        FUSIONCATCHER_PREFIX = expand('.')
        if os.getuid() == 0:
            FUSIONCATCHER_PREFIX = expand('/opt')

    if installdir:
        FUSIONCATCHER_PATH = expand(installdir)
        FUSIONCATCHER_PREFIX = expand(os.path.dirname(installdir.rstrip(os.sep)))
    else:
        FUSIONCATCHER_PATH = expand(FUSIONCATCHER_PREFIX,'fusioncatcher')
    
    FUSIONCATCHER_BIN = expand(FUSIONCATCHER_PATH,'bin')
    FUSIONCATCHER_URL = 'http://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v1.00.zip'
    FUSIONCATCHER_VERSION = "1.00"
    FUSIONCATCHER_DATA = expand(FUSIONCATCHER_PATH,'data')
    FUSIONCATCHER_CURRENT = expand(FUSIONCATCHER_DATA,'current')
    FUSIONCATCHER_ORGANISM = 'homo_sapiens'
    FUSIONCATCHER_THREADS = '1'
    FUSIONCATCHER_TOOLS = expand(FUSIONCATCHER_PATH,'tools')
    FUSIONCATCHER_CONFIGURATION = expand(FUSIONCATCHER_BIN,'..','etc','configuration.cfg')
    # numpy
    NUMPY_PATH = os.path.join(FUSIONCATCHER_TOOLS,'numpy')
    NUMPY_URL = 'http://github.com/numpy/numpy/releases/download/v1.13.1/numpy-1.13.1.tar.gz'
    # biopython
    BIOPYTHON_PATH = os.path.join(FUSIONCATCHER_TOOLS,'biopython')
    BIOPYTHON_URL = 'https://pypi.python.org/packages/72/04/73a4bb22fed40eed26c7e1a673ab51778c577afc3d5dd6f1256424a62c35/biopython-1.70.tar.gz'
    # xlrd python
    XLRD_PATH = os.path.join(FUSIONCATCHER_TOOLS,'xlrd')
    XLRD_URL = 'https://pypi.python.org/packages/42/85/25caf967c2d496067489e0bb32df069a8361e1fd96a7e9f35408e56b3aab/xlrd-1.0.0.tar.gz'
    # openpyxl python
    OPENPYXL_PATH = os.path.join(FUSIONCATCHER_TOOLS,'openpyxl')
    OPENPYXL_URL = 'https://pypi.python.org/packages/5a/8b/798a853ef87d505392227b91d598fd0bdfc8552e64020092e262b1ea7d5f/openpyxl-2.5.0a2.tar.gz'
    # setuptools python
    SETUPTOOLS_PATH = os.path.join(FUSIONCATCHER_TOOLS,'setuptools')
    SETUPTOOLS_URL = 'https://pypi.python.org/packages/07/a0/11d3d76df54b9701c0f7bf23ea9b00c61c5e14eb7962bb29aed866a5844e/setuptools-36.2.7.zip'
    # BOWTIE
    BOWTIE_PATH = os.path.join(FUSIONCATCHER_TOOLS,'bowtie')
#    BOWTIE_URL = 'http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.0/bowtie-1.2-linux-x86_64.zip'
#    BOWTIE_VERSION = ('1.2',)
    BOWTIE_URL = 'http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip'
    BOWTIE_VERSION = ('1.1.2',)
    # BOWTIE
    BOWTIE_OLD_PATH = os.path.join(FUSIONCATCHER_TOOLS,'bowtie-old')
    BOWTIE_OLD_URL = 'http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip'
    BOWTIE_OLD_VERSION = ('1.1.2',)
    # BOWTIE2
    BOWTIE2_PATH = os.path.join(FUSIONCATCHER_TOOLS,'bowtie2')
    BOWTIE2_URL = 'http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip'
    BOWTIE2_VERSION = ('2.2.9',)
    # BLAT
    BLAT_PATH = os.path.join(FUSIONCATCHER_TOOLS,'blat')
    BLAT_URL = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat'
    BLAT_VERSION = ('35x1',)
    # STAR
    STAR_PATH = os.path.join(FUSIONCATCHER_TOOLS,'star')
    STAR_URL = 'http://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz'
    STAR_VERSION = ('STAR_2.5.2b',)
   # BWA
    BWA_PATH = os.path.join(FUSIONCATCHER_TOOLS,'bwa')
    BWA_URL = 'http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2'
    BWA_VERSION = ('0.7.10-r789',)
   # BBMAP
    BBMAP_PATH = os.path.join(FUSIONCATCHER_TOOLS,'bbmap')
    BBMAP_URL = 'https://sourceforge.net/projects/bbmap/files/BBMap_37.68.tar.gz'
    BBMAP_VERSION = ('37','37.68')
    # faToTwoBit
    FATOTWOBIT_PATH = os.path.join(FUSIONCATCHER_TOOLS,'fatotwobit')
    FATOTWOBIT_URL = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit'
    # SRATOOLKIT
    SRATOOLKIT_PATH = os.path.join(FUSIONCATCHER_TOOLS,'sratoolkit')
    SRATOOLKIT_URL = 'http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-centos_linux64.tar.gz'
    SRATOOLKIT_VERSION = ('2.3.5-2','2.4.2','2.5.1','2.6.2','2.8.0','2.8.1-3','2.8.2-1')
    VELVET_PATH = os.path.join(FUSIONCATCHER_TOOLS,'velvet')
    VELVET_URL = 'http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz'
    VELVET_VERSION = ('1.2.09','1.2.10')
    # ENSEMBL version
    ENSEMBL_VERSION = ensembl_version(internet = internet)
    # LZO
    LZO_PATH = os.path.join(FUSIONCATCHER_TOOLS,'lzo')
    LZO_URL = 'http://www.oberhumer.com/opensource/lzo/download/lzo-2.08.tar.gz'
    LZO_VERSION = ('v2.08',)
    # LZOP
    LZOP_PATH = os.path.join(FUSIONCATCHER_TOOLS,'lzop')
    LZOP_URL = 'http://www.lzop.org/download/lzop-1.03.tar.gz'
    LZOP_VERSION = ('v1.03',)
    # COREUTILS (for SORT parallel)
    COREUTILS_PATH = os.path.join(FUSIONCATCHER_TOOLS,'coreutils')
    COREUTILS_URL = 'http://ftp.gnu.org/gnu/coreutils/coreutils-8.27.tar.xz'
    COREUTILS_VERSION = ('v8.27',)
    # PIGZ (GZIP parallel)
    PIGZ_PATH = os.path.join(FUSIONCATCHER_TOOLS,'pigz')
    PIGZ_URL = 'http://http.debian.net/debian/pool/main/p/pigz/pigz_2.3.1.orig.tar.gz' #'http://zlib.net/pigz/pigz-2.3.3.tar.gz'
    PIGZ_VERSION = ('2.3','2.3.1','2.3.3')
    # PXZ (XZ parallel)
    PXZ_PATH = os.path.join(FUSIONCATCHER_TOOLS,'pxz')
    PXZ_URL = 'http://jnovy.fedorapeople.org/pxz/pxz-4.999.9beta.20091201git.tar.xz'
    PXZ_VERSION = ('4.999.9beta',)
    # GNU PARALLEL
    PARALLEL_PATH = os.path.join(FUSIONCATCHER_TOOLS,'parallel')
    PARALLEL_URL = 'http://ftp.gnu.org/gnu/parallel/parallel-20170522.tar.bz2'
    PARALLEL_VERSION = ('20170522',)
    # samtools
    SAMTOOLS_PATH = os.path.join(FUSIONCATCHER_TOOLS,'samtools')
    SAMTOOLS_URL = 'http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2'
    SAMTOOLS_VERSION = ('0.1.19-44428cd',)
    # SEQTK
    SEQTK_PATH = os.path.join(FUSIONCATCHER_TOOLS,'seqtk')
    SEQTK_URL = 'http://github.com/ndaniel/seqtk/archive/1.2-r101c.tar.gz'
    SEQTK_VERSION = ('1.2-r101c-dirty',)
    # 'http://github.com/lh3/seqtk/archive/master.zip'
    #SEQTK_URL = 'http://github.com/lh3/seqtk/archive/1.0.tar.gz'
    # PICARD
    PICARD_PATH = os.path.join(FUSIONCATCHER_TOOLS,'picard')
    #PICARD_URL = 'http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip'
    PICARD_URL = 'http://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar'
    PICARD_VERSION = ('2.9.4',)
    # LiftOver
    LIFTOVER_PATH = os.path.join(FUSIONCATCHER_TOOLS,'liftover')
    LIFTOVER_URL = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver'
    LIFTOVER_VERSION = ('move annotations',)
    # JAVA
    JAVA_PATH = "/usr/bin/"



################################################################################
################################################################################
################################################################################

#############################################
# expand path
#############################################
def expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))

#############################################
# get ensembl version
#############################################
def ensembl_version(internet = True):
    last_version = "unknown"
    if internet:
        print "Checking latest version of Ensembl database that is available..."
        list_files = []
        try:
            ftp = ftplib.FTP("ftp.ensembl.org",timeout=10)
            if ftp:
                ftp.login()
                ftp.cwd("pub")
                list_files = [int(el.lstrip('release').lstrip('-')) for el in ftp.nlst() if el.lower().startswith('release-') and el.lstrip('release').lstrip('-').isdigit()]
        except:
            pass

        if list_files:
            last_version = "v"+str(max(list_files))
        else:
            # try again!
            try:
                import subprocess
                p = subprocess.Popen("wget -nv ftp://ftp.ensembl.org/pub/ -O -", stdout=subprocess.PIPE, stderr=None, shell=True)
                result = p.communicate()[0].split()
                result = [el.split('>release-')[1].split('<')[0] for el in result if el.lower().find('>release-')!=-1]
                result = [int(el) for el in result]
                if result:
                    last_version = "v"+str(max(result))
            except:
                pass
        if last_version == 'unknown':
            print "   * Not found! (WARNING: Is the internet connection working?)"
        else:
            print "   * Version %s found!" % (last_version,)
    return last_version




#############################################
# test Python modules
#############################################
def test_module(module, name = "", package = "", web = "", verbose = False, exit = False):
    """ Test is a given module is installed
        Example:
            module = 'Bio'
            name = 'BioPython'
            package = 'python-biopython'
            description = '<http://pypi.python.org/pypi/biopython>'
    """
    if verbose:
        print "Checking if the Python module named '%s' is installed..." % (name,)
    flag = True
    try:
        __import__(module)
    except:
        flag = False
        if verbose:
            if os.getuid() == 0:
                print >>sys.stderr, "  * WARNING: The Python module '%s' is not installed!\n" % (name,)
                print >>sys.stderr, "             Please, install the Python module: %s (see %s for more info)," % (name,web)
                print >>sys.stderr, "             like for example using the commands: "
                print >>sys.stderr, "               sudo apt-get install %s" % (package,)
                print >>sys.stderr, "             or"
                print >>sys.stderr, "               sudo yum install %s" % (package,)
                print >>sys.stderr, "             or"
                print >>sys.stderr, "               sudo easy_install %s" % (package.lstrip('python').lstrip('-'),)
                print >>sys.stderr, ""
            else:
                print >>sys.stderr, "  * WARNING: The Python module '%s' is not installed!\n" % (name,)
                print >>sys.stderr, "             Please, install the Python module: %s (see %s for more info). It is recommended" % (name,web)
                print >>sys.stderr, "             that YOU the admin/root install this module, like for example using the commands: "
                print >>sys.stderr, "               apt-get install %s" % (package,)
                print >>sys.stderr, "             or"
                print >>sys.stderr, "               yum install %s" % (package,)
                print >>sys.stderr, "             or"
                print >>sys.stderr, "               easy_install %s" % (package.lstrip('python').lstrip('-'),)
                print >>sys.stderr, ""

            print >>sys.stderr, "  * HINTS: Please, also make sure that the following are installed also, before installing the above Python library:"
            print >>sys.stderr, "     - Building tools:"
            print >>sys.stderr, "         sudo apt-get install build-essential"
            print >>sys.stderr, "         or"
            print >>sys.stderr, '         sudo yum groupinstall "Development Tools"'
            print >>sys.stderr, '         sudo yum install ncurses'
            print >>sys.stderr, '         sudo yum install ncurses-devel'
            print >>sys.stderr, "           or"
            print >>sys.stderr, '         sudo zypper install --type pattern Basis-Devel'
            print >>sys.stderr, '         sudo zypper in ncurses'
            print >>sys.stderr, '         sudo zypper in ncurses-devel'
            print >>sys.stderr, "     - Python development:"
            print >>sys.stderr, "         sudo apt-get install python-dev"
            print >>sys.stderr, "           or"
            print >>sys.stderr, "         sudo yum install python-devel"
            print >>sys.stderr, "     - GCC:"
            print >>sys.stderr, "         sudo apt-get install gcc"
            print >>sys.stderr, "           or"
            print >>sys.stderr, "         sudo yum install gcc"
            print >>sys.stderr, "     - ZLIB development:"
            print >>sys.stderr, "         sudo apt-get install zlib-dev"
            print >>sys.stderr, "           or"
            print >>sys.stderr, "         sudo yum install zlib-devel"
            print >>sys.stderr, "     - NumPy library:"
            print >>sys.stderr, "         sudo apt-get install python-numpy"
            print >>sys.stderr, "           or"
            print >>sys.stderr, "         sudo yum install python-numpy"
            print >>sys.stderr, "     - BioPython library:"
            print >>sys.stderr, "         sudo apt-get install python-biopython"
            print >>sys.stderr, "           or"
            print >>sys.stderr, "         sudo yum install python-biopython"
            print >>sys.stderr, "     - TBB:"
            print >>sys.stderr, "         sudo apt-get install libtbb-dev"
            print >>sys.stderr, "           or"
            print >>sys.stderr, "         sudo yum install libtbb-dev"
            print >>sys.stderr, ""

        if exit:
            sys.exit(1)
    module_path = None
    try:
        module_path = getattr(__import__(module),'__path__')[0]
    except:
        if verbose:
            print >>sys.stderr, "  * WARNING: Cannot find the path of the Python module '%s'!" % (name,)

    if verbose:
        if flag:
            if module_path:
                print "  * Ok! Python module '%s' found at '%s'!" % (name,module_path)
            else:
                print "  * Ok! Python module '%s' found!" % (name,)
        else:
            print >>sys.stderr,"  * WARNING! Python module '%s' not found!" % (name,)
    return (flag,module_path)


#############################################
# test a given program/tools/software
#############################################
def test_tool(name = "",
              exe = "",
              web = "",
              param = "",
              verbose = False,
              versions = None,
              version_word = 'version',
              exit = False):
    """ Test if a given tools/program is installed
    """
    flag = False
    if verbose:
        print "Checking if '%s' is installed..." % (name,)

    p = which(exe, cwd = False)
    if p and versions:
        #p = os.path.dirname(expand(p))
        p = expand(p)
        if verbose:
            print "  * Found at '%s'!" % (p,)
            print "  * Test running:  '%s %s'" % (p,param)
        flag,r = cmd([[[p,param],False]], exit = False, verbose = False)
        r = [line for line in r if line.lower().find(version_word.lower()) != -1]
        f = False
        v = None
        for el in versions:
            for line in r:
                if line.lower().find(el.lower()) != -1:
                    f = True
                    v = el
                    break
        if f:
            flag = True
            if verbose:
                print "  * Found supported version '%s'!" % (v,)
        elif verbose:
            print >>sys.stderr,"  * WARNING: Unsupported version found!"
            p = None
            flag = False
    elif verbose:
        print >>sys.stderr,"  * WARNING: Not found!"

    return (flag,p)

#############################################
# temp file
#############################################
def give_me_temp_filename(tmp_dir = None):
    if tmp_dir and (not os.path.isdir(tmp_dir)) and (not os.path.islink(tmp_dir)):
        os.makedirs(tmp_dir)
    (ft,ft_name) = tempfile.mkstemp(dir = tmp_dir)
    os.close(ft)
    return ft_name

#############################################
# simulates which
#############################################
def which(program, cwd = True):
    """
    Simulates which from Linux
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

##############################################
# execute command line commands
##############################################
def cmd(cmds = [],
        exit=True,
        verbose = True):
    """
    stderr = os.devnull
    stderr = subprocess.STDOUT
    """
    f = True
    r = []
    rr = []
    for c in cmds:
        if not c:
            continue
        c0 = [el for el in c[0] if el and el.strip()]
        c1 = c[1]
        if verbose:
            print "    # " + ' '.join(c0)

        if c[0][0] == "cd" and len(c[0]) == 2:
            exit = os.chdir(c[0][1])
            if verbose and exit:
                print >>sys.stderr,''
                print >>sys.stderr, "  * ERROR: Unable to change the directory: '%s'!" % (c[0][1],)
            if exit:
                sys.exit(1)
        else:
            anerr = False
            p = 0
            try:
                p = subprocess.Popen(c0, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = c1)
                r = p.communicate()[0].splitlines()
            except:
                anerr = True

            if anerr:
                if verbose:
                    print >>sys.stderr,''
                    print >>sys.stderr, "  * ERROR: Unable to execute: '%s' (shell = %s)!" % (' '.join(c0),str(c1))
                if exit:
                    sys.exit(1)
                break
            elif p and p.returncode != 0:
                f = False
                if verbose:
                    print >>sys.stderr,''
                    print >>sys.stderr, "  * ERROR: Unable to execute: '%s' (shell = %s)!" % (' '.join(c0),str(c1))
                    for line in r:
                        print >>sys.stderr, line
                if exit:
                    sys.exit(1)
                rr.extend(r)
                break
            rr.extend(r)

    return (f,rr)

#############################################
# install a Python modules from source
#############################################
def install_module(package, url, path, exe = '', pythonpath = '', root_aptget_install = False, verbose = True, exit = True):
    """
    module = module's name
    url =
    path = where to be installed
    """
    # wget http://pypi.python.org/packages/source/o/openpyxl/openpyxl-1.6.2.tar.gz
    # tar zxvf openpyxl-1.6.2.tar.gz
    # cd openpyxl-1.6.2
    # python setup.py build
    afile = os.path.basename(url)
    path = os.path.abspath(os.path.expanduser(path))
    short_path = os.path.dirname(path.rstrip(os.path.sep))
    adir = afile.rstrip('.tar.gz').rstrip('.tgz')
    apython = exe
    if not exe:
        apython = 'python'
    pypath = ''
    if pythonpath:
        pypath = 'PYTHONPATH=%s' % (pythonpath,)
    cmds = []
    if os.getuid() == 0 and root_aptget_install:
        if verbose:
            print "Installing Python package '%s' as root..." % (package,)
        f,r = cmd([(['apt-get','--help'],False)], verbose = False, exit = False)
        if f:
            cmds = cmds + [(['apt-get','--yes','install',package],False)]
        else:
            f,r = cmd([(['yum','--help'],False)], verbose = False, exit = False)
            if f:
                cmds = cmds +[(['yum','-y','install',package],False)]
            else:
                f,r = cmds([(['easy_install','--version'],False)], verbose = False, exit = False)
                if f:
                    cmds = cmds +[([apython,'easy_install','-y',package.lstrip('python').lstrip('-')],False)]
                else:
                    f,r = cmds([(['pip','--version'],False)], verbose = False, exit = False)
                    if f:
                        cmds = cmds +[([apython,'pip','-y',package.lstrip('python').lstrip('-')],False)]
                    else:
                        if verbose:
                            print >>sys.stderr, "  * ERROR: No idea how to install the Python module '%s' using other package " % (module,)
                            print >>sys.stderr, "           managers than 'apt-get', 'yum', 'easy_install', or 'pip'!"
                        if exit:
                            sys.exit(1)
    else:
        if verbose:
            print "Installing Python module '%s' locally '%s'..." % (package,path)
        cmds = [(['mkdir','-p', short_path],False),
                (["rm",'-rf', os.path.join(short_path,afile)],False)]
        if url.startswith('http:') or url.startswith('https:') or url.startswith('ftp:'):
            cmds = cmds + [(['wget', url, '-O', os.path.join(short_path,afile),'--no-check-certificate'],False)]
        else:
            cmds = cmds + [(['cp', url, short_path],False)]
        cmds = cmds + [(["rm",'-rf', path],False),
                (["tar", "-xvzf", os.path.join(short_path,afile), "-C", short_path],False),
                (["ln",'-s',os.path.join(short_path,adir),path],False),
                #(["cd",path,";",pypath,apython,"setup.py","build"],True),
                (["cd",os.path.join(short_path,adir),";",pypath,apython,"setup.py","build"],True)
                #([pypath,apython,os.path.join(short_path,adir,'setup.py'),"build"],True)
             ]
    cmd(cmds,
        verbose = verbose,
        exit = exit
        )
    if verbose:
        print "  * Done!"

#############################################
# install a tool (containing executables)
#############################################
def install_tool(name, url, path, verbose = True, exit = True, env_configure = [], custom_install = []):
    """
    url =
    path = where to be installed
    """
    # cd /apps/fusioncatcher/tools
    # wget http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.08.tgz
    # tar zxvf velvet_1.2.08.tgz
    # cd velvet_1.2.08
    # make
    # cd ..
    # ln -s velvet_1.2.08 velvet
    if verbose:
        print "Installing tool '%s' at '%s' from '%s'" % (name,path,url)
    afile = os.path.basename(url)
    path = os.path.abspath(os.path.expanduser(path))
    short_path = os.path.dirname(path.rstrip(os.path.sep))
    adir = afile.replace('.tar.gz','').replace('.tgz','').replace('.zip','').replace('.tar.xz','').replace('.tar.bz2','').replace('.tbz2','')

    decompress = []

    if afile.endswith('.tar.gz') or afile.endswith('.tgz'):
        decompress = [(["tar","--overwrite", "-xvzf", os.path.join(short_path,afile), "-C", short_path],False)]
    elif afile.endswith('.tar.bz2') or afile.endswith('.tbz2'):
        decompress = [(["tar","--overwrite", "-xvjf", os.path.join(short_path,afile), "-C", short_path],False)]
    elif afile.endswith('.zip'):
        decompress = [(["unzip","-o", os.path.join(short_path,afile), "-d", short_path],False)]
    elif afile.endswith('.tar.xz'):
        cwd = os.path.abspath(os.path.expanduser(os.getcwd()))
        decompress = [(["cd",short_path],False),
                      (["unxz","-f","--keep",os.path.join(short_path,afile)],False),
                      (["tar", "xvf", os.path.join(short_path,afile[:-3]), "-C", short_path],False),
                      (["cd",cwd],False),
                      ]
        #unxz -c coreutils-8.22.tar.xz | tar xv

    if decompress:

        cmds = [(['mkdir','-p',short_path],False),
                (["rm",'-rf',os.path.join(short_path,afile)],False)]
        if url.startswith('http:') or url.startswith('https:') or url.startswith('ftp:'):
            cmds = cmds + [(['wget', url, '-O', os.path.join(short_path,afile),'--no-check-certificate'],False)]
        else:
            cmds = cmds + [(['cp', url, short_path],False)]
        cmds = cmds + [(["rm",'-rf',path],False),
                       (["rm",'-rf',os.path.join(short_path,adir)],False)
                      ]
        cmd(cmds,
            verbose = verbose,
            exit = exit
            )

        listdir = set([os.path.join(short_path,el) for el in os.listdir(short_path) if (not el.startswith('.')) and os.path.isdir(os.path.join(short_path,el))])
        #cmds = [ decompress  ]
        cmds = decompress
        cmd(cmds,
            verbose = verbose,
            exit = exit
            )
        newdir = [os.path.join(short_path,el) for el in os.listdir(short_path) if ( (not el.startswith('.')) and
                                                                                    os.path.isdir(os.path.join(short_path,el)) and
                                                                                    not os.path.join(short_path,el) in listdir)]

        if newdir and len(newdir) == 1:
            newdir = newdir[0]
        else:
            if verbose:
                print >>sys.stderr, "ERROR: Cannot detect the directory where the archive has been decompressed! A potential solution is to delete entirely the previously FusionCatcher installation or to install FusionCatcher in a different path!"
            if exit:
                sys.exit(1)

        cmds = []
        if not os.path.isdir(os.path.join(short_path,adir)):
            cmds.append((["rm","-rf",os.path.join(short_path,adir)],False))
            cmds.append((["mv",newdir,os.path.join(short_path,adir)],False))
            cmds.append((["ln",'-s',os.path.join(short_path,adir),path],False))
            newdir = os.path.join(short_path,adir)
        else:
            cmds.append((["ln",'-s',newdir,path],False))

        cmd(cmds,
            verbose = verbose,
            exit = exit
            )
        cmds = []

        if custom_install:
            # custom commands are passed to fix/install the tool
            file(os.path.join(newdir,'custom_install.sh'),'w').writelines([line.rstrip('\r\n')+'\n' for line in custom_install])
            cwd = os.path.abspath(os.path.expanduser(os.getcwd()))
            cmds.append((["cd",newdir],True))
            cmds.append((["chmod","+x","custom_install.sh"],False))
            if env_configure:
                x = env_configure[:]
                x.insert(0,"env")
                cmds.append((x,True))
            else:
                cmds.append((["./custom_install.sh"],False))
            cmds.append((["cd",cwd],False))

            cmd(cmds,
                verbose = verbose,
                exit = exit
                )
            cmds = []
        else:
            if os.path.isfile(os.path.join(newdir,'configure')):
                cwd = os.path.abspath(os.path.expanduser(os.getcwd()))
                cmds.append((["cd",newdir],True))
                if env_configure:
                    x = env_configure[:]
                    x.insert(0,"env")
                    x.append("./configure")
                    cmds.append((x,True))
                else:
                    cmds.append((["./configure"],False))
                cmds.append((["cd",cwd],False))

            cmd(cmds,
                verbose = verbose,
                exit = exit
                )
            cmds = []

            if os.path.isfile(os.path.join(newdir,'Makefile')):
                cmds.append((["make","-C",newdir],False))
            cmds.append((['chmod','-R','+x',path],False))

            cmd(cmds,
                verbose = verbose,
                exit = exit
                )
    else:
        # it is just an executable
        cmds = [(['mkdir','-p',path],False),
                (["rm",'-rf',os.path.join(path,afile)],False)]
        if url.startswith('http:') or url.startswith('https:') or url.startswith('ftp:'):
            cmds = cmds + [(['wget', url, '-O',os.path.join(path,afile),'--no-check-certificate'],False)]
        else:
            cmds = cmds + [(['cp', url, path],False)]
        cmds = cmds + [(['chmod','+x',os.path.join(path,afile)],False)
                ]
        cmd(cmds,
            verbose = verbose,
            exit = exit
            )
    if verbose:
        print "  * Done!"

#############################################
# User input
#############################################
def accept(question = "", yes = True, no = False, exit = False, force = False, skip = False):
    """
    force => forces all answers on yes!
    """
    flag = yes
    if not force:
        if yes and no:
            print >>sys.stderr, "ERROR: YES and NO answer cannot be both default answers!"
            sys.exit(1)
        y = 'y'
        if yes: # default is YES
            y = 'Y'
        n = 'n'
        if no: # default is NO
            n = 'N'
        s = 's'
        while True:
            if skip:
                text = raw_input("%s [%s/%s/%s] (Y = YES, N = NO, S = SKIP): " % (question, y, n, s))
            else:
                text = raw_input("%s [%s/%s]: " % (question, y, n))
            if not text:
                if yes:
                    text = 'y'
                elif no:
                    text = 'n'
                else:
                    continue

            if text.lower() == 'n':
                flag = False
                break
            elif text.lower() == 'y':
                flag = True
                break
            else:
                continue
        if exit and not flag:
            sys.exit(1)
    return flag

#############################################
# test and install Python module
#############################################
def module(module,
           name,
           package,
           web,
           force,
           url,
           path,
           install = False, # forced install
           pythonpath = '',
           root_aptget_install = False):
    r = False
    thepath = None
    if not install:
        (r, p) = test_module(module = module,
                             name = name,
                             package = package,
                             web = web,
                             verbose = True,
                             exit = False)
    if r:
        thepath = p
    else:
        r = accept(question = "  Do you accept to install the Python module '%s' here '%s'?" % (name,path),
                   yes = True,
                   no = False,
                   exit = False,
                   force = force)
        if not r:
            p = expand(raw_input("  Type new path: "))
            if p:
                path = p
            thepath = path
        install_module(package = package,
                       url = url,
                       path = path,
                       exe = PYTHON_EXE,
                       pythonpath = pythonpath,
                       root_aptget_install = root_aptget_install
                      )
    return thepath

#############################################
# test and install software tool module
#############################################
def tool(name,
         exe,
         param,
         web,
         versions,
         force,
         url,
         path,
         install = False,
         version_word = None,
         env_configure = [],
         custom_install = [],
         skip = False):

    r = False
    thepath = None
    p = None
    if not install:
        if version_word:
            (r, p) = test_tool(name = name,
                               exe = exe,
                               param = param,
                               web = web,
                               verbose = True,
                               versions = versions,
                               version_word = version_word,
                               exit = False)
        else:
            (r, p) = test_tool(name = name,
                               exe = exe,
                               param = param,
                               web = web,
                               verbose = True,
                               versions = versions,
                               exit = False)

    if r:
        thepath = p
    else:
        r = accept(question = "  Do you accept to install the %s here: %s ?" % (name,expand(path)),
                   yes = True,
                   no = False,
                   exit = False,
                   force = force)
        if not r:
            p = ''
            if skip:
                p = expand(raw_input("  Type new path (type S for SKIP): "))
            else:
                p = expand(raw_input("  Type new path: "))
                if p:
                    path = p
                thepath = path
        else:
            thepath = path

        if skip and p and p.lower() == 's':
            pass
        else:
            install_tool(name = name,
                         url = url,
                         path = expand(path),
                         verbose = True,
                         env_configure = env_configure,
                         custom_install = custom_install
                        )
    return thepath

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':


    # initializing the PATHS
    internet = True
    test = [el for el in sys.argv if el in ('-l','--local','-f','--local-fusioncatcher')]
    if test:
        internet = False
    PATHS(internet = internet)

    #command line parsing

    usage = "%prog [options]"
    description = ("A bootstrap script to automatically install FusionCatcher\n"+
                  "<http://github.com/ndaniel/fusioncatcher>. It only needs\n"+
                  "to have pre-installed: (i) Python version >=2.6.0 and < 3.0,\n"+
                  "and (ii) NumPy <http://pypi.python.org/pypi/numpy>.")
    version = "%prog 1.00"

    parser = optparse.OptionParser(usage = usage,
                                   description = description,
                                   version = version)

    parser.add_option("-i","--installation-path",
                      action = "store",
                      type = "string",
                      dest = "installation_directory",
                      default = FUSIONCATCHER_PATH,
                      help = """The directory where FusionCatcher will be installed. Default is '%default'.""")

    parser.add_option("-p","--prefix-path",
                      action = "store",
                      type = "string",
                      dest = "prefix_directory",
                      help = """The FusionCatcher will be installed in 'prefix_directory/fusioncatcher'.""")

    parser.add_option("-a","--install-all",
                      action = "store_true",
                      default = False,
                      dest = "install_all",
                      help = """It forcibly installs (1) all the software tools, and (2) the Python modules, which are needed (even if they are already installed).""")

    parser.add_option("-m","--install-all-py",
                      action = "store_true",
                      default = False,
                      dest = "install_all_py",
                      help = """It forcibly only installs all the Python modules needed (even if they are already installed).""")

    parser.add_option("-t","--install-all-tools",
                      action = "store_true",
                      default = False,
                      dest = "install_all_tools",
                      help = """It forcibly only installs all the software tools needed (even if they are already installed).""")

    parser.add_option("-s","--skip-dependencies",
                      action = "store_true",
                      default = False,
                      dest = "skip_install_all",
                      help = """It skips installing and testing all the dependencies, which are the (1) software tools, and (2) Python modules. Only the FusionCatcher scripts will be installed. Use this when there the internet connection is broken or not available.""")

    parser.add_option("-k","--skip-blat",
                      action = "store_true",
                      default = False,
                      dest = "skip_blat",
                      help = """It skips installing (and using) the BLAT aligner.""")


    parser.add_option("-y","--yes",
                      action = "store_true",
                      default = False,
                      dest = "force_yes",
                      help = """It answers automatically all questions with yes.""")

    parser.add_option("-e","--list-dependencies",
                      action = "store_true",
                      default = False,
                      dest = "list_dependencies",
                      help = """It lists all the needed dependencies. No installation is done!""")

    parser.add_option("-b","--build",
                      action = "store_true",
                      default = False,
                      dest = "build",
                      help = """It builds (and also some download is required) the build files for human organism, which are needed to run FusionCatcher. Default value is '%default'.""")

    parser.add_option("-d","--download",
                      action = "store_true",
                      default = False,
                      dest = "download",
                      help = """It downloads from <http://sourceforge.net/projects/fusioncatcher/files/> the build files for human organism, which are needed to run FusionCatcher. Default value is '%default'.""")

    parser.add_option("-r","--root-apt-get-install",
                      action = "store_true",
                      default = False,
                      dest = "root_aptget_install",
                      help = """If specified and 'bootstrap.py' is run as root then 'bootstrap.py' will run 'apt-get install' in order to install some Python libraries. Default value is '%default'.""")


    parser.add_option("-x","--extra",
                      action = "store_true",
                      default = False,
                      dest = "extra",
                      help = """It installs the last version of SORT (which allows several CPUs to be used and compression of temporary files), LZOP compression, PICARD, VELVET, etc. Default value is '%default'.""")

    parser.add_option("-l","--local",
                      action = "store",
                      type = "string",
                      dest = "local",
                      help = """By default all the software tools and Python modules are downloaded using internet. In case that one wishes to proceed with the install by using all the locally pre-downloaded software tools and Python modules this option should be used. It specifies the local path where all the software tools and Python modules are available and no internet connection will be used.""")

    parser.add_option("-f","--local-fusioncatcher",
                      action = "store",
                      type = "string",
                      dest = "local_fusioncatcher",
                      help = """By default scripts belonging to FusionCatcher are downloaded using internet. In case that one wishes to proceed with the installation by using ONLY the FusionCatcher ZIP archive this option should be used. It specifies the local path where FusionCatcher archive is available and no internet connection will be used.""")

    (options, args) = parser.parse_args()


################################################################################
################################################################################
################################################################################

    hints = """
================================================================================
NOTE: On a Ubuntu running these before installing FusionCatcher might make the installation go smoother:
    
sudo apt-get -y install \\
build-essential\\
libncurses5-dev \\
default-jdk \\
gawk \\
gcc \\
g++ \\
bzip2 \\
make \\
cmake \\
automake \\
gzip \\
zip \\
unzip \\
zlib1g-dev \\
zlib1g \\
wget \\
curl \\
pigz \\
tar \\
parallel \\
libtbb-dev \\
libtbb2 \\
python \\
python-dev \\
python-numpy \\
python-biopython \\
python-xlrd \\
python-openpyxl
================================================================================
"""
    print hints
    time.sleep(5) # wait 5 seconds

    os.system("set +e") # make sure that the shell scripts are still executed if there are errors
    v = "human_v90"
    ############################################################################
    # List all dependencies
    ############################################################################
    if options.list_dependencies:
        print "FusionCatcher [REQUIRED]: ",FUSIONCATCHER_URL
        print "NumPy [REQUIRED but strongly recommended to be installed by root]: ",NUMPY_URL
        print "BioPython [REQUIRED but strongly recommended to be installed by root]: ",BIOPYTHON_URL
        print "Python module XLRD [OPTIONAL; needed only one plans to build database indexes from scratch instead of downloading them]: ",XLRD_URL
        print "Python module OpenPyXL [OPTIONAL; needed only one plans to build database indexes from scratch instead of downloading them]: ",OPENPYXL_URL
        print "Python SETUPTOOLS [OPTIONAL; needed to install XLRD and/or OpenPyXL]: ",SETUPTOOLS_URL
        print "Bowtie: [REQUIRED]: ",BOWTIE_URL
        print "Bowtie2: [REQUIRED]: ",BOWTIE2_URL
        print "Blat [REQUIRED]: ",BLAT_URL
        print "LiftOver [REQUIRED]: ",LIFTOVER_URL
        print "FaToTwoBit (from Blat toolbox) [REQUIRED]: ",FATOTWOBIT_URL
        print "SAMTools [REQUIRED]: ",SAMTOOLS_URL
        print "SRAToolKit (from NCBI) [REQUIRED]: ",SRATOOLKIT_URL
        print "STAR [REQUIRED]: ",STAR_URL
        print "BWA [REQUIRED]: ",BWA_URL
        print "SeqTK [REQUIRED]: ",SEQTK_URL
        print "BBMap [REQUIRED]: ",BBMAP_URL
        print "Velvet (de novo assembler) [OPTIONAL]: ",VELVET_URL
        print "Picard (Java-based SAM tools) [OPTIONAL]: ",PICARD_URL
        print "GNU Parallel (shell tool for executing jobs in parallel) [OPTIONAL]: ",PARALLEL_URL
        print "Pre-built database indexes for human [REQUIRED unless one wants to build them from scratch]:"
        print "  * http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.aa" % (v,)
        print "  * http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.ab" % (v,)
        print "  * http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.ac" % (v,)
        print "  * http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.ad" % (v,)

        sys.exit(0)

    ############################################################################
    # Current working directory
    ############################################################################
    print "Current working directory: '%s'" % (os.getcwd(),)

    ############################################################################
    # Absolute path to the Python executable
    ############################################################################
    print "Obtaining the absolute path of the Python executable..."
    PYTHON_EXE = expand(sys.executable)
    print "  * Ok! '%s' found!" % (PYTHON_EXE,)

    ############################################################################
    # Give option to exit to user if it is not the right Python
    ############################################################################
    print "Python used for installation of FusionCatcher: '%s'" % (PYTHON_EXE,)
    r = accept(question = "  Do you accept this Python?",
               yes = True,
               no = False,
               exit = False,
               force = options.force_yes)
    if not r:
        p = expand(raw_input("  Type new path and filename of the Python binary: "))
        print >>sys.stderr,"Now the boostrap.py will be re-launched!"
        x = '"%s" "%s" %s' % (p,expand(__file__),' '.join(sys.argv[1:]))
        xx = os.system(x)
        sys.exit(0)


    ############################################################################
    # Python version
    ############################################################################
    print "Checking Python version..."
    version = sys.version_info
    if version >= (2,6) and version < (3,0):
        print "  * Ok! Found Python version: %s.%s" % (version[0],version[1])
    else:
        print >>sys.stderr, "  * ERROR: Found Python version: %s.%s !\n" % (version[0],version[1])
        print >>sys.stderr, "           The Python version should be >=2.6.0 and < 3.0 . If there is another"
        print >>sys.stderr, "           Python version installed you could run again this script using that"
        print >>sys.stderr, "           Python version, for example: '/some/other/pythonXYZ bootstrap.py' !"
        sys.exit(1)

    ############################################################################
    # Test 64-bit environment
    ############################################################################
    print "Checking if this environment is a 64-bit environment..."
    import struct
    if struct.calcsize("P") * 8 >= 64:
        print "  * Ok! 64-bit environment found."
    else:
        print >>sys.stderr, "  * ERROR: Not a 64-bit environment! 64-bit environment is needed!"
        sys.exit(1)

    ############################################################################
    # Check if WGET is available
    ############################################################################
    r = os.system('wget --help > /dev/null 2>&1')
    if r != 0:
        print >>sys.stderr, "  * ERROR: 'wget' is not available! Please, install 'wget'!"
        sys.exit(1)

    print ""
    print "Installing FusionCatcher from <http://github.com/ndaniel/fusioncatcher>"
    print "------------------------------------------------------------------------"
    print ""

    # validate options
    if options.skip_install_all and (options.install_all_tools or
       options.install_all_py or options.install_all):
        parser.error("Incompatible command line options (which cannot be used simultaneously)!")

    # validate options
    if options.local and options.local_fusioncatcher:
        parser.error("Incompatible command line options (which cannot be used simultaneously)!")

    if options.prefix_directory:
        options.installation_directory = expand(options.prefix_directory,'fusioncatcher')
        PATHS(exe=PYTHON_EXE, installdir = options.installation_directory)


    ############################################################################
    # Download using URL or copy from a local path?
    ############################################################################
    if options.local:
        if os.path.isdir(options.local):
            # modify all URLS
            ol = expand(options.local)
            FUSIONCATCHER_URL = os.path.join(ol, os.path.basename(FUSIONCATCHER_URL))
            NUMPY_URL = os.path.join(ol, os.path.basename(NUMPY_URL))
            BIOPYTHON_URL = os.path.join(ol, os.path.basename(BIOPYTHON_URL))
            XLRD_URL = os.path.join(ol, os.path.basename(XLRD_URL))
            OPENPYXL_URL = os.path.join(ol, os.path.basename(OPENPYXL_URL))
            SETUPTOOLS_URL = os.path.join(ol, os.path.basename(SETUPTOOLS_URL))
            SAMTOOLS_URL = os.path.join(ol, os.path.basename(SAMTOOLS_URL))
            BOWTIE_URL = os.path.join(ol, os.path.basename(BOWTIE_URL))
            BOWTIE2_URL = os.path.join(ol, os.path.basename(BOWTIE2_URL))
            BLAT_URL = os.path.join(ol, os.path.basename(BLAT_URL))
            BBMAP_URL = os.path.join(ol, os.path.basename(BBMAP_URL))
            LIFTOVER_URL = os.path.join(ol, os.path.basename(LIFTOVER_URL))
            FATOTWOBIT_URL = os.path.join(ol, os.path.basename(FATOTWOBIT_URL))
            SRATOOLKIT_URL = os.path.join(ol, os.path.basename(SRATOOLKIT_URL))
            STAR_URL = os.path.join(ol, os.path.basename(STAR_URL))
            BWA_URL = os.path.join(ol, os.path.basename(BWA_URL))
            SEQTK_URL = os.path.join(ol, os.path.basename(SEQTK_URL))
            VELVET_URL = os.path.join(ol, os.path.basename(VELVET_URL))
            PICARD_URL = os.path.join(ol, os.path.basename(PICARD_URL))
            PARALLEL_URL = os.path.join(ol, os.path.basename(PARALLEL_URL))
            # remove options.local from the PATH variable in order tooid conflicts
            ps = []
            for p in os.environ["PATH"].split(os.pathsep):
                if expand(p) == ol:
                    continue
                else:
                    ps.append(p)
            os.environ["PATH"] = os.pathsep.join(ps)
        else:
            print >> sys.stderr,"ERROR: '%s' should be a directory!" % (options.local,)
            sys.exit(1)
    elif options.local_fusioncatcher:
        FUSIONCATCHER_URL = options.local_fusioncatcher

    ############################################################################
    # Detect silently JAVA
    ############################################################################
    p = which("java")
    if p:
        p = os.path.dirname(expand(p))
        JAVA_PATH = p

    ############################################################################
    # Installation path of FusionCatcher
    ############################################################################
    print "Path for installation of FusionCatcher: '%s'" % (expand(options.installation_directory),)
    r = accept(question = "  Do you accept this path (WARNING: some files/directories within this path may be deleted/replaced/updated without warning)?",
               yes = True,
               no = False,
               exit = False,
               force = options.force_yes)
    if r:
        FUSIONCATCHER_PATH = expand(options.installation_directory)
        PATHS(exe = PYTHON_EXE, installdir = FUSIONCATCHER_PATH)
    else:
        p = expand(raw_input("  Type new path: "))
        PATHS(exe = PYTHON_EXE, installdir = p)

    ############################################################################
    # Number of threads for FusionCatcher
    ############################################################################
    print "Default number of threads/CPUs to be used by FusionCatcher (use 0 for using the number of CPUs detected at the runtime): '%s'" % (FUSIONCATCHER_THREADS,)
    r = accept(question = "  Do you accept?",
               yes = True,
               no = False,
               exit = False,
               force = options.force_yes)
    if r:
        FUSIONCATCHER_THREADS = '1'
    else:
        p = raw_input("  Type the new default for number of threads: ")
        FUSIONCATCHER_THREADS = p

    ############################################################################
    # FusionCatcher
    ############################################################################
    install_tool(name = 'FusionCatcher (fusion genes finder in RNA-seq data)',
                 url = FUSIONCATCHER_URL,
                 path = FUSIONCATCHER_BIN,
                 verbose = True,
                 custom_install = ["#!/usr/bin/env bash",
                                   "rm ../bin",
#                                   "rm -rf ../bin",
#                                   "rm -rf ../etc",
#                                   "rm -rf ../test",
#                                   "rm -rf ../doc",
#                                   "rm -rf ../tools",
#                                   "rm -rf ../data",
#                                   "rm -rf ../docker",
#                                   "rm -rf ../VERSION",
#                                   "rm -rf ../README",
#                                   "rm -rf ../README.md",
#                                   "rm -rf ../NEWS",
#                                   "rm -rf ../LICENSE",
#                                   "rm -rf ../DEPENDENCIES",
#                                   "mkdir -p ../bin",
#                                   "mkdir -p ../etc",
#                                   "mkdir -p ../test",
#                                   "mkdir -p ../doc",
#                                   "mkdir -p ../tools",
#                                   "mkdir -p ../data",
#                                   "mkdir -p ../docker",
#                                   "cp -f $(pwd)/etc/* ../etc/",
#                                   "cp -f $(pwd)/bin/* ../bin/",
#                                   "cp -f $(pwd)/test/* ../test/",
#                                   "cp -f $(pwd)/doc/* ../doc/",
#                                   "cp -f $(pwd)/docker/* ../docker/",
                                   "cp -R -f $(pwd)/* ../",
                                   "chmod +x ../bin/*.py",
                                   "chmod +x ../bin/*.sh",
                                   "chmod +r ../etc/configuration.cfg",
                                   "chmod +x ../bin/fusioncatcher*",
                                   "chmod +x ../bin/FC",
                                   "chmod +x ../test/*.sh",
                                   "chmod -R +r ../test/*",
                                   "deleteme=$(pwd)",
                                   "cd ..",
                                   'rm -rf "$deleteme"'
                                   ]
                )


    if os.getuid() == 0: # root
        cmd([
             [["rm","-f","/usr/bin/fusioncatcher"],False],
             [["ln","-s",os.path.join(FUSIONCATCHER_BIN,'fusioncatcher'),"/usr/bin/fusioncatcher"],False],
             [["rm","-f","/usr/bin/fusioncatcher-build"],False],
             [["ln","-s",os.path.join(FUSIONCATCHER_BIN,'fusioncatcher-build'),"/usr/bin/fusioncatcher-build"],False]
            ],
            exit=True,
            verbose = True)

    if not options.skip_install_all:
        ############################################################################
        # NumPy
        ############################################################################
        (r, p) = test_module(module = "numpy",
                             name = "NumPy",
                             package = "python-numpy",
                             web = "<http://pypi.python.org/pypi/numpy>",
                             verbose = True,
                             exit = True)
        if r:
            NUMPY_PATH = p

        ############################################################################
        # BioPython
        ############################################################################
        r = module(module = "Bio",
                   name = "BioPython",
                   package = "python-biopython",
                   web = "<http://pypi.python.org/pypi/biopython>",
                   force = options.force_yes,
                   url = BIOPYTHON_URL,
                   path = BIOPYTHON_PATH,
                   install = options.install_all or options.install_all_py,
                   root_aptget_install = options.root_aptget_install
                   )
        if r:
            BIOPYTHON_PATH = r

        ############################################################################
        # Python module: XLRD
        ############################################################################
        r = module(module = "xlrd",
                   name = "Xlrd",
                   package = "python-xlrd",
                   web = "<http://pypi.python.org/pypi/xlrd>",
                   force = options.force_yes,
                   url = XLRD_URL,
                   path = XLRD_PATH,
                   install = options.install_all or options.install_all_py,
                   root_aptget_install = options.root_aptget_install
                   )
        if r:
            XLRD_PATH = r

        ############################################################################
        # Python module: OPENPYXL & SETUPTOOLS
        ############################################################################
        r = False
        if (not options.install_all) and (not options.install_all_py):
            (r, p) = test_module(module = "openpyxl",
                                 name = "OpenPyXL",
                                 package = "python-openpyxl",
                                 web = "<http://pypi.python.org/pypi/openpyxl>",
                                 verbose = False,
                                 exit = False)
        if r:
            OPENPYXL_PATH = p
        else:
            r = module(module = "setuptools",
                       name = "SetupTools",
                       package = "python-setuptools",
                       web = "<http://pypi.python.org/pypi/setuptools>",
                       force = options.force_yes,
                       url = SETUPTOOLS_URL,
                       path = SETUPTOOLS_PATH,
                       install = options.install_all or options.install_all_py,
                       root_aptget_install = options.root_aptget_install
                      )
            if r:
                SETUPTOOLS_PATH = r

            r = module(module = "openpyxl",
                       name = "OpenPyXL",
                       package = "python-openpyxl",
                       web = "<http://pypi.python.org/pypi/openpyxl>",
                       force = options.force_yes,
                       url = OPENPYXL_URL,
                       path = OPENPYXL_PATH,
                       install = options.install_all or options.install_all_py,
                       pythonpath = SETUPTOOLS_PATH,
                       root_aptget_install = options.root_aptget_install)
            if r:
                OPENPYXL_PATH = r


        ############################################################################
        # BOWTIE
        ############################################################################
        r = tool(name = "BOWTIE (short read aligner)",
                 exe = "bowtie",
                 param = "--version",
                 web = "<http://bowtie-bio.sourceforge.net/index.shtml>",
                 versions = BOWTIE_VERSION,
                 force = options.force_yes,
                 url = BOWTIE_URL,
                 path = BOWTIE_PATH,
                 install = options.install_all or options.install_all_tools)
        if r:
            BOWTIE_PATH = r

        ############################################################################
        # BOWTIE (old version)
        ############################################################################
        r = tool(name = "BOWTIE (short read aligner) -- older version",
                 exe = "bowtie",
                 param = "--version",
                 web = "<http://bowtie-bio.sourceforge.net/index.shtml>",
                 versions = BOWTIE_OLD_VERSION,
                 force = options.force_yes,
                 url = BOWTIE_OLD_URL,
                 path = BOWTIE_OLD_PATH,
                 install = options.install_all or options.install_all_tools)
        if r:
            BOWTIE_OLD_PATH = r

        ############################################################################
        # BOWTIE2
        ############################################################################
        r = tool(name = "BOWTIE2 (short read aligner)",
                 exe = "bowtie2",
                 param = "--version",
                 web = "<http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>",
                 versions = BOWTIE2_VERSION,
                 force = options.force_yes,
                 url = BOWTIE2_URL,
                 path = BOWTIE2_PATH,
                 install = options.install_all or options.install_all_tools)
        if r:
            BOWTIE2_PATH = r


        ############################################################################
        # SRATOOLKIT
        ############################################################################
        r = tool(name = "NCBI SRA Toolkit (SRA System Development Kit)",
                 exe = "fastq-dump",
                 param = "",
                 web = "<http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>",
                 versions = SRATOOLKIT_VERSION,
                 version_word = 'fastq-dump',
                 force = options.force_yes,
                 url = SRATOOLKIT_URL,
                 path = SRATOOLKIT_PATH,
                 install = options.install_all or options.install_all_tools)
        if r:
            SRATOOLKIT_PATH = r

        ############################################################################
        # LiftOver
        ############################################################################
        r = tool(name = "LiftOver (Batch Coordinate Conversion)",
                 exe = "liftOver",
                 param = "",
                 web = "<http://genome.ucsc.edu/cgi-bin/hgLiftOver>",
                 versions = LIFTOVER_VERSION,
                 version_word = 'liftOver',
                 force = options.force_yes,
                 url = LIFTOVER_URL,
                 path = LIFTOVER_PATH,
                 install = options.install_all or options.install_all_tools)
        if r:
            LIFTOVER_PATH = r

        ############################################################################
        # BLAT
        ############################################################################
        if options.skip_blat:
            BLAT_PATH = ""
            FATOTWOBIT_PATH = ""
            DEFAULT_ALIGNERS = "star,bowtie2"
        else:
            blat = False
            if (not options.install_all) and (not options.install_all_py):
                (b,p) = test_tool(name = "BLAT (alignment tool)",
                                  exe = "blat",
                                  web = "<http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>",
                                  param = "",
                                  versions = BLAT_VERSION,
                                  version_word = 'blat v.',
                                  verbose = False,
                                  exit = False)
                if not b:
                    blat = True
            r = tool(name = "BLAT (alignment tool)",
                     exe = "blat",
                     param = "",
                     web = "<http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>",
                     versions = BLAT_VERSION,
                     version_word = 'blat v.',
                     force = options.force_yes,
                     url = BLAT_URL,
                     path = BLAT_PATH,
                     install = options.install_all or options.install_all_tools)
            if r:
                BLAT_PATH = r

            r = tool(name = "FaToTwoBit (companion of BLAT alignment tool)",
                     exe = "faToTwoBit",
                     param = "",
                     web = "<http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>",
                     versions = ('in.fa',),
                     version_word = "faToTwoBit",
                     force = options.force_yes,
                     url = FATOTWOBIT_URL,
                     path = FATOTWOBIT_PATH,
                     install = options.install_all or options.install_all_tools or blat)
            if r:
                FATOTWOBIT_PATH = r

        ############################################################################
        # SEQTK
        ############################################################################
        r = tool(name = "SEQTK (Toolkit for processing sequences in FASTA/Q formats)",
                 exe = "seqtk",
                 param = "",
                 web = "<http://github.com/ndaniel/seqtk/>",
                 versions = SEQTK_VERSION,
                 version_word = "Version:",
                 force = options.force_yes,
                 url = SEQTK_URL,
                 path = SEQTK_PATH,
                 install = options.install_all or options.install_all_tools)
        if r:
            SEQTK_PATH = r

        ############################################################################
        # STAR
        ############################################################################
        r = tool(name = "STAR (alignment tool)",
                 exe = "STAR",
                 param = "--version",
                 web = "<https://github.com/alexdobin/STAR>",
                 versions = STAR_VERSION,
                 version_word = 'STAR_',
                 force = options.force_yes,
                 url = STAR_URL,
                 path = STAR_PATH,
                 install = options.install_all or options.install_all_tools,
                 custom_install = ["#!/usr/bin/env bash",
#                                   "rm -f source/STAR",
                                   "cp -f bin/Linux_x86_64_static/STAR source/STAR",
#                                   "cd source",
#                                   "make",
#                                   "if ! ./STAR --version; then",
#                                   "    rm -f STAR",
#                                   "    cp ../bin/Linux_x86_64_static/STAR .",
#                                   "fi",
                                   "exit 0"])
        if r:
            STAR_PATH = r

        # BBMAP
        r = tool(name = "BBMap short read aligner, and other bioinformatic tools.",
                 exe = "bbmap.sh",
                 param = "--version",
                 web = "<https://sourceforge.net/projects/bbmap/>",
                 versions = BBMAP_VERSION,
                 force = options.force_yes,
                 url = BBMAP_URL,
                 path = BBMAP_PATH,
                 install = options.install_all or options.install_all_tools,
                 skip = True)
        if r:
            BBMAP_PATH = r

        # PICARD (Java-based command-line utilities that manipulate SAM files)
        r = tool(name = "PICARD (Java-based command-line utilities that manipulate SAM files)",
                 exe = "java -jar picard.jar SamToFastq",
                 param = "--version",
                 web = "<http://github.com/broadinstitute/picard/>",
                 versions = PICARD_VERSION,
                 force = options.force_yes,
                 url = PICARD_URL,
                 path = PICARD_PATH,
                 install = options.install_all or options.install_all_tools,
                 skip = True)
        if r:
            PICARD_PATH = r

        ############################################################################
        # EXTRA (SORT and LZOP)
        ############################################################################
        if options.extra:


            ############################################################################
            # BWA
            ############################################################################
            r = tool(name = "BWA (alignment tool)",
                     exe = "bwa",
                     param = "",
                     web = "<http://bio-bwa.sourceforge.net/>",
                     versions = BWA_VERSION,
                     version_word = 'Version',
                     force = options.force_yes,
                     url = BWA_URL,
                     path = BWA_PATH,
                     install = options.install_all or options.install_all_tools)
            if r:
                BWA_PATH = r


            ############################################################################
            # SAMTOOLS
            ############################################################################
            r = tool(name = "SAMTOOLS (tools for manipulating alignments in the SAM format)",
                     exe = "samtools",
                     param = "",
                     web = "<http://samtools.sourceforge.net/>",
                     versions = SAMTOOLS_VERSION,
                     version_word = 'Version:',
                     force = options.force_yes,
                     url = SAMTOOLS_URL,
                     path = SAMTOOLS_PATH,
                     install = options.install_all or options.install_all_tools)
            if r:
                SAMTOOLS_PATH = r

            # Velvet
            r = tool(name = "VELVET (sequence assembler for short reads)",
                     exe = "velveth",
                     param = "--version",
                     web = "<http://www.ebi.ac.uk/~zerbino/velvet/>",
                     versions = VELVET_VERSION,
                     force = options.force_yes,
                     url = VELVET_URL,
                     path = VELVET_PATH,
                     install = options.install_all or options.install_all_tools,
                     skip = True)
            if r:
                VELVET_PATH = r

            # LZO library
            r = tool(name = "LZO (LZO library for LZOP compression)",
                     exe = "lzop",
                     param = "-V",
                     web = "<http://www.oberhumer.com/opensource/lzo/>",
                     versions = LZO_VERSION,
                     version_word = 'lzo',
                     force = options.force_yes,
                     url = LZO_URL,
                     path = LZO_PATH,
                     install = options.install_all or options.install_all_tools,
                     skip = True)
            if r:
                LZO_PATH = r

            # LZOP executable
            r = tool(name = "LZOP compression",
                     exe = "lzop",
                     param = "-V",
                     web = "<http://www.lzop.org/>",
                     versions = LZOP_VERSION,
                     version_word = 'lzop',
                     force = options.force_yes,
                     url = LZOP_URL,
                     path = LZOP_PATH,
                     install = options.install_all or options.install_all_tools,
                     env_configure = ['CPPFLAGS="-I%s"' % (expand(LZO_PATH,'include','lzo'),),
                                      'LDFLAGS="-L%s"' % (expand(LZO_PATH,'src','.libs'),)],
                     skip = True)
            if r:
                LZOP_PATH = r

            # COREUTILS (for new SORT)
            r = tool(name = "COREUTILS (for latest SORT)",
                     exe = "sort",
                     param = "--version",
                     web = "<http://ftp.gnu.org/gnu/coreutils>",
                     versions = COREUTILS_VERSION,
                     version_word = 'coreutils',
                     force = options.force_yes,
                     url = COREUTILS_URL,
                     path = COREUTILS_PATH,
                     install = options.install_all or options.install_all_tools,
                     skip = True)
            if r:
                COREUTILS_PATH = r

            # PIGZ (GZIP parallel)
            r = tool(name = "PIGZ (GZIP parallel)",
                     exe = "pigz",
                     param = "--version",
                     web = "<http://zlib.net/pigz/>",
                     versions = PIGZ_VERSION,
                     version_word = 'pigz',
                     force = options.force_yes,
                     url = PIGZ_URL,
                     path = PIGZ_PATH,
                     install = options.install_all or options.install_all_tools,
                     skip = True)
            if r:
                PIGZ_PATH = r

            # PXZ (XZ parallel)
#            r = tool(name = "PXZ (XZ parallel)",
#                     exe = "pxz",
#                     param = "--version",
#                     web = "<http://jnovy.fedorapeople.org/pxz/>",
#                     versions = PXZ_VERSION,
#                     version_word = 'pxz',
#                     force = options.force_yes,
#                     url = PXZ_URL,
#                     path = PXZ_PATH,
#                     install = options.install_all or options.install_all_tools,
#                     skip = True)
#            if r:
#                PXZ_PATH = r


            # GNU PARALLEL (shell tool for executing jobs in parallel)
            r = tool(name = "GNU PARALLEL (shell tool for executing jobs in parallel)",
                     exe = "parallel",
                     param = "--version",
                     web = "<http://www.gnu.org/software/parallel/>",
                     versions = PARALLEL_VERSION,
                     force = options.force_yes,
                     url = PARALLEL_URL,
                     path = PARALLEL_PATH,
                     install = options.install_all or options.install_all_tools,
                     skip = True)
            if r:
                PARALLEL_PATH = r


    ############################################################################
    # PYTHON SHEBANG
    ############################################################################
    print "Checking the shebang of FusionCatcher Python scripts..."
#    if which('python') != PYTHON_EXE or os.getuid() == 0: # root
    r = accept(question = "  Shall the shebang of all Python scripts, belonging to FusionCatcher, be set/hardcoded to use this '%s'?" % (PYTHON_EXE,),
               yes = True,
               no = False,
               exit = False,
               force = options.force_yes)
    if not r:
        print "  * The shebang of all Python scripts, belonging to FusionCatcher, is '#!/usr/bin/env python'!"
    else:
        print "  * Updating the SHEBANG of Python scripts with '%s'" % (PYTHON_EXE,)
        pies = [os.path.join(FUSIONCATCHER_BIN,f) for f in os.listdir(FUSIONCATCHER_BIN) if f.endswith('.py') and os.path.isfile(os.path.join(FUSIONCATCHER_BIN,f)) and not f.startswith('.')]
        for f in pies:
            d = file(f,'r').readlines()
            if d:
                d[0] = "#!%s\n" % (PYTHON_EXE,)
            file(f,'w').writelines(d)
        print "  * Done!"
    #    else:
    #        print "  * Ok!"

    ############################################################################
    # FUSIONCATCHER CONFIGURATION
    ############################################################################
    print "Updating the configuration file of FusionCatcher..."
    print "  * configuration file '%s'" % (FUSIONCATCHER_CONFIGURATION,)
    
    def update_path(SOME_PATH,executable,subdir='src'):
        # update the SOME_PATH with subdir
        some_var = os.path.join(SOME_PATH,subdir)
        if (SOME_PATH and
            (not os.path.isfile(os.path.join(SOME_PATH,subdir,executable))) and
            os.path.isfile(os.path.join(SOME_PATH,executable))
           ):
            some_var = SOME_PATH
        return some_var

    sra = update_path(SRATOOLKIT_PATH, 'fastq-dump','bin')
    coreutils = update_path(COREUTILS_PATH, 'sort','src')
    lzop = update_path(LZOP_PATH,'lzop','src')
    star = update_path(STAR_PATH,'star','source')
    parallel2 = update_path(PARALLEL_PATH,'parallel','src')
        
#    # update the SRATOOLKIT with 'bin'
#    sra = os.path.join(SRATOOLKIT_PATH,'bin')
#    if (SRATOOLKIT_PATH and
#        (not os.path.isfile(os.path.join(SRATOOLKIT_PATH,'bin','fastq-dump'))) and
#        os.path.isfile(os.path.join(SRATOOLKIT_PATH,'fastq-dump'))
#       ):
#        sra = SRATOOLKIT_PATH
#    # update the COREUTILS with 'src'
#    coreutils = os.path.join(COREUTILS_PATH,'src')
#    if (COREUTILS_PATH and
#        (not os.path.isfile(os.path.join(COREUTILS_PATH,'src','sort'))) and
#        os.path.isfile(os.path.join(COREUTILS_PATH,'sort'))
#       ):
#        coreutils = COREUTILS_PATH
#    # update the LZOP with 'src'
#    lzop = os.path.join(LZOP_PATH,'src')
#    if (LZOP_PATH and
#        (not os.path.isfile(os.path.join(LZOP_PATH,'src','lzop'))) and
#        os.path.isfile(os.path.join(LZOP_PATH,'lzop'))
#       ):
#        lzop = LZOP_PATH
#    # update the STAR with 'source'
#    star = os.path.join(STAR_PATH,'source')
#    if (STAR_PATH and
#        (not os.path.isfile(os.path.join(STAR_PATH,'source','STAR'))) and
#        os.path.isfile(os.path.join(STAR_PATH,'STAR'))
#       ):
#        star = STAR_PATH
#    # update the PARALLEL with 'src'
#    parallel2 = os.path.join(PARALLEL_PATH,'src')
#    if (PARALLEL_PATH and
#        (not os.path.isfile(os.path.join(PARALLEL_PATH,'src','parallel'))) and
#        os.path.isfile(os.path.join(PARALLEL_PATH,'parallel'))
#       ):
#        parallel2 = PARALLEL_PATH

    config_file = FUSIONCATCHER_CONFIGURATION
    data = []
    data.append("[paths]\n")
    data.append("python = %s\n" %(os.path.dirname(PYTHON_EXE),))
    data.append("data = %s\n"%(FUSIONCATCHER_CURRENT,))
    data.append("scripts = %s\n"%(FUSIONCATCHER_BIN,))
    data.append("bowtie = %s\n"%(BOWTIE_PATH,))
    data.append("bowtie2 = %s\n"%(BOWTIE2_PATH,))
    data.append("bwa = %s\n"%(BWA_PATH,))
    data.append("blat = %s\n"%(BLAT_PATH,))
    data.append("liftover = %s\n"%(LIFTOVER_PATH,))
    data.append("star = %s\n"%(star,))
    data.append("velvet = %s\n"%(VELVET_PATH,))
    data.append("fatotwobit = %s\n"%(FATOTWOBIT_PATH,))
    data.append("samtools = %s\n"%(SAMTOOLS_PATH,))
    data.append("seqtk = %s\n"%(SEQTK_PATH,))
    data.append("sra = %s\n"%(sra,))
    data.append("numpy = %s\n"%(NUMPY_PATH,))
    data.append("biopython = %s\n"%(BIOPYTHON_PATH,))
    data.append("xlrd = %s\n"%(XLRD_PATH,))
    data.append("openpyxl = %s\n"%(OPENPYXL_PATH,))
    data.append("lzo = %s\n"%(LZO_PATH,))
    data.append("lzop = %s\n"%(lzop,))
    data.append("coreutils = %s\n"%(coreutils,))
    data.append("pigz = %s\n"%(PIGZ_PATH,))
    data.append("pxz = %s\n"%(PXZ_PATH,))
    data.append("picard = %s\n"%(PICARD_PATH,))
    data.append("bbmap = %s\n"%(BBMAP_PATH,))
    data.append("parallel = %s\n"%(parallel2,))
    data.append("java = %s\n"%(JAVA_PATH,))
    data.append("\n")
    data.append("[parameters]\n")
    data.append("threads = %s\n" % (FUSIONCATCHER_THREADS,))
    data.append("aligners = %s\n" % (DEFAULT_ALIGNERS,))
    data.append("\n")
    data.append("[versions]\n")
    data.append("fusioncatcher = %s\n"%(FUSIONCATCHER_VERSION,))
    data.append("\n")
    
    file(config_file,'w').writelines(data)

    ############################################################################
    ############################################################################
    print "-------------------------------------------------------------------------------"
    print "FusionCatcher is installed here:\n  %s" % (FUSIONCATCHER_PATH,)
    print "FusionCatcher's scripts are here:\n  %s" % (os.path.join(FUSIONCATCHER_BIN),)
    print "FusionCatcher's dependencies and tools are installed here:\n  %s" % (os.path.join(FUSIONCATCHER_TOOLS),)
    print "FusionCatcher's organism data is here:\n  %s" % (os.path.join(FUSIONCATCHER_DATA),)
    print "Run FusionCatcher as following:\n  %s" % (os.path.join(FUSIONCATCHER_BIN,'fusioncatcher'),)
    print "In order to download and build the files for FusionCatcher run the following:\n  %s" % (os.path.join(FUSIONCATCHER_PATH,FUSIONCATCHER_BIN,'fusioncatcher-build'),)
    print ""
    print HIGHLIGHT+"=== Installed successfully! ==="+ENDC
    print ""
    if options.build:
        time.sleep(5)
        cmd([
             [["rm","-rf",os.path.join(FUSIONCATCHER_CURRENT)],False],
             [["rm","-rf",os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION)],False],
             [["ln","-s",os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION),os.path.join(FUSIONCATCHER_CURRENT)],False]
            ],
            exit=True,
            verbose = True)
        c = [os.path.join(FUSIONCATCHER_BIN,'fusioncatcher-build'),"-g",FUSIONCATCHER_ORGANISM,"-o",os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION)]
        print "  # %s" % (' '.join([el.replace(' ','\\ ') for el in c]))
        time.sleep(5)
        subprocess.call(c)
    else:
        print "*****************************************************************"
        print "*  DON'T FORGET to download (or build) the organism's data needed by FusionCatcher to run!"
        print "*****************************************************************"
        print ""
        print "Several options to get the data needed by FusionCatcher are shown below (please try them in this order)!"
        ########################################################################
        print ""
        print "---------------------------------------------------------------------------"
        print "*  OPTION 1: Download the data needed by FusionCatcher from SOURCEFORGE!"
        print "             THIS IS HIGHLY RECOMMENDED"
        print "---------------------------------------------------------------------------"
        print "In order to download the latest human data files needed by FusionCatcher, please run these (it will take several hours):"
        print ""
        txt = []
        txt.append("rm -rf %s" % (FUSIONCATCHER_CURRENT.replace(" ","\\ "),))
        txt.append("rm -f %s.tar.gz.*" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
        txt.append("rm -rf %s" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
        txt.append("rm -f %s/checksums.md5" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
        txt.append("mkdir -p %s" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
        txt.append("ln -s %s %s" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),FUSIONCATCHER_CURRENT.replace(" ","\\ ")))
        txt.append("wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.aa -O %s.tar.gz.aa" % (v,os.path.join(FUSIONCATCHER_DATA.replace(" ","\\ "),v)))
        txt.append("wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.ab -O %s.tar.gz.ab" % (v,os.path.join(FUSIONCATCHER_DATA.replace(" ","\\ "),v)))
        txt.append("wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.ac -O %s.tar.gz.ac" % (v,os.path.join(FUSIONCATCHER_DATA.replace(" ","\\ "),v)))
        txt.append("wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/%s.tar.gz.ad -O %s.tar.gz.ad" % (v,os.path.join(FUSIONCATCHER_DATA.replace(" ","\\ "),v)))
        txt.append("wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/checksums.md5 -O %s/checksums.md5" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
        txt.append("cd %s" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
        txt.append("md5sum -c %s/checksums.md5" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
        txt.append('if [ "$?" -ne "0" ]; then')
        txt.append('  echo -e "\\n\\n\\n\\033[33;7m   ERROR: The downloaded files from above have errors! MD5 checksums do not match! Please, download them again or re-run this script again!   \\033[0m\\n"')
        txt.append('  exit 1')
        txt.append('fi')
        #txt.append("tar  zxvf  %s -C %s" % (v,FUSIONCATCHER_DATA.replace(" ","\\ ")))
        txt.append("cat %s.tar.gz.* > %s.tar.gz" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ ")))
        txt.append("rm -f %s.tar.gz.*" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
        txt.append("if ! tar -xzf %s.tar.gz -C %s; then" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),FUSIONCATCHER_DATA.replace(" ","\\ ")))
        txt.append('    echo -e "\\n\\n\\n\\033[33;7m   ERROR: The downloaded files are corrupted!   \\033[0m\\n"')
        txt.append("    exit 1")
        txt.append("fi")
        txt.append("rm -f %s.tar.gz" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
        txt.append("rm -f %s/checksums.md5" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
        for t in txt:
            print t
        print ""
        f = os.path.join(FUSIONCATCHER_BIN,"download.sh")
        print "All these commands are saved in '%s' file! You may execute '%s' shell script or copy/paste all the previous commands and run them manually in the terminal!" % (f,f)
        print ""
        txt.append("exit 0")
        txt.insert(0,'#!/usr/bin/env bash')
        file(f,'w').writelines([el+'\n' for el in txt])
        os.system('chmod +rx "%s"' % (f,))
        file_download = f
        ########################################################################
#        file_download = f
##        print ""
##        print "---------------------------------------------------------------------------"
##        print "*  OPTION 2: Download the data needed by FusionCatcher from MEGA.CO.NZ!"
##        print "             TRY THIS ONLY IF OPTION 1 DID NOT WORK!"
##        print "---------------------------------------------------------------------------"
##        print "In order to download the latest human data files needed by FusionCatcher, please run these (it will take several hours):"
##        print ""
#        txt = []
#        txt.append("rm  -rf  %s" % (FUSIONCATCHER_CURRENT.replace(" ","\\ "),))
#        txt.append("rm -f %s.tar.gz.*" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
#        txt.append("rm -rf %s" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
#        txt.append("mkdir  -p  %s" % (FUSIONCATCHER_DATA.replace(" ","\\ "),))
#        txt.append("ln  -s  %s  %s" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),FUSIONCATCHER_CURRENT.replace(" ","\\ ")))
#        txt.append("#  ====>>> 1. Download manually the file '%s.tar.gz.aa' to here '%s' using your favourite Internet browser from here:" % (v,FUSIONCATCHER_DATA))
#        txt.append("#  ====>>>            https://mega.co.nz/#!HYMDGQYK!r7nFk27bu-E-hbfE_eQuz0D0Y8g_MMU2rfLF2-LovtQ")
#        txt.append("#  ====>>> 2. Download manually the file '%s.tar.gz.ab' to here '%s' using your favourite Internet browser from here:" % (v,FUSIONCATCHER_DATA))
#        txt.append("#  ====>>>            https://mega.co.nz/#!2BtgwQJa!Xr7M6hn4WLxomHsvLPyF4nyVk7PeFdZTop6EsB8CYMo")
#        txt.append("#  ====>>> 3. Download manually the file '%s.tar.gz.ac' to here '%s' using your favourite Internet browser from here:" % (v,FUSIONCATCHER_DATA))
#        txt.append("#  ====>>>            https://mega.co.nz/#!Dc1k0IKC!VNBIVe6YuuPRkVfWEMsNZRpYKeVohGl1ws2xAnIUEvY")
#        txt.append("#  ====>>> 4. Download manually the file '%s.tar.gz.ad' to here '%s' using your favourite Internet browser from here:" % (v,FUSIONCATCHER_DATA))
#        txt.append("#  ====>>>            https://mega.co.nz/#!CdlggKqS!c8vkDFS-sNTsWBeqeq8sSWLoupr8-56xiBrNhJYkbeA")
#        txt.append("cat %s.tar.gz.* | tar xz -C %s" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),FUSIONCATCHER_DATA.replace(" ","\\ ")))
#        txt.append("rm -f %s.tar.gz.*" % (os.path.join(FUSIONCATCHER_DATA,v).replace(" ","\\ "),))
##        for t in txt:
##            print t
#        txt.append("exit 0")
#        #txt.insert(0,'#!/usr/bin/env bash')
#        f = os.path.join(FUSIONCATCHER_BIN,"mega.sh")
##        print ""
##        print "All these commands are saved in '%s' file! You shall copy/paste all the previous commands (except the URLS which need to be downloaded manually) and run them manually in the terminal!" % (f,)
##        print ""
#        file(f,'w').writelines([el+'\n' for el in txt])
#        os.system('chmod +rx "%s"' % (f,))
        ########################################################################
        print ""
        print "---------------------------------------------------------------------------"
        print "*  OPTION 2: Build yourself the data needed by FusionCatcher!"
        print "             TRY THIS ONLY IF OPTION 1!"
        print "---------------------------------------------------------------------------"
        print "In order to build yourself the latest human data files needed by FusionCatcher, please run these (it will take several hours):"
        print ""
        txt = []
        txt.append("rm  -rf  %s" % (FUSIONCATCHER_CURRENT.replace(" ","\\ "),))
        txt.append("rm  -rf  %s" % (os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION).replace(" ","\\ "),))
        txt.append("mkdir  -p  %s" % (os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION).replace(" ","\\ "),))
        txt.append("ln  -s  %s  %s" % (os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION).replace(" ","\\ "),FUSIONCATCHER_CURRENT.replace(" ","\\ ")))
        txt.append("%s  -g  homo_sapiens  -o %s" % (os.path.join(FUSIONCATCHER_BIN,'fusioncatcher-build').replace(" ","\\ "),os.path.join(FUSIONCATCHER_DATA,ENSEMBL_VERSION).replace(" ","\\ ")))
        for t in txt:
            print t
        txt.append("exit 0")
        txt.insert(0,'#!/usr/bin/env bash')
        f = os.path.join(FUSIONCATCHER_BIN,"build.sh")
        print ""
        print "All these commands are saved in '%s' file! You may execute '%s' shell script or copy/paste all the previous commands and run them manually in the terminal!" % (f,f)
        print ""
        file(f,'w').writelines([el+'\n' for el in txt])
        os.system('chmod +rx "%s"' % (f,))
        file_build = f
        ########################################################################
        print "---------------------------------------------------------------------------"

        forget = False
        if options.download:
            print "---------------------------------------------------------------------------"
            print "Downloading and installing the databases required by FusionCatcher"
            print "---------------------------------------------------------------------------"
            r = os.system(file_download)
            if r:
                print HIGHLIGHT+"ERROR found!"+ENDC
                print >>sys.stderr,"ERROR: Something went wrong during the execution of '%s'. Exit code %d." % (file_download,r)
                sys.exit(1)
            print "--> DONE!"
            forget = True
        elif options.build:
            print "---------------------------------------------------------------------------"
            print "Building, downloading, and installing the databases required by FusionCatcher"
            print "---------------------------------------------------------------------------"
            r = os.system(file_build)
            if r:
                print HIGHLIGHT+"ERROR found!"+ENDC
                print >>sys.stderr,"ERROR: Something went wrong during the execution of '%s'. Exit code %d." % (file_build,r)
                sys.exit(1)
            print "--> DONE!"
            forget = True
        print ""
        print HIGHLIGHT+"--------------> THE END! <---------------------------"+ENDC
        print ""
        print ""
        if not forget:
            print "*****************************************************************"
            print "*  DON'T FORGET to download (or build) the organism's data needed"
            print "   by FusionCatcher to run (see above for options)!"
            print "*****************************************************************"
    #time.sleep(1)

####
# to add LZOP (for sort???)
# wget http://www.oberhumer.com/opensource/lzo/download/lzo-2.06.tar.gz
# tar zxvf lzo-2.06.tar.gz
# ./configure
# make
#
# wget http://www.lzop.org/download/lzop-1.03.tar.gz
# tar zxvf lzop-1.03.tar.gz
# env CPPFLAGS="-I/apps/fusioncatcher/tools/lzo/include/lzo" LDFLAGS="-L/apps/fusioncatcher/tools/lzo/src/.libs" ./configure
# make

####
# install newer version of SORT
# wget http://ftp.gnu.org/gnu/coreutils/coreutils-8.22.tar.xz
# unxz -c coreutils-8.22.tar.xz | tar xv
# ./configure
# make

####
# zlib from source for PIGZ (GZIP parallel)
# CFLAGS=-O3 -Wall -Wextra -I../zlib-1.2.5/ -L../zlib-1.2.5/ (in Makefile of pigz)
#
# wget http://zlib.net/pigz/pigz-2.3.1.tar.gz
# tar zxvf pigz-2.3.1.tar.gz
# make
#
