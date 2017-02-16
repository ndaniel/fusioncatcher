#!/usr/bin/env
# -*- coding: utf-8 -*-
"""
It Reading the configuration file: "configuration.cfg".



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
import ConfigParser


def _expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))


def _get_config(c,k,v):
    r = None
    if c.has_option(k,v):
        r = c.get(k,v)
    return r


def _pythonpath(p,c,k,v,skip=False,last=False):
    r = _get_config(c,k,v)
    if r:
        p[v.upper()] = _expand(r)
        sys.path.append(_expand(r))
        if not skip:
            ep = os.getenv('PYTHONPATH')
            if ep:
                if last:
                    os.environ["PYTHONPATH"] = ep + os.pathsep + _expand(r)
                else:
                    os.environ["PYTHONPATH"] = _expand(r) + os.pathsep + ep
            else:
                os.environ["PYTHONPATH"] = _expand(r)


def _envpath(p,c,k,v):
    r = _get_config(c,k,v)
    if r:
        p[v.upper()] = _expand(r)
        ep = os.getenv("PATH")
        if ep:
            os.environ["PATH"] = _expand(r) + os.pathsep + ep
        else:
            os.environ["PATH"] = _expand(r)


def _versions(p,c,k,v):
    r = _get_config(c,k,v)
    if r:
        r = r.strip()
    p[v.upper()] = r


def _parameters(p,c,k,v):
    r = _get_config(c,k,v)
    if r:
        r = r.strip()
    p[v.upper()] = r

def _path(p,c,k,v):
    r = _get_config(c,k,v)
    if r:
        p[v.upper()] = _expand(r)


def manage(configuration_filename, skip_python = []):
    #
    CONF = dict()
    if (not os.path.isfile(configuration_filename)) and (not os.path.islink(configuration_filename)):
        print >> sys.stderr,"WARNING: Configuration file '%s' not found!  Moving on..." % (configuration_filename,)
    else:

        config = ConfigParser.ConfigParser()
        config.read(configuration_filename)

        if "openpyxl" not in skip_python:
            _pythonpath(CONF,config,"paths","openpyxl",last=True)
        if "xlrd" not in skip_python:
            _pythonpath(CONF,config,"paths","xlrd")
        _pythonpath(CONF,config,"paths","numpy",skip=True)
        _pythonpath(CONF,config,"paths","biopython")
        _pythonpath(CONF,config,"paths","scripts")

        _envpath(CONF,config,"paths","java")
        _envpath(CONF,config,"paths","fatotwobit")
        _envpath(CONF,config,"paths","velvet")
        _envpath(CONF,config,"paths","sratoolkit")
        _envpath(CONF,config,"paths","blat")
        _envpath(CONF,config,"paths","liftover")
        _envpath(CONF,config,"paths","python")
        _envpath(CONF,config,"paths","sra")
        _envpath(CONF,config,"paths","lzo")
        _envpath(CONF,config,"paths","lzop")
        _envpath(CONF,config,"paths","coreutils")
        _envpath(CONF,config,"paths","parallel")
        _envpath(CONF,config,"paths","pigz")
        _envpath(CONF,config,"paths","pxz")
        _envpath(CONF,config,"paths","scripts")
        _envpath(CONF,config,"paths","samtools")
        _envpath(CONF,config,"paths","picard")
        _envpath(CONF,config,"paths","bwa")
        _envpath(CONF,config,"paths","bowtie2")
        _envpath(CONF,config,"paths","star")
        _envpath(CONF,config,"paths","bowtie")
        _envpath(CONF,config,"paths","seqtk")
        
        _path(CONF,config,"paths","data")
        _path(CONF,config,"paths","picard")

        _versions(CONF,config,"versions","fusioncatcher")

        _parameters(CONF,config,"parameters","threads")
        _parameters(CONF,config,"parameters","aligners")

    return CONF
#
