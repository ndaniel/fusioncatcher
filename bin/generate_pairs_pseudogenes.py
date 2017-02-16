#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It builds a list of pairs of pseudogenes.



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
import optparse

def clean(x):
    return x.replace('    ',' ').replace('   ',' ').replace('  ',' ').strip()

def remove_last_parantheses(x,p1='(',p2=')'):
    y = clean(x)
    if y.endswith(p2):
        w = y.rfind(p1)
        if w != -1:
            y = clean(y[:w])
    return y

#
#
#
if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It builds a list of pairs of pseudogenes based on the description of genes."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("--input",
                      action="store",
                      type="string",
                      dest="input_genes_descriptions",
                      help="""Input file with genes positions.""")

    parser.add_option("--paralogs",
                      action="store",
                      type="string",
                      dest="input_paralogs",
                      help="""Input file containing the paralogs genes. It is optional and by using the pairs of paralog genes are combined with the pairs of pseudogenes.""")

    parser.add_option("--output",
                      action="store",
                      type="string",
                      dest="output_directory",
                      default = '.',
                      help="""The output directory where the list of pairs of pseudogenes are written. Default is '%default'.""")


    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_genes_descriptions and
            options.output_directory
            ):
        parser.print_help()
        parser.error("Missing argument(s)!")
        sys.exit(1)


    #
    #
    # read the input genes descriptions
    data = [line.rstrip('\r\n').split('\t') for line in file(options.input_genes_descriptions,'r').readlines() if line.rstrip('\r\n') and line.startswith('ENS')]

    #
    #
    #
    paralogs = []
    original_paralogs = set()
    if options.input_paralogs:
        paralogs = [line.rstrip('\r\n').split('\t') for line in file(options.input_paralogs,'r').readlines() if line.rstrip('\r\n')]
        original_paralogs = set(['\t'.join(sorted(line)) for line in paralogs])
        para = dict()
        for line in paralogs:
            if not para.has_key(line[0]):
                para[line[0]] = set()
            if not para.has_key(line[1]):
                para[line[1]] = set()
            para[line[0]].add(line[1])
            para[line[1]].add(line[0])
        paralogs = para

    # dictionary
    d = dict()
    flag = dict()
    for line in data:
        k = line[0] # ensembl gene id
        if not line[1]:
            continue
        v = clean(line[1].replace(',',' ').lower())
        v = remove_last_parantheses(v, p1 = '[', p2=']')
        v = remove_last_parantheses(v, p1 = '(', p2=')')
        if line[1].find('pseudogene') != -1:
            v = clean(v.split('pseudogene')[0])
            v = remove_last_parantheses(v,p1 = '(',p2=')')
            flag[v] = True
        elif v and (not flag.has_key(v)):
            flag[v] = False

        if not v:
            continue
        if not d.has_key(v):
            d[v] = set()
        d[v].add(k)

    # output
    out = []
    for k,v in d.iteritems():
        if v and len(v) == 1 and flag[k]:
            vx = list(v)
            print "Found orphan pseudogene '%s' (%s)" % (k,vx[0])
        if v and len(v) != 1 and flag[k]:
            x = sorted(set(v))
            if paralogs:
                y = set()
                for e in x:
                    if paralogs.has_key(e):
                        y.update(paralogs[e])
                if y:
                    x = sorted(y.union(set(x)))
            # do all the combinations between the paralogs and pairs of pseudogenes
            n = len(x)
            for i in xrange(n-1):
                for j in xrange(i+1,n):
                    if original_paralogs:
                        p = '\t'.join(sorted([x[i],x[j]]))
                        if p in original_paralogs:
                            continue
                    out.append('%s\t%s\n' % (x[i],x[j]))
#        elif n!= 1 and flag[k] == False:
          #print k, v.pop()
    file(os.path.join(options.output_directory,'pairs_pseudogenes.txt'),'w').writelines(out)
