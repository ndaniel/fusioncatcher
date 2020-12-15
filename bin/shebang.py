#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It changes the shebang for all Python scripts belonging to FusionCatcher.

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
"""

import os
import sys
import optparse


#############################################
# expand path
#############################################
def expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))


################################################################################
################################################################################
################################################################################


def main():

    # this is needed to overides the newlines formatter in OptionParser
    #The default format_epilog strips the newlines (uses textwrap),
    # so you would need to override format_epilog in your parser like this.
    class MyOptionParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    #command line parsing
    usage = "%prog [options]"

    description = """It changes the default shebang to a new one in all Python scripts belonging to FusionCatcher."""

    epilog = """

Author: Daniel Nicorici
Email: Daniel.Nicorici@gmail.com
Copyright (c) 2009-2020 Daniel Nicorici

"""



    version = "%prog 0.01 beta"

    parser = MyOptionParser(usage       = usage,
                            epilog      = epilog,
                            description = description,
                            version     = version
                            )



    parser.add_option("-p","--path",
                      action = "store",
                      type = "string",
                      dest = "path",
                      default = os.path.dirname(expand(__file__)),
                      help = """The path where are FusionCatcher's Python scripts.""")

    parser.add_option("-s","--shebang",
                      action = "store",
                      type = "string",
                      dest = "shebang",
                      help = """The new shebang for Python scripts belonging to FusionCatcher.""")



    parser.add_option("-a","--author",
                      action = "store",
                      type = "string",
                      dest = "author",
                      help = """Daniel Nicorici, E-mail: Daniel.Nicorici@gmail.com""")


    ( options , args ) = parser.parse_args()

    # validate options
    if not (options.shebang and
            options.path
            ):
        parser.print_help()
        parser.error("One of the arguments has not been specified.")


    #
    for p in os.listdir(options.path):
        f = os.path.join(options.path,p)
        if p.endswith('.py') and os.path.isfile(f):
            d = file(f,"r").readlines()
            if d:
                d[0] = options.shebang.strip()+'\n'
                file(f,"w").writelines(d)
                print "Shebang changed for:",f


if __name__ == '__main__':
    main()

#
