#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
==============================================================================
FusionCatcher-build
==============================================================================
FusionCatcher-build downloads and builds the data necessary for FusionCatcher.


Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2018 Daniel Nicorici

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

This file is executing by default the BLAT aligner.
"""
import sys
import os
import optparse
import multiprocessing
import subprocess
import shutil
import socket
import locale
import workflow
import configuration

bowtie_error = ("Please, check if 'Bowtie' (from <http://bowtie-bio.sourceforge.net/index.shtml>) "+
               "is installed correctly and it is in the corresponding PATH "+
               "(or if 'configuration.cfg' file is set up correctly)!")
bowtie2_error = ("Please, check if 'Bowtie2' (from <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>) "+
               "is installed correctly and it is in the corresponding PATH "+
               "(or if 'configuration.cfg' file is set up correctly)!")

# for sort in linux
locale.setlocale(locale.LC_ALL, 'C')

def expand(*p):
    return os.path.abspath(os.path.expanduser(os.path.join(*p)))

# get the path of this script
pipeline_path = os.path.dirname(expand(sys.argv[0]))

def outdir(*more_paths):
    global out_dir
    return os.path.join(out_dir,*more_paths)

# make sure that a directory ends with path separator such that workflow can
# recognize it as directory
def adir(a_dir):
    if (not a_dir.endswith('\\')) and (not a_dir.endswith('/')):
        a_dir = a_dir + os.sep
    return a_dir

#
#
def islink(alink = None):
    """
    Wrapper for: os.path.islink()
    """
    f = False
    if alink:
        alink = alink[:-1] if alink.endswith(os.sep) else alink
        if os.path.islink(alink):
            f = True
    return f

#
#
def empty(a_file):
    f = True
    if (os.path.isfile(a_file) or islink(a_file)):
        s = os.path.getsize(a_file)
        if s < 1000:
            d = [line for line in file(a_file,'r').readlines() if line.rstrip('\r\n')]
            if d:
                f = False
        else:
            f = False
    return f

#
# command line parsing
#

class MyOptionParser(optparse.OptionParser):
    def format_epilog(self, formatter):
        return self.epilog


if __name__ == '__main__':

    #command line parsing
    usage = "%prog [options] arg"

    epilog = ("\n" +
             "Author: Daniel Nicorici \n" +
             "Email: Daniel.Nicorici@gmail.com \n" +
             "Copyright (c) 2009-2018 Daniel Nicorici \n " +
             "\n")

    description = ("FusionCatcher-build downloads data from Ensembl database and "+
                   "builds the data necessary for FusionCatcher. FusionCather-build "+
                   "needs to be run only once for each organism or when the Ensembl " +
                   "database is updated to a new version. FusionCatcher is able to "+
                   "reuse the data created and built here. FusionCatcher-build is "+
                   "should be run before FusionCatcher "+
                   "Please see the file "+
                   "'version.txt' for information regarding the Ensembl database "+
                   "version, genome version, and organism name used here."
                  )

    version = "%prog 1.20"

    parser = MyOptionParser(
                usage       = usage,
                epilog      = epilog,
                description = description,
                version     = version
             )



    parser.add_option("--output","-o",
                      action = "store",
                      type = "string",
                      dest = "output_directory",
                      help = "The output directory where all the outputs files "+
                             " and directories will be written.")

    parser.add_option("--organism","-g",
                      action = "store",
                      type = "string",
                      dest = "organism",
                      default = 'homo_sapiens',
                      help = "Organism for which the data is downloaded from Ensembl "+
                             "database and built, for example: 'homo_sapiens', "+
                             "'mus_musculus', 'rattus_norvegicus', 'canis_familiaris', " +
                             "etc. "+
                             "Default is '%default'.")

    parser.add_option("--config","-c",
                      action = "store",
                      type = "string",
                      dest = "configuration_filename",
                      default = os.path.join(pipeline_path,"..","etc","configuration.cfg") + "," + os.path.join(pipeline_path,"configuration.cfg"),
                      help = "Configuration file containing the paths to external "+
                             "tools (e.g. Bowtie, etc.) in case that they are not "+
                             "in PATH! "+
                             "Default is '%default'.")

    parser.add_option("--force-paths","-F",
                      action = "store_true",
                      dest = "force_paths",
                      default = False,
                      help = "If it is specified then all external tools and all Python tools "+
                             "will be executed by FusionCatcher by using their corresponding absolute paths, "+
                             "which will be obined from the fusioncatcher/bin/configuration.cfg file. "+
                             "By default no paths are specified when executing tools/scripts. "+
                             "Default is '%default'. ")

    parser.add_option("--web","-w",
                      action = "store",
                      type = "string",
                      dest = "web_ensembl",
                      default = 'www.ensembl.org',
                      help = "Ensembl database web site from where the data is downloaded. "+
                             " e.g. 'www.ensembl.org', 'uswest.ensembl.org', "+
                             "'useast.ensembl.org', 'asia.ensembl.org', etc. "+
                             "Default is '%default'.")

    parser.add_option("--ftp-ensembl","-e",
                      action = "store",
                      type = "string",
                      dest = "ftp_ensembl",
                      default = 'ftp.ensembl.org',
                      help = "Ensembl database FTP site from where the data is downloaded. "+
                             "Default is '%default'.")

    parser.add_option("--ftp-ensembl-path",
                      action = "store",
                      type = "string",
                      dest = "ftp_ensembl_path",
                      help = "The path for Ensembl database FTP site from where the data is downloaded.")

    parser.add_option("--ftp-ucsc","-x",
                      action = "store",
                      type = "string",
                      dest = "ftp_ucsc",
                      default = 'hgdownload.cse.ucsc.edu',
                      help = "UCSC database FTP site from where the data is downloaded. "+
                             "Default is '%default'.")

    parser.add_option("--ftp-ncbi","-n",
                      action = "store",
                      type = "string",
                      dest = "ftp_ncbi",
                      default = 'ftp.ncbi.nlm.nih.gov',
                      help = "NCBI database FTP site from where the data is downloaded. "+
                             "Default is '%default'.")

    parser.add_option("--skip-blat",
                      action = "store_true",
                      dest = "skip_blat",
                      default = False,
                      help = "If it is set then the pipeline will NOT use the BLAT "+
                             "aligner and all options and methods which make use of "+
                             "BLAT will be disabled. "+
                             "BLAT aligner is used by default. Please, note "+
                             "that BLAT license does not allow BLAT to be used for "+
                             "commercial activities. Fore more information "+
                             "regarding BLAT please see its license: "+
                             "<http://users.soe.ucsc.edu/~kent/src/>. "+
                             "Default is '%default'.")

    parser.add_option("--enlarge-genes",
                      action = "store_true",
                      dest = "enlarge_genes",
                      default = False,
                      help = "If it is set then the genes are enlarged (i.e. their introns include also in the transcriptome). "+
                             "Default is '%default'.")

    parser.add_option("--threads","-p",
                      action = "store",
                      type = "int",
                      dest = "processes",
                      default = 0,
                      help = "Number or processes/threads to be used. "+
                             "Default is '%default'.")

    choices = ('cosmic','conjoing','chimerdb2','chimerdb3','ticdb','cgp','cacg')
    parser.add_option("--skip-database",
                      action = "store",
                      type = "string",
                      dest = "skip_database",
                      default = '',
                      help = "If it is set then the pipeline will skip the specified database(s). "+
                             "The choices are ['"+"','".join(choices)+"']. "+
                             "If several databases should be skipped, then their names shall be separated by comma. "
                             "Default is '%default'.")

    #
    # debugging arguments
    #
    parser.add_option("--start","-s",
                      action = "store",
                      type = "int",
                      dest = "start_step",
                      default = 0,
                      help = "It starts executing the workflow from the given "+
                             "step number. This can be used when the pipeline "+
                             "has crashed/stopped and one wants to re-run it from "+
                             "from the step where it stopped without re-running " +
                             "from the beginning the entire pipeline. "+
                             "0 is for restarting automatically and 1 is the first step. "+
                             "This is intended to be used for debugging. "+
                             "Default is '%default'.")

    choices = ('no','crc32','md5','adler32','sha512','sha256')
    parser.add_option("--hash","-l",
                      action = "store",
                      type = "choice",
                      choices = choices,
                      dest = "hash",
                      default = "no",
                      help = "Hash to be used for computing checksum. The choices "+
                             "are ['"+"','".join(choices)+"']. "+
                             "If it is set up to 'no' then no checksum is used and "+
                             "the entire pipeline is executed as a normal shell "+
                             "script. For more information see 'hash_library' in "+
                             "'workflow.py'. "+
                             "This is intended to be used for debugging. "+
                             "Default is '%default'.")

    parser.add_option("--keep","-k",
                      action = "store_true",
                      dest = "keep_temporary_files",
                      default = False,
                      help = "Preserve intermediate files produced during the run. "+
                             "By default, they are NOT deleted upon exit. "+
                             "This is intended to be used for debugging. "+
                             "Default value is '%default'.")

    parser.add_option("--checksums","-u",
                      action = "store",
                      type = "string",
                      dest = "checksums_filename",
                      default = 'checksums.txt',
                      help = "The name of the checksums file. "+
                             "This is intended to be used for debugging. "+
                             "Default value is '%default'. ")



    #command line parsing
    (options, args) = parser.parse_args()

    #
    # validate options
    #
    if (
        (not options.output_directory)
        ):
        parser.print_help()
        print "EXAMPLE:"
        print "========"
        print ""
        print ""
        print "fusioncatcher-build -g homo_sapiens -o /some/data/directory/"
        print ""
        print ""
        print "NOTE:"
        print "'fusioncatcher-build' needs to be run only once (for each organism"
        print "or when the Ensembl database is updated) and 'fusioncatcher'"
        print "will reuse the '/some/data/directory/'."
        print ""
        print ""
        print "ERROR: output directory is not specified!"
        sys.exit(1)

    #
    #
    #
    out_dir = adir(expand(options.output_directory))
    log_file = expand(outdir('fusioncatcher-build.log'))

    # create the ouput directory
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


    # deal with temporary files flag
    temp_flag = 'yes'
    if options.keep_temporary_files or (options.hash != '' and options.hash != 'no'):
        temp_flag = 'no'

    if not options.processes:
        options.processes = multiprocessing.cpu_count()

    #
    # Reading the configuration file: "configuration.cfg"
    #
    config_files = [el for el in options.configuration_filename.split(",") if el and (os.path.isfile(el) or os.path.islink(el))]
    configfile = ''
    if config_files:
        configfile = config_files[0] # first one has priority
    confs = configuration.manage(configfile)
    # check if version of fusioncatcher.py matches the configuration.cfg file
    p = confs.get("FUSIONCATCHER",None)
    if p:
        t = parser.get_version()
        t = t.lower().split(".py")
        if t and len(t) == 2 and t[1].strip() == p.lower():
            pass
        else:
            print >>sys.stderr,"................................................................................"
            print >>sys.stderr,"ERROR: The version of configuration.cfg file does not match the version of the fusioncatcher-build.py!"
            print >>sys.stderr,"Please, fix this!"
            print >>sys.stderr,"................................................................................"
            sys.exit(1)

    # getting absolute paths for the tools and scripts from configuration.cfg
    _B2_ = confs.get("BOWTIE2").rstrip("/")+"/" if options.force_paths else ''
    _BA_ = confs.get("BWA").rstrip("/")+"/" if options.force_paths else ''
    _BE_ = confs.get("BOWTIE").rstrip("/")+"/" if options.force_paths else ''
    _BT_ = confs.get("BLAT").rstrip("/")+"/" if options.force_paths else ''
    _BP_ = confs.get("BBMAP").rstrip("/")+"/" if options.force_paths else ''
    _FC_ = confs.get("SCRIPTS").rstrip("/")+"/" if options.force_paths else ''
    _FT_ = confs.get("FATOTWOBIT").rstrip("/")+"/" if options.force_paths else ''
    _JA_ = confs.get("JAVA").rstrip("/")+"/" if options.force_paths else ''
    _LR_ = confs.get("LIFTOVER").rstrip("/")+"/" if options.force_paths else ''
    _OS_ = confs.get("OASES").rstrip("/")+"/" if options.force_paths else ''
    _PD_ = confs.get("PICARD").rstrip("/")+"/" if options.force_paths else ''
    _PL_ = confs.get("PARALLEL").rstrip("/")+"/" if options.force_paths else ''
    _PZ_ = confs.get("PIGZ").rstrip("/")+"/" if options.force_paths else ''
    _SA_ = confs.get("SRA").rstrip("/")+"/" if options.force_paths else ''
    _SK_ = confs.get("SEQTK").rstrip("/")+"/" if options.force_paths else ''
    _SS_ = confs.get("SAMTOOLS").rstrip("/")+"/" if options.force_paths else ''
    _SR_ = confs.get("STAR").rstrip("/")+"/" if options.force_paths else ''
    _VT_ = confs.get("VELVET").rstrip("/")+"/" if options.force_paths else ''

    # skipped databases
    skip_database = set([el for el in options.skip_database.lower().split(',') if el])

    # initialize the pipeline
    job = workflow.pipeline(
            log_filename       = log_file,
            checksums_filename = options.checksums_filename,
            hash_library       = options.hash,
            start_step         = options.start_step)


    ##############################################################################
    # START
    ##############################################################################

    os.system("set +e") # make sure that the shell scripts are still executed if there are errors

    os.system(_BE_+"bowtie --version | head -1 > '%s'" % (outdir('bowtie_version.txt'),))
    last_line = file(outdir('bowtie_version.txt'),'r').readline().lower().rstrip("\r\n")
    #correct_versions = set(['bowtie-align version 1.2.1','bowtie-align version 1.2.1.1','bowtie-align version 1.2','bowtie version 1.1.2'])
    correct_versions = set(['version 1.2','version 1.1.2','version 1.2.2','version 1.2.3'])
    bowtie121 = False
    if last_line.find("1.2.") != -1:
        bowtie121 = True
    if (last_line not in correct_versions) and not [1 for el in correct_versions if last_line.lower().endswith(el)]:
        print last_line
        job.close()
        os.system("which bowtie > '%s'" % (outdir('bowtie_path.txt'),))
        bowtie_path = file(outdir('bowtie_path.txt'),'r').readline().rstrip("\r\n")
        print >>sys.stderr,("\n\n\nERROR: Wrong version of BOWTIE found ("+bowtie_path+")! It should be: "+', '.join(sorted(correct_versions))+".")
        print >>sys.stderr,("\nERROR: One may specify the path to the correct version in 'fusioncatcher/etc/configuration.cfg',")
        print >>sys.stderr,("\nERROR: like for example change manually the line fusioncatcher/tools/bowtie to fusioncatcher/tools/bowtie-old")
        print >>sys.stderr,("\nERROR: Also it may be that some of Bowtie's dependencies are missing.")
        print >>sys.stderr,("\nERROR: Therefore also make sure that Bowtie's dependencies are installed, like for example:")
        print >>sys.stderr,("\nERROR:    sudo apt-get install libtbb-dev libtbb2 libc6-dev")
        print >>sys.stderr,("\nERROR: or")
        print >>sys.stderr,("\nERROR:    sudo yum install libtbb-devel libtbb2 libc6-devel")
        sys.exit(1)
    os.remove(outdir('bowtie_version.txt'))


    # save version of Python used to analyze this data
    job.add(_FC_+'python_version.py',kind='program')
    job.run(error_message=("Please, check if 'Python' (from "+
        "<http://python.org/>) is installed (or if 'configuration.cfg' file "+
        "is set up correctly)!"))

    # save version of BioPython used to analyze this data
    job.add(_FC_+'biopython_version.py',kind='program')
    job.run(error_message=("Please, check if 'BioPython' (from "+
        "<http://biopython.org/>) is installed (or if 'configuration.cfg' file "+
        "is set up correctly)!"))

    job.add('printf',kind='program')
    job.add('"%s"' % (options.organism.lower()), kind='parameter')
    job.add('>',outdir('organism.txt'),kind='output')
    job.run()

    job.add(_FC_+'get_genome.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.ftp_ensembl,kind='parameter')
    if options.ftp_ensembl_path:
        job.add('--server-path',options.ftp_ensembl_path+'/fasta/',kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('genome.fa'),kind='output',command_line='no')
    job.add('',outdir('mt.fa'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_chrom_lens.py',kind='program')
    job.add('--input_genome',outdir('genome.fa'),kind='output')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('chromosomes_lengths.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_gtf.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.ftp_ensembl,kind='parameter')
    if options.ftp_ensembl_path:
        job.add('--server-path',options.ftp_ensembl_path+'/gtf/',kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('organism.gtf'),kind='output',command_line='no')
    job.run()

#    job.add('get_phix174.py',kind='program')
#    job.add('--server',options.ftp_ncbi,kind='parameter')
#    job.add('--output',out_dir,kind='output',checksum='no')
#    job.add('',outdir('phix174.fa'),kind='output',command_line='no')
#    job.run()

    job.add(_FC_+'get_viruses.py',kind='program')
    job.add('--server',options.ftp_ncbi,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('viruses.fa'),kind='output',command_line='no')
    job.run()

    job.add('sed',kind='program')
    job.add("'s/ /_/g'",kind='parameter')
    job.add('',outdir('viruses.fa'),kind='input')
    job.add('>',outdir('viruses-noblanks.fa'),kind='output')
    job.run()

    job.add(_FC_+'get_synonyms.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.ftp_ensembl,kind='parameter')
    if options.ftp_ensembl_path:
        job.add('--server-path',options.ftp_ensembl_path+'/mysql/',kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('synonyms.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_genes_symbols.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('genes_symbols.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_genes_descriptions.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('descriptions.txt'),kind='output',command_line='no')
    job.run()

    if options.organism == 'homo_sapiens':
        job.add('sed',kind='program')
        job.add('-i',kind='parameter')
        job.add("""'s/glyceraldehyde-3-phosphate\ dehydrogenase/glyceraldehyde\ 3\ phosphate\ dehydrogenase/g'""",kind='parameter')
        job.add('',outdir('descriptions.txt'),kind='output',checksum='no')
        job.run()



    job.add(_FC_+'get_exons_positions.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--threshold_length','150',kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('exons.txt'),kind='output',command_line='no')
    job.add('',outdir('genes.txt'),kind='output',command_line='no')
    job.run()

    job.link(outdir('genes.txt'),outdir('genes_backup.txt'),kind='copy')

    job.add(_FC_+'get_biotypes.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('biotypes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_paralogs.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('paralogs.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_mtrna.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('mtrna.fa'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_rrna.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('rrna.fa'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_trna.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('trna.fa'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_ucsc.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.ftp_ucsc,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('ucsc_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_refseq.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.ftp_ucsc,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('refseq_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_refseq_ensembl.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('refseq_ids.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_gencode.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
#    job.add('--server',options.ftp_ucsc,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('gencode_genes.txt'),kind='output',command_line='no')
    job.run()

    if options.organism == 'homo_sapiens':
        job.add('wget',kind='program')
        job.add('ftp://%s/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' %(options.ftp_ucsc,),kind='parameter')
        job.add('-O',outdir('hg38ToHg19.over.chain.gz'),kind='output')
        job.run()
    else:
        job.add('echo',kind='program')
        job.add('-n','""',kind='parameter')
        job.add('>',outdir('hg38ToHg19.over.chain.gz'),kind='output')
        job.run()


#    if 'hla' not in skip_database:
#        job.add(_FC_+'get_hla.py',kind='program')
#        job.add('--organism',options.organism,kind='parameter')
#        job.add('--output',out_dir,kind='output',checksum='no')
#        job.add('',outdir('hla.fa'),kind='output',command_line='no')
#        job.run()
#        job.run(error_message = ("If this steps fails to run for whatever reason "+
#                "then it can be skipped by re-running fusioncatcher-build with "+
#                "command line option '--skip-database hla'! This database is optional."))
#    else:
#        job.add('printf',kind='program')
#        job.add('">mock-hla\nACGTGGG%sA\n"' % (500*'A',), kind='parameter')
#        job.add('>',outdir('hla.fa'),kind='output')
#        job.run()

    # adds some missing genes to Ensembl Database
    job.add(_FC_+'add_custom_gene.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('exons.txt'),kind='output',command_line='no')
    job.add('',outdir('genes.txt'),kind='output',command_line='no')
    job.add('',outdir('genes_symbols.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_overlapping_genes.py',kind='program')
    job.add('--input_genes',outdir('genes.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--head','ensembl-1_',kind='parameter')
    job.add('',outdir('ensembl-1_fully_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('ensembl-1_partially_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('ensembl-1_same_strand_overlapping_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'remove_custom_feature.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('exons.txt'),kind='output',command_line='no')
    job.add('',outdir('genes.txt'),kind='output',command_line='no')
    job.add('',outdir('organism.gtf'),kind='output',command_line='no')
    job.run()



    # again
    job.add(_FC_+'generate_overlapping_genes.py',kind='program')
    job.add('--input_genes',outdir('genes.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--head','ensembl_',kind='parameter')
    job.add('',outdir('ensembl_fully_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('ensembl_partially_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('ensembl_same_strand_overlapping_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add('LC_ALL=C',kind='program')
    job.add('comm',kind='parameter')
    job.add('-23',kind='parameter')
    job.add('',outdir('ensembl-1_fully_overlapping_genes.txt'),kind='input',temp_path='yes')
    job.add('',outdir('ensembl_fully_overlapping_genes.txt'),kind='input')
    job.add('>',outdir('ensembl_fully_overlapping_genes.dif'),kind='output')
    job.run()

    job.add('LC_ALL=C',kind='program')
    job.add('comm',kind='parameter')
    job.add('-23',kind='parameter')
    job.add('',outdir('ensembl-1_partially_overlapping_genes.txt'),kind='input',temp_path='yes')
    job.add('',outdir('ensembl_partially_overlapping_genes.txt'),kind='input')
    job.add('>',outdir('ensembl_partially_overlapping_genes.dif'),kind='output')
    job.run()

    job.add('LC_ALL=C',kind='program')
    job.add('comm',kind='parameter')
    job.add('-23',kind='parameter')
    job.add('',outdir('ensembl-1_same_strand_overlapping_genes.txt'),kind='input',temp_path='yes')
    job.add('',outdir('ensembl_same_strand_overlapping_genes.txt'),kind='input')
    job.add('>',outdir('ensembl_same_strand_overlapping_genes.dif'),kind='output')
    job.run()


    if options.organism == 'homo_sapiens':
        job.add('printf',kind='program')
        job.add('"ENSG00000047932\tENSG00000047936\nENSG00000134853\tENSG00000145216\n"',kind='parameter') # GOPC-ROS1 & FIP1L1-PDGFRA
        job.add('>',outdir('extra_overlapping-genes.txt'),kind='output')
        job.run()
    else:
        job.add('echo',kind='program')
        job.add('-n','""',kind='parameter')
        job.add('>',outdir('extra_overlapping-genes.txt'),kind='output')
        job.run()


    job.add('cat',kind='program')
    job.add('',outdir('ensembl_same_strand_overlapping_genes.dif'),kind='input',temp_path='yes')
    job.add('',outdir('ensembl_partially_overlapping_genes.dif'),kind='input',temp_path='yes')
    job.add('',outdir('ensembl_fully_overlapping_genes.dif'),kind='input',temp_path='yes')
    job.add('',outdir('extra_overlapping-genes.txt'),kind='input',temp_path='yes')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('ensembl_removed-overlapping-genes.txt'),kind='output')
    job.run()


    job.add(_FC_+'generate_overlapping_genes.py',kind='program')
    job.add('--input_genes',outdir('ucsc_genes.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--head','ucsc-temp_',kind='parameter')
    job.add('',outdir('ucsc-temp_fully_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('ucsc-temp_partially_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('ucsc-temp_same_strand_overlapping_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('ucsc-temp_fully_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('ucsc_fully_overlapping_genes.txt'),kind='output')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('ucsc-temp_partially_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('ucsc_partially_overlapping_genes.txt'),kind='output')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('ucsc-temp_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('ucsc_same_strand_overlapping_genes.txt'),kind='output')
    job.run()



    job.add(_FC_+'generate_overlapping_genes.py',kind='program')
    job.add('--input_genes',outdir('refseq_genes.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--head','refseq-temp_',kind='parameter')
    job.add('',outdir('refseq-temp_fully_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('refseq-temp_partially_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('refseq-temp_same_strand_overlapping_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('refseq-temp_fully_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('refseq_fully_overlapping_genes.txt'),kind='output')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('refseq-temp_partially_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('refseq_partially_overlapping_genes.txt'),kind='output')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('refseq-temp_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('refseq_same_strand_overlapping_genes.txt'),kind='output')
    job.run()



    job.add(_FC_+'generate_overlapping_genes.py',kind='program')
    job.add('--input_genes',outdir('gencode_genes.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--head','gencode-temp_',kind='parameter')
    job.add('',outdir('gencode-temp_fully_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('gencode-temp_partially_overlapping_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('gencode-temp_same_strand_overlapping_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('gencode-temp_fully_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('gencode_fully_overlapping_genes.txt'),kind='output')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('gencode-temp_partially_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('gencode_partially_overlapping_genes.txt'),kind='output')
    job.run()

    job.add(_FC_+'convert_pairs_symbols.py',kind='program')
    job.add('--input',outdir('gencode-temp_same_strand_overlapping_genes.txt'),kind='input')
    job.add('--filter',outdir('ensembl_removed-overlapping-genes.txt'),kind='input')
    job.add('--output',outdir('gencode_same_strand_overlapping_genes.txt'),kind='output')
    job.run()


    job.add('cat',kind='program')
    job.add('',outdir('ensembl_partially_overlapping_genes.txt'),kind='input')
    job.add('',outdir('ensembl_fully_overlapping_genes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('ensembl_overlapping_genes.txt'),kind='output')
    job.run()

    if options.organism == 'homo_sapiens':
        job.add('cut',kind='program')
        job.add('-f1,2',kind='parameter')
        job.add('',outdir('genes_symbols.txt'),kind='input')
        job.add('|',kind='parameter')
        job.add('grep',kind='parameter')
        job.add('','"\tIG[K|L|H|J]"',kind='parameter')
        job.add('|',kind='parameter')
        job.add('cut',kind='parameter')
        job.add('-f1',kind='parameter')
        job.add('|',kind='parameter')
        job.add('LC_ALL=C',kind='parameter')
        job.add('sort',kind='parameter')
        job.add('|',kind='parameter')
        job.add('uniq',kind='parameter')
        job.add('>',outdir('ig_loci.txt'),kind='output')
        job.run()

    else:
        job.add('echo',kind='program')
        job.add('-n','""',kind='parameter')
        job.add('>',outdir('ig_loci.txt'),kind='output')
        job.run()

    job.add(_FC_+'add_gap_gene.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('gap_fusions.txt'),kind='output',command_line='no')
    job.run()

#    if options.organism == 'homo_sapiens':
#        job.add('grep',kind='program')
#        job.add('-i',kind='parameter')
#        job.add('-e','"\tEPOR"',kind='parameter')
#        job.add('-e','"\tCRLF2"',kind='parameter')
#        job.add('',outdir('genes_symbols.txt'),kind='input')
#        job.add('|',kind='parameter')
#        job.add('cut',kind='parameter')
#        job.add('-f1',kind='parameter')
#        job.add('|',kind='parameter')
#        job.add('uniq',kind='parameter')
#        job.add('>',outdir('eporcrlf2.txt'),kind='output')
#        job.run()

#    else:
#        job.add('echo',kind='program')
#        job.add('-n','""',kind='parameter')
#        job.add('>',outdir('eporcrlf2.txt'),kind='output')
#        job.run()

    if options.enlarge_genes:
        job.add(_FC_+'enlarge_genes.py',kind='program')
        job.add('--enlargement-size','5000',kind='parameter') # 5000
        job.add('--full-cover','10000000',kind='parameter')
        job.add('--gene-short','250',kind='parameter')
        #job.add('--genes',outdir('ig_loci.txt'),kind='input') # enlarge only IG loci
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('exons.txt'),kind='output',command_line='no')
        job.add('',outdir('genes.txt'),kind='output',command_line='no')
        job.run()
    elif options.organism == 'homo_sapiens':
        job.add(_FC_+'generate_enlarge.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('enlarge.txt'),kind='output',command_line='no')
        job.run()

        job.add(_FC_+'enlarge_genes.py',kind='program')
        job.add('--enlargement-size','5000',kind='parameter') # 5000
        job.add('--full-cover','10000000',kind='parameter')
        job.add('--gene-short','250',kind='parameter')
        job.add('--genes',outdir('enlarge.txt'),kind='input') # enlarge only the targeted genes
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('exons.txt'),kind='output',command_line='no')
        job.add('',outdir('genes.txt'),kind='output',command_line='no')
        job.run()

    job.add(_FC_+'generate_known.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('known.txt'),kind='output',command_line='no')
    job.run()

    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','pseudogene',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('pseudogenes.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','ribosomal',kind='parameter')
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','protein',kind='parameter')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('ribosomal_proteins.txt'),kind='output')
    job.run()

    job.add(_FC_+'shield_genes.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--read-len','500',kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--pseudo-genes-check',kind='parameter')
    job.add('',outdir('shielded_genes.bed'),kind='output',command_line='no')
    job.add('',outdir('shielded_reads_similarity.bed'),kind='output',command_line='no')
    job.run()

    if not empty(outdir('shield_erase-regions.bed')):

        job.add(_SK_+'seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-n','A',kind='parameter')
        #job.add('-M',outdir('shield_against_pseudo-genes.bed'),kind='input') # this erases only the pseudogenes associated with a manually given list of genes
        job.add('-M',outdir('shield_pseudogenes-predicted-to-be-erased.bed'),kind='input') # # this erases ALL the pseudogenes associated with genes (based on shielding..._report.txt file)
        job.add('',outdir('genome.fa'),kind='input')
        job.add('>',outdir('genome2.fa'),kind='output')
        job.run()
        
        job.add(_SK_+'seqtk',kind='program')
        job.add('seq',kind='parameter')
        job.add('-n','A',kind='parameter')
        job.add('-M',outdir('shield_erase-regions.bed'),kind='input')
        job.add('',outdir('genome2.fa'),kind='input',temp_path='yes')
        job.add('>',outdir('genome3.fa'),kind='output')
        job.run()

        job.link(outdir('genome3.fa'),outdir('genome.fa'),temp_path='yes')

#        job.add('seqtk',kind='program')
#        job.add('seq',kind='parameter')
#        job.add('-n','A',kind='parameter')
#        job.add('-M',outdir('shielded_genes.bed'),kind='input')
#        job.add('',outdir('genome.fa'),kind='input')
#        job.add('>',outdir('genome2.fa'),kind='output')
#        job.run()

#        job.add('seqtk',kind='program')
#        job.add('subseq',kind='parameter')
#        job.add('',outdir('genome.fa'),kind='input')
#        job.add('',outdir('shielded_reads_similarity.bed'),kind='input',temp_path='yes')
#        job.add('>',outdir('shielded_reads.fa'),kind='output')
#        job.run()

#        job.add('bowtie-build',kind='program')
#        job.add('-f',kind='parameter')
#        job.add('--quiet',kind='parameter')
#        job.add('--offrate','1',kind='parameter')
#        job.add('--ftabchars','7',kind='parameter')
#        job.add('',outdir('genome2.fa'),kind='input',temp_path='yes')
#        job.add('',outdir('genome2_index/'),kind='output')
#        job.run(error_message = bowtie_error)

#        job.add('bowtie',kind='program')
#        job.add('-p',options.processes,kind='parameter',checksum='no')
#        job.add('-f',kind='parameter') # fasta
#        #job.add('-a',kind='parameter')
#        job.add('-k','1000',kind='parameter')
#        job.add('-v','0',kind='parameter') # mismatches
#        #job.add('--best',kind='parameter')
#        #job.add('--sam',kind='parameter')
#        job.add('--chunkmbs','128',kind='parameter',checksum='no')
#        job.add('--suppress','1,2,6,7,8',kind='parameter') # original
#        if os.path.isfile(outdir('genome2_index','.1.ebwtl')):
#            job.add('--large-index',kind='parameter')
#        job.add('',outdir('genome2_index/'),kind='input',temp_path='yes')
#        job.add('',outdir('shielded_reads.fa'),kind='input',temp_path='yes')
#        job.add('',outdir('shield.map'),kind='output')
#        job.run()

#        job.add('awk',kind='program')
#        job.add("""'{print $1"\\t"$2"\\t"$2+length($3)}'""",kind='parameter')
#        job.add('',outdir('shield.map'),kind='input',temp_path='yes')
#        job.add('|',kind='parameter') # XXX
#        job.add('LC_ALL=C',kind='parameter')
#        job.add('sort',kind='parameter')
#        job.add('-t',"'\t'",kind='parameter')
#        job.add('-k','1,1',kind='parameter')
#        job.add('-k','2,2n',kind='parameter')
#        job.add('--buffer-size',"80%",kind='parameter',checksum='no')
#        job.add('|',kind='parameter') # XXX
#        job.add('uniq',kind='parameter') # XXX
#        job.add('|',kind='parameter') # XXX
#        job.add('clean_bed.py',kind='parameter') # XXX
#        job.add('-i','-',kind='parameter') # XXX
#        job.add('-o',outdir('shield_mask.bed'),kind='output')
#        job.run()

#        job.add('seqtk',kind='program')
#        job.add('seq',kind='parameter')
#        job.add('-n','A',kind='parameter')
#        job.add('-M',outdir('shield_mask.bed'),kind='input')
#        job.add('',outdir('genome.fa'),kind='input')
#        job.add('>',outdir('genome3.fa'),kind='output')
#        job.run()
#        
#        job.link(outdir('genome3.fa'),outdir('genome.fa'),temp_path='yes')


    job.add(_FC_+'generate_rrna_unit.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('rrna_unit.fa'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_banned.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('banned.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_healthy.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('healthy.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_nontumor.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('non-tumor_cells.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_noncancer.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('non-cancer_tissues.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_cancer-genes.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('cancer_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_oncogenes.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('oncogenes_more.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_tumor-genes.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('tumor_genes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_exons.py',kind='program')
    job.add('--input_exons',outdir('exons.txt'),kind='input')
    job.add('--input_genome',outdir('genome.fa'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('exons.fa'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_genes.py',kind='program')
    job.add('--input_genes',outdir('genes.txt'),kind='input')
    job.add('--input_genome',outdir('genome.fa'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('genes.fa'),kind='output',command_line='no')
    job.run()


    job.add(_FC_+'generate_adjacent_genes.py',kind='program')
    job.add('--input_genes',outdir('genes.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('adjacent_genes.txt'),kind='output',command_line='no')
    job.add('',outdir('readthroughs.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_pairs_pseudogenes.py',kind='program')
    job.add('--input',outdir('descriptions.txt'),kind='input')
    job.add('--paralogs',outdir('paralogs.txt'),kind='input')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('pairs_pseudogenes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_genes_with_no_proteins.py',kind='program')
    job.add('--input',outdir('exons.txt'),kind='input')
    job.add('--output',outdir('genes_with_no_proteins.txt'),kind='output')
    job.add('',outdir('genes_with_no_proteins.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_tcga.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('tcga.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_tcga.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('tcga-normal.txt'),kind='output',command_line='no')
    job.add('',outdir('tcga-cancer.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_tcga2.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('tcga2.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_bodymap2.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('bodymap2.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_cortex.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('cortex.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_hpa.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('hpa.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_1000genomes.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('1000genomes.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_18cancers.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('18cancers.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_oesophagus.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('oesophagus.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'generate_gliomas.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('--skip-filter-overlap',out_dir,kind='parameter')
    job.add('',outdir('gliomas.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_celllines.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('celllines.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_prostate_cancer.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('prostate_cancer.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_pancreases.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('pancreases.txt'),kind='output',command_line='no')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case','\tHLA-',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('hla2.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case','\tHLA-',kind='parameter')
    job.add('',outdir('synonyms.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('hla_2.txt'),kind='output')
    job.run()

    job.add('cat',kind='program')
    job.add('',outdir('hla2.txt'),kind='input')
    job.add('',outdir('hla_2.txt'),kind='input',temp_path='yes')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1,1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('hla.txt'),kind='output')
    job.run()

    job.add(_FC_+'get_hla2.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--server',options.web_ensembl,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('hla2.fa'),kind='output',command_line='no')
    job.run()

    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','antisense',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('antisenses.txt'),kind='output')
    job.run()

    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','rRNA',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('rrnas.txt'),kind='output')
    job.run()

    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','tRNA',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('trnas.txt'),kind='output')
    job.run()

    job.add('cut',kind='program')
    job.add('-f1,5',kind='parameter')
    job.add('',outdir('genes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','MT',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('mt.txt'),kind='output')
    job.run()

    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','lncRNA',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('lncrnas.txt'),kind='output')
    job.run()


    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','microrna',kind='parameter') #
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('>',outdir('temp.txt'),kind='output')
    job.run()
    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','mirna',kind='parameter') #
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('>>',outdir('temp.txt'),kind='output')
    job.run()
    job.add('cat',kind='program')
    job.add('',outdir('temp.txt'),kind='input',temp_path='yes')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('mirnas.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','"small nuclear rna"',kind='parameter') #
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('>',outdir('temp.txt'),kind='output')
    job.run()
    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','snrna',kind='parameter') #
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('>>',outdir('temp.txt'),kind='output')
    job.run()
    job.add('cat',kind='program')
    job.add('',outdir('temp.txt'),kind='input',temp_path='yes')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('snrnas.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','"small nucleolar rna"',kind='parameter') #
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('>',outdir('temp.txt'),kind='output')
    job.run()
    job.add('cut',kind='program')
    job.add('-f1,2',kind='parameter')
    job.add('',outdir('biotypes.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('--ignore-case',kind='parameter')
    job.add('','snorna',kind='parameter') #
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('>>',outdir('temp.txt'),kind='output')
    job.run()
    job.add('cat',kind='program')
    job.add('',outdir('temp.txt'),kind='input',temp_path='yes')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('snornas.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','"\ty rna"',kind='parameter')
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('rnas_y.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','"7SK "',kind='parameter')
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('7skrnas.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','oncogene',kind='parameter')
    job.add('',outdir('descriptions.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('oncogenes.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('','RP11-',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('rp11.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('--ignore-case',kind='parameter')
    job.add('','metazoa',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('metazoa.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('','CTA-',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('cta.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('','CTB-',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('ctb.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('','CTD-',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('ctd.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('','CTC-',kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('ctc.txt'),kind='output')
    job.run()

    job.add('grep',kind='program')
    job.add('-E',"'\\bRP[0-9]*-\\b'",kind='parameter')
    job.add('',outdir('genes_symbols.txt'),kind='input')
    job.add('|',kind='parameter')
    job.add('grep',kind='parameter')
    job.add('','ENS',kind='parameter')
    job.add('|',kind='parameter')
    job.add('cut',kind='parameter')
    job.add('-f1',kind='parameter')
    job.add('|',kind='parameter')
    job.add('LC_ALL=C',kind='parameter')
    job.add('sort',kind='parameter')
    job.add('|',kind='parameter')
    job.add('uniq',kind='parameter')
    job.add('>',outdir('rp.txt'),kind='output')
    job.run()

    job.add(_FC_+'generate_transcripts.py',kind='program')
    job.add('--input_fasta_exons',outdir('exons.fa'),kind='input')
    job.add('--input_database',outdir('exons.txt'),kind='input')
    job.add('--skip',outdir('mirnas.txt'),kind='input')
    job.add('--threshold_length','150',kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('transcripts.fa'),kind='output',command_line='no')
    job.run()

#    job.add('extract_transcripts.py',kind='program')
#    job.add('--input_genes',outdir('hla_12.txt'),kind='input')
#    job.add('--input_transcriptome',outdir('transcripts.fa'),kind='input')
#    job.add('--output',outdir('hla_2.fa'),kind='output')
#    job.run()

#    job.add('concatenate.py',kind='program')
#    job.add('',outdir('hla.fa'),kind='input')
#    job.add('',outdir('hla_2.fa'),kind='input')
#    job.add('',outdir('hla_12.fa'),kind='output')
#    job.run()

#    job.add('sed',kind='program') # replace blanks with underscores
#    job.add("'s/ /_/g'",kind='parameter')
#    job.add('',outdir('hla_12.fa'),kind='input')
#    job.add('>',outdir('hla-noblanks.fa'),kind='output')
#    job.run()

    if 'cosmic' not in skip_database:
        job.add(_FC_+'get_cosmic.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('cosmic.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database cosmic'! This database is optional."))
    elif job.run():
        file(outdir('cosmic.txt'),'w').write('')

    if 'chimerdb2' not in skip_database:
        job.add(_FC_+'get_chimerdb2.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('chimerdb2.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database chimerdb2'! This database is optional."))
    elif job.run():
        file(outdir('chimerdb2.txt'),'w').write('')

    if 'chimerdb3' not in skip_database:
        job.add(_FC_+'get_chimerdb3.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('chimerdb3kb.txt'),kind='output',command_line='no')
        job.add('',outdir('chimerdb3pub.txt'),kind='output',command_line='no')
        job.add('',outdir('chimerdb3seq.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database chimerdb3'! This database is optional."))
    elif job.run():
        file(outdir('chimerdb3kb.txt'),'w').write('')
        file(outdir('chimerdb3pub.txt'),'w').write('')
        file(outdir('chimerdb3seq.txt'),'w').write('')


    if 'ticdb' not in skip_database:
        job.add(_FC_+'get_ticdb.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('ticdb.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database ticdb'! This database is optional."))
    elif job.run():
        file(outdir('ticdb.txt'),'w').write('')

    if 'conjoing' not in skip_database:
        job.add(_FC_+'get_conjoing.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('conjoing.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database conjoing'! This database is optional."))

        if options.organism == 'homo_sapiens' and (not empty(outdir('conjoing.txt'))):
            job.add('grep',kind='program')
            job.add('-v',kind='parameter')
            job.add("'ENSG00000047932\tENSG00000047936'",kind='parameter') # remove GOPC--ROS1
            job.add('',outdir('conjoing.txt'),kind='input')
            job.add('|',kind='parameter')
            job.add('grep',kind='parameter')
            job.add('-v',kind='parameter')
            job.add("'ENSG00000134853\tENSG00000145216'",kind='parameter') # remove PDGFRA--FIP1L1
            job.add('>',outdir('conjoing_.txt'),kind='output')
            job.run()

            job.add('mv',kind='program')
            job.add('',outdir('conjoing_.txt'),kind='output')
            job.add('',outdir('conjoing.txt'),kind='output')
            job.run()


    elif job.run():
        file(outdir('conjoing.txt'),'w').write('')

    if 'cgp' not in skip_database:
        job.add(_FC_+'get_cgp.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('cgp.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database cgp'! This database is optional."))
    elif job.run():
        file(outdir('cgp.txt'),'w').write('')

    if 'cacg' not in skip_database:
        job.add(_FC_+'get_cacg.py',kind='program')
        job.add('--organism',options.organism,kind='parameter')
        job.add('--output',out_dir,kind='output',checksum='no')
        job.add('',outdir('cacg.txt'),kind='output',command_line='no')
        job.run(error_message = ("If this steps fails to run for whatever reason "+
                "then it can be skipped by re-running fusioncatcher-build with "+
                "command line option '--skip-database cacg'! This database is optional."))
    elif job.run():
        file(outdir('cacg.txt'),'w').write('')

    job.add(_FC_+'get_dgd.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('dgd.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_mitelman.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('mitelman.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_gtex.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('gtex.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_oncokb.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('oncokb.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'get_rtcircrna.py',kind='program')
    job.add('--organism',options.organism,kind='parameter')
    job.add('--output',out_dir,kind='output',checksum='no')
    job.add('',outdir('rtcircrnas.txt'),kind='output',command_line='no')
    job.run()

    job.add(_FC_+'concatenate.py',kind='program')
    job.add('',outdir('trna.fa'),kind='input')
    job.add('',outdir('rrna.fa'),kind='input')
    job.add('',outdir('rrna_unit.fa'),kind='input')
    job.add('',outdir('mt.fa'),kind='input')
    job.add('',outdir('mtrna.fa'),kind='input')
    job.add('',outdir('rtrna_mt.fa'),kind='output')
    job.run()

    job.add(_FC_+'concatenate.py',kind='program')
    job.add('',outdir('trna.fa'),kind='input')
    job.add('',outdir('rrna.fa'),kind='input')
    job.add('',outdir('rrna_unit.fa'),kind='input')
    job.add('',outdir('mt.fa'),kind='input')
    job.add('',outdir('mtrna.fa'),kind='input')
    job.add('',outdir('hla2.fa'),kind='input')
    job.add('',outdir('rtrna_hla_mt.fa'),kind='output')
    job.run()

    job.add(_FC_+'concatenate.py',kind='program')
    job.add('',outdir('trna.fa'),kind='input')
    job.add('',outdir('rrna.fa'),kind='input')
    job.add('',outdir('rrna_unit.fa'),kind='input')
    job.add('',outdir('rtrna.fa'),kind='output')
    job.run()

    job.add(_FC_+'generate_labels_descriptions.py',kind='program')
    job.add('--output',out_dir,kind='output')
    job.add('',outdir('final-list_candidate-fusion-genes.caption.md.txt'),kind='output',command_line='no')
    job.run()


    job.add(_BE_+'bowtie-build',kind='program')
    if bowtie121:
        job.add('--threads',options.processes,kind='parameter')    
    job.add('-f',kind='parameter')
#    job.add('--ntoa',kind='parameter')
    job.add('--quiet',kind='parameter')
    job.add('--offrate','1',kind='parameter')
    job.add('--ftabchars','7',kind='parameter')
    job.add('',outdir('rtrna_hla_mt.fa'),kind='input')
    job.add('',outdir('rtrna_hla_mt_index/'),kind='output')
    job.run(error_message = bowtie_error)

    job.add(_BE_+'bowtie-build',kind='program')
    if bowtie121:
        job.add('--threads',options.processes,kind='parameter')
    job.add('-f',kind='parameter')
#    job.add('--ntoa',kind='parameter')
    job.add('--quiet',kind='parameter')
    job.add('--offrate','1',kind='parameter')
    job.add('--ftabchars','7',kind='parameter')
    job.add('',outdir('rtrna_mt.fa'),kind='input')
    job.add('',outdir('rtrna_mt_index/'),kind='output')
    job.run(error_message = bowtie_error)



    job.add(_BE_+'bowtie-build',kind='program')
    if bowtie121:
        job.add('--threads',options.processes,kind='parameter')
    job.add('-f',kind='parameter')
#    job.add('--ntoa',kind='parameter')
    job.add('--quiet',kind='parameter')
    job.add('--offrate','1',kind='parameter')
    job.add('--ftabchars','7',kind='parameter')
    job.add('',outdir('rtrna.fa'),kind='input')
    job.add('',outdir('rtrna_index/'),kind='output')
    job.run(error_message = bowtie_error)



#    job.add(_BE_+'bowtie-build',kind='program')
#    job.add('-f',kind='parameter')
#    job.add('--quiet',kind='parameter')
##    job.add('--ntoa',kind='parameter')
#    job.add('--offrate','1',kind='parameter')
#    job.add('--ftabchars','7',kind='parameter')
#    job.add('',outdir('hla-noblanks.fa'),kind='input')
#    job.add('',outdir('hla_index/'),kind='output')
#    job.run(error_message = bowtie_error)


#    job.add(_BE_+'bowtie-build',kind='program')
#    job.add('-f',kind='parameter')
#    job.add('--quiet',kind='parameter')
##    job.add('--ntoa',kind='parameter')
#    job.add('--offrate','1',kind='parameter')
#    job.add('--ftabchars','7',kind='parameter')
#    job.add('',outdir('phix174.fa'),kind='input')
#    job.add('',outdir('phix174_index/'),kind='output')
#    job.run(error_message = bowtie_error)

    job.add(_BE_+'bowtie-build',kind='program')
    if bowtie121:
        job.add('--threads',options.processes,kind='parameter')
    job.add('-f',kind='parameter')
    job.add('--quiet',kind='parameter')
#    job.add('--ntoa',kind='parameter')
    job.add('--offrate','1',kind='parameter')
    job.add('--ftabchars','7',kind='parameter')
    job.add('',outdir('transcripts.fa'),kind='input')
    job.add('',outdir('transcripts_index/'),kind='output')
    job.run(error_message = bowtie_error)

#    job.add(_BE_+'bowtie-build',kind='program')
#    if bowtie121:
#        job.add('--threads',options.processes,kind='parameter')
#    job.add('-f',kind='parameter')
#    job.add('--quiet',kind='parameter')
#    job.add('--offrate','1',kind='parameter')
#    job.add('--ftabchars','7',kind='parameter')
#    job.add('',outdir('genome.fa'),kind='input')
#    job.add('',outdir('genome_index/'),kind='output')
#    job.run(error_message = bowtie_error)

    job.add(_B2_+'bowtie2-build',kind='program')
    job.add('-f',kind='parameter')
    job.add('--threads',options.processes,kind='parameter')
    job.add('--quiet',kind='parameter')
    job.add('--offrate','1',kind='parameter')
    job.add('--ftabchars','7',kind='parameter')
    job.add('',outdir('genome.fa'),kind='input')
    job.add('',outdir('genome_index2/index'),kind='output',command_line='no')
    job.add('',outdir('genome_index2/index'),kind='output',checksum='no')
    job.run(error_message = bowtie2_error)


    job.add(_BE_+'bowtie-build',kind='program')
    if bowtie121:
        job.add('--threads',options.processes,kind='parameter')
    job.add('-f',kind='parameter')
    job.add('--quiet',kind='parameter')
    job.add('--offrate','1',kind='parameter')
    job.add('--ftabchars','7',kind='parameter')
    job.add('',outdir('viruses-noblanks.fa'),kind='input')
    job.add('',outdir('viruses_index/'),kind='output')
    job.run(error_message = bowtie_error)

    if not options.skip_blat:
        job.add(_FT_+'faToTwoBit',kind='program')
        job.add('',outdir('genome.fa'),kind='input')
        job.add('',outdir('genome.2bit'),kind='output')
        job.add('-noMask',kind='parameter')
        job.run(error_message = ("Please, check if BLAT (from "+
            "<http://users.soe.ucsc.edu/~kent/src/> and "+
            "<http://hgdownload.cse.ucsc.edu/admin/exe/>) is installed correctly and it "+
            "is in the corresponding PATH (or if 'configuration.cfg' file is "+
            "set up correctly)!\n If there is no wish to use BLAT aligner then please "+
            "(re)run FusionCatcher using command line option '--skip-blat'.\n"+
            "Please, also read its commercial "+
            "license <http://www.kentinformatics.com/> if this applies in your case!"))

    job.clean(outdir('genome.fa'))
    job.clean(outdir('rtrna_mt.fa'))
    job.clean(outdir('rtrna.fa'))
    job.clean(outdir('trna.fa'))
    job.clean(outdir('rrna.fa'))
    job.clean(outdir('mtrna.fa'))
    job.clean(outdir('rrna_unit.fa'))
    job.clean(outdir('mt.fa'))
    job.clean(outdir('hla.fa'))
    job.clean(outdir('hla-noblanks.fa'))
    job.clean(outdir('hla_2.fa'))
#    job.clean(outdir('hla_12.fa'))
    job.clean(outdir('hla_1.txt'))
    job.clean(outdir('hla_2.txt'))
    job.clean(outdir('hla_12.txt'))
    job.clean(outdir('phix174.fa'))
    job.clean(outdir('viruses.fa'))
    job.clean(outdir('viruses-noblanks.fa'))
    job.clean(outdir('ucsc-temp_fully_overlapping_genes.txt'))
    job.clean(outdir('ucsc-temp_partially_overlapping_genes.txt'))
    job.clean(outdir('ucsc-temp_same_strand_overlapping_genes.txt'))
    job.clean(outdir('refseq-temp_fully_overlapping_genes.txt'))
    job.clean(outdir('refseq-temp_partially_overlapping_genes.txt'))
    job.clean(outdir('refseq-temp_same_strand_overlapping_genes.txt'))
    job.clean(outdir('gencode-temp_fully_overlapping_genes.txt'))
    job.clean(outdir('gencode-temp_partially_overlapping_genes.txt'))
    job.clean(outdir('gencode-temp_same_strand_overlapping_genes.txt'))

    job.close()

    ###########
    d = []
    if os.path.isfile(outdir('version.txt')):
        d = [line.rstrip('\r\n') for line in file(outdir('version.txt'),'r') if line.rstrip('\r\n')]
    d.insert(0,parser.get_version())
    file(outdir('version.txt'),'w').writelines([line+'\n' for line in d])
    for line in d:
        print line
    

    #
