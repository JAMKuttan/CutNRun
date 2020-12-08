#!/usr/bin/env python3

#
# * --------------------------------------------------------------------------
# * Licensed under MIT (https://github.com/JAMKuttan/CutNRun/LICENSE.md)
# * --------------------------------------------------------------------------
#

'''Trim low quality reads and remove sequences less than 35 base pairs.'''

import subprocess
import argparse
import shutil
import os
import logging

EPILOG = '''
For more details:
        %(prog)s --help
'''

# SETTINGS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)


def get_args():
    '''Define arguments.'''

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-f', '--fastq',
                        help="The fastq file to run triming on.",
                        nargs='+',
                        required=True)

    parser.add_argument('-s', '--sample',
                        help="The name of the sample.",
                        required=True)

    parser.add_argument('-p', '--paired',
                        help="True/False if paired-end or single end.",
                        default=False,
                        action='store_true')

    args = parser.parse_args()
    return args


def check_tools():
    '''Checks for required componenets on user system.'''

    logger.info('Checking for required libraries and components on this system')

    trimgalore_path = shutil.which("trim_galore")
    if trimgalore_path:
        logger.info('Found trimgalore: %s', trimgalore_path)

        # Get Version
        trim_version_command = "trim_galore --version"
        trimgalore_version = subprocess.check_output(trim_version_command, shell=True)

        # Write to file
        trimgalore_file = open("version_trimgalore.txt", "wb")
        trimgalore_file.write(trimgalore_version)
        trimgalore_file.close()
    else:
        logger.error('Missing trimgalore')
        raise Exception('Missing trimgalore')

    cutadapt_path = shutil.which("cutadapt")
    if cutadapt_path:
        logger.info('Found cutadapt: %s', cutadapt_path)

        # Get Version
        cutadapt_version_command = "cutadapt --version"
        cutadapt_version = subprocess.check_output(cutadapt_version_command, shell=True)

        # Write to file
        cutadapt_file = open("version_cutadapt.txt", "wb")
        cutadapt_file.write(b"Version %s" % (cutadapt_version))
        cutadapt_file.close()
    else:
        logger.error('Missing cutadapt')
        raise Exception('Missing cutadapt')


def rename_reads(fastq, sample, paired):
    '''Rename fastq files by sample name.'''

    # Get current directory to build paths
    cwd = os.getcwd()

    renamed_fastq = []

    if paired:  # paired-end data
        # Set file names
        renamed_fastq.append(cwd + '/' + sample + '_R1.fastq.gz')
        renamed_fastq.append(cwd + '/' + sample + '_R2.fastq.gz')

        # Great symbolic links
        os.symlink(fastq[0], renamed_fastq[0])
        os.symlink(fastq[1], renamed_fastq[1])
    else:
        # Set file names
        renamed_fastq.append(cwd + '/' + sample + '_R1.fastq.gz')

        # Great symbolic links
        os.symlink(fastq[0], renamed_fastq[0])

    return renamed_fastq


def trim_reads(fastq, paired):
    '''Run trim_galore on 1 or 2 files.'''

    if paired:  # paired-end data
        trim_params = '--paired -q 25 --illumina --gzip --length 35'
        trim_command = "trim_galore %s %s %s " \
                    % (trim_params, fastq[0], fastq[1])
    else:
        trim_params = '-q 25 --illumina --gzip --length 35'
        trim_command = "trim_galore %s %s " \
                    % (trim_params, fastq[0])

    logger.info("Running trim_galore with %s", trim_command)

    trim = subprocess.Popen(trim_command, shell=True)
    out, err = trim.communicate()


def main():
    args = get_args()
    fastq = args.fastq
    sample = args.sample
    paired = args.paired

    # Create a file handler
    handler = logging.FileHandler('trim.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Rename fastq files by sample
    fastq_rename = rename_reads(fastq, sample, paired)

    # Run trim_reads
    trim_reads(fastq_rename, paired)


if __name__ == '__main__':
    main()
