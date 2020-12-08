#!/usr/bin/env python3

#
# * --------------------------------------------------------------------------
# * Licensed under MIT (https://github.com/JAMKuttan/CutNRun/LICENSE.md)
# * --------------------------------------------------------------------------
#

'''Align reads to reference genome.'''

import os
import subprocess
import argparse
import shutil
import shlex
import logging
from multiprocessing import cpu_count
import utils

EPILOG = '''
For more details:
        %(prog)s --help
'''

# SETTINGS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

# the order of this list is important.
# strip_extensions strips from the right inward, so
# the expected right-most extensions should appear first (like .gz)
# Modified from J. Seth Strattan
STRIP_EXTENSIONS = ['.gz', '.fq', '.fastq', '_trimmed']


def get_args():
    '''Define arguments.'''

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-f', '--fastq',
                        help="The fastq file to run triming on.",
                        nargs='+',
                        required=True)

    parser.add_argument('-r', '--reference',
                        help="The bwa index of the reference genome.",
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


# Functions


def check_tools():
    '''Checks for required componenets on user system'''

    logger.info('Checking for required libraries and components on this system')

    bowtie2_path = shutil.which("bowtie2")
    if bowtie2_path:
        logger.info('Found bowtie2: %s', bowtie2_path)

        # Get Version
        bowtie2_version_command = "bowtie2"
        try:
            subprocess.check_output(bowtie2_version_command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            bowtie2_version = e.output

        # Write to file
        bowtie2_file = open("version_bowtie2.txt", "wb")
        bowtie2_file.write(bowtie2_version)
        bowtie2_file.close()
    else:
        logger.error('Missing bowtie2')
        raise Exception('Missing bowtie2')

    samtools_path = shutil.which("samtools")
    if samtools_path:
        logger.info('Found samtools: %s', samtools_path)

        # Get Version
        samtools_version_command = "samtools --version"
        samtools_version = subprocess.check_output(samtools_version_command, shell=True)

        # Write to file
        samtools_file = open("version_samtools.txt", "wb")
        samtools_file.write(samtools_version)
        samtools_file.close()
    else:
        logger.error('Missing samtools')
        raise Exception('Missing samtools')


def align_se(fastq, reference, fastq_basename):
    '''Use bowtie2 to align SE data.'''

    bam_filename = '%s.bam' % (fastq_basename)

    steps = [
        "bowtie2 -p %d --very-sensitive -x %s -1 %s"
        % (cpu_count(), reference, fastq[0])
        "samtools view -@%d -Su -" % (cpu_count()),
        "samtools sort -@%d -o %s"
        % (cpu_count(), bam_filename)]

    out, err = utils.run_pipe(steps)
    if err:
        logger.error("samse/samtools error: %s", err)

    return bam_filename


def align_pe(fastq, reference, fastq_basename):
    '''Use bowtie2 to align PE data.'''

    bam_filename = '%s.bam' % (fastq_basename)

    steps = [
        "bowtie2 -p %d --very-sensitive -x %s -1 %s -2 %s"
        % (cpu_count(), reference, fastq[0], fastq[1])
        "samtools view -@%d -Su -" % (cpu_count()),
        "samtools sort -@%d -o %s"
        % (cpu_count(), bam_filename)]

    out, err = utils.run_pipe(steps)
    if err:
        logger.error("samtools error: %s", err)

    return bam_filename


def main():
    args = get_args()
    paired = args.paired
    fastq = args.fastq
    reference = args.reference
    sample = args.sample

    # Create a file handler
    handler = logging.FileHandler('map.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Make file basename
    fastq_basename = sample

    # Run alignment for either PE or SE
    if paired:  # paired-end data
        bam_filename = align_pe(fastq, reference, fastq_basename)

    else:
        bam_filename = align_se(fastq, reference, fastq_basename)

    bam_mapstats_filename = '%s.flagstat.qc' % (fastq_basename)
    with open(bam_mapstats_filename, 'w') as temp_file:
        subprocess.check_call(
            shlex.split("samtools flagstat %s" % (bam_filename)),
            stdout=temp_file)

    #Genome/Bad fastq File Check
    file_check = open(bam_mapstats_filename).readlines()
    percent = file_check[4].split('(')[1]
    percent = percent.split('%')[0]
    if float(percent) < 10:
        raise Exception ('Mapped Genes too low: Check for correct Genotype')


if __name__ == '__main__':
    main()
