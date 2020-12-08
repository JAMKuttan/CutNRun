#!/usr/bin/env python3

#
# * --------------------------------------------------------------------------
# * Licensed under MIT (https://github.com/JAMKuttan/CutNRun/LICENSE.md)
# * --------------------------------------------------------------------------
#

'''Remove duplicates and filter unmapped reads.'''

import os
import subprocess
import argparse
import shutil
import shlex
import logging
from multiprocessing import cpu_count
import utils
import pandas as pd

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
STRIP_EXTENSIONS = ['.bam', '.srt']


def get_args():
    '''Define arguments.'''

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--bam',
                        help="The bam file to run filtering and qc on.",
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


def filter_mapped_pe(bam, bam_basename):
    '''Use samtools to filter unmapped reads for PE data.'''

    filt_bam_prefix = bam_basename + ".filt.srt"
    filt_bam_filename = filt_bam_prefix + ".bam"

    out, err = utils.run_pipe([
        "samtools view -F 1804 -f 3 -q 30 -hu %s" % (bam),
        "samtools sort -n -@ %d -o %s" % (cpu_count(), filt_bam_filename)])
    if err:
        logger.error("samtools filter error: %s", err)

    filter_index_command = \
        "samtools index -@ %d %s" % (cpu_count(), filt_bam_filename)
    logger.info(filter_index_command)
    subprocess.check_output(shlex.split(filter_index_command))

    return filt_bam_filename


def filter_mapped_se(bam, bam_basename):
    '''Use samtools to filter unmapped reads for SE data.'''

    filt_bam_prefix = bam_basename + ".filt.srt"
    filt_bam_filename = filt_bam_prefix + ".bam"

    out, err = utils.run_pipe([
        "samtools view -F 1804 -q 30 -hu %s" % (bam),
        "samtools sort -n -@ %d -o %s" % (cpu_count(), filt_bam_filename)])
    if err:
        logger.error("samtools filter error: %s", err)

    filter_index_command = \
        "samtools index -@ %d %s" % (cpu_count(), filt_bam_filename)
    logger.info(filter_index_command)
    subprocess.check_output(shlex.split(filter_index_command))

    return filt_bam_filename


def dedup_mapped(bam, bam_basename, paired):
    '''Use picard and samtools to remove duplicates.'''

    # Markduplicates
    dup_filename = bam_basename + ".dups.txt"
    dedup_filename = bam_basename + ".dedup.bam"
    dedup_sorted_filename = bam_basename + ".dedup.srt.bam"
    picard_params = "REMOVE_DUPLICATES=true"

    out, err = utils.run_pipe([
        "java -jar $PICARDJAR MarkDuplicates I=%s O=%s M=%s %s" % (bam, dedup_filename, dup_filename, picard_params),
        "samtools sort -n -@ %d -o %s" % (cpu_count(), dedup_sorted_filename)
        "samtools index %s" % (dedup_sorted_filename)])

    # Remove duplicates
    final_bam_prefix = bam_basename + ".dedup"
    final_bam_filename = final_bam_prefix + ".bam"

    if paired:  # paired-end data
        samtools_dedupe_command = \
            "samtools view -F 1804 -f 2 -b %s" % (tmp_dup_mark_filename)
    else:
        samtools_dedupe_command = \
            "samtools view -F 1804 -b %s" % (tmp_dup_mark_filename)

    with open(final_bam_filename, 'w') as temp_file:
        logger.info(samtools_dedupe_command)
        subprocess.check_call(
            shlex.split(samtools_dedupe_command),
            stdout=temp_file)

    # Index final bam file
    sambamba_index_command = \
        "samtools index -@ %d %s" % (cpu_count(), final_bam_filename)
    logger.info(sambamba_index_command)
    subprocess.check_output(shlex.split(sambamba_index_command))

    # Generate mapping statistics
    mapstats_filename = final_bam_prefix + ".flagstat.qc"
    with open(mapstats_filename, 'w') as temp_file:
        flagstat_command = "sambamba flagstat -t %d %s" \
                            % (cpu_count(), final_bam_filename)
        logger.info(flagstat_command)
        subprocess.check_call(shlex.split(flagstat_command), stdout=temp_file)

    os.remove(bam)
    return tmp_dup_mark_filename


def compute_complexity(bam, paired, bam_basename):
    '''Calculate library complexity .'''

    pbc_file_qc_filename = bam_basename + ".pbc.qc"
    tmp_pbc_file_qc_filename = "tmp.%s" % (pbc_file_qc_filename)

    # Sort by name
    # convert to bedPE and obtain fragment coordinates
    # sort by position and strand
    # Obtain unique count statistics

    # PBC File output
    # Sample Name[tab]
    # TotalReadPairs [tab]
    # DistinctReadPairs [tab]
    # OneReadPair [tab]
    # TwoReadPairs [tab]
    # NRF=Distinct/Total [tab]
    # PBC1=OnePair/Distinct [tab]
    # PBC2=OnePair/TwoPair
    pbc_headers = ['TotalReadPairs',
                   'DistinctReadPairs',
                   'OneReadPair',
                   'TwoReadPairs',
                   'NRF',
                   'PBC1',
                   'PBC2']

    if paired:
        steps = [
            "samtools sort -@%d -n %s" % (cpu_count(), bam),
            "bamToBed -bedpe -i stdin",
            r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}'"""]
    else:
        steps = [
            "bamToBed -i %s" % (bam),
            r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}'"""]
    steps.extend([
        "grep -v 'chrM'",
        "sort",
        "uniq -c",
        r"""awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'"""
        ])
    out, err = utils.run_pipe(steps, tmp_pbc_file_qc_filename)
    if err:
        logger.error("PBC file error: %s", err)

    # Add Sample Name and headers
    pbc_file = pd.read_csv(tmp_pbc_file_qc_filename, sep='\t', header=None,
                           names=pbc_headers)
    pbc_file['Sample'] = bam_basename
    pbc_headers_new = list(pbc_file)
    pbc_headers_new.insert(0, pbc_headers_new.pop(pbc_headers_new.index('Sample')))
    pbc_file = pbc_file[pbc_headers_new]
    pbc_file.to_csv(pbc_file_qc_filename, header=True, sep='\t', index=False)
    os.remove(bam)
    os.remove(bam + '.bai')
    os.remove(tmp_pbc_file_qc_filename)


def main():
    args = get_args()
    paired = args.paired
    bam = args.bam

    # Create a file handler
    handler = logging.FileHandler('map_qc.log')
    logger.addHandler(handler)

    # Check if tools are present
    check_tools()

    # Run filtering for either PE or SE
    bam_basename = os.path.basename(
        utils.strip_extensions(bam, STRIP_EXTENSIONS))
    if paired:  # paired-end data
        filter_bam_filename = filter_mapped_pe(bam, bam_basename)

    else:
        filter_bam_filename = filter_mapped_se(bam, bam_basename)

    # Remove duplicates
    markdup_bam_filename = dedup_mapped(filter_bam_filename, bam_basename, paired)

    # Compute library complexity
    compute_complexity(markdup_bam_filename, paired, bam_basename)


if __name__ == '__main__':
    main()

