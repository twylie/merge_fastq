#!/usr/bin/python3.10

# Project     : merge_fastq
# File Name   : eval_fastq_counts.py
# Description : Evaluate and compare read counts for merged FASTQ files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Tue Oct 01 14:08:59 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

"""Eval FASTQ Counts

This script will evaluate and compare read counts for merged FASTQ
files, specifically the read count outputs as provided by the
merge_fastq.py script.

This script should be executed after running the merge_fastq.py script
for a set of FASTQ files.
"""

import mergefastq  # type: ignore
import argparse
import sys
from pathlib import Path


# FUNCTIONS ###################################################################

def collect_cli_arguments(version: str) -> argparse.Namespace:
    """Collect the command line arguments.

    Parameters
    ----------
    version : str
        A version id in semantic version format.

    Raises
    ------
    None

    Returns
    -------
    argparse.Namespace
        Returns an argparse object with argument information.
    """
    parser = argparse.ArgumentParser(
        description=('Evaluate and compare read counts for merged '
                     'FASTQ files.'),
        prog='eval_fastq_counts',
        add_help=False
    )

    # Optional arguments.

    parser.add_argument(
        '-h',
        '--help',
        action='help',
        help='Display the extended usage statement.'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=version,
        help='Display the software version number.'
    )

    parser.set_defaults(lsf_dry=True)

    # Required arguments.

    required_group = parser.add_argument_group('required')

    required_group.add_argument(
        '--merged-samplemap',
        metavar='FILE',
        action='store',
        help='A merged_samplemap.tsv file.',
        required=True
    )

    required_group.add_argument(
        '--gtac-counts',
        metavar='FILE',
        action='store',
        help='A gtac_read_counts.tsv file.',
        required=True
    )

    required_group.add_argument(
        '--outdir',
        metavar='DIR',
        action='store',
        help='Output directory path to write results.',
        required=True
    )
    return parser.parse_args()


def eval_cli_arguments(args: argparse.Namespace) -> None:
    """Evaluate the individual command line arguments.

    We can perform a few simple evaluations on reagent input files and
    directories prior to moving forward with downstream functions.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments object as provided by argparse.

    Raises
    ------
    FileNotFoundError
        The --merged-samplemap input file does not exist.

    FileNotFoundError
        The --gtac-counts input file does not exist.

    NotADirectoryError
        The --outdir does not exist.

    Returns
    -------
    None
    """
    if Path(args.merged_samplemap).is_file() is False:
        raise FileNotFoundError(
            'The --merged-samplemap input file does not exist.',
            args.merged_samplemap
        )
    if Path(args.gtac_counts).is_file() is False:
        raise FileNotFoundError(
            'The --gtac-counts input file does not exist.',
            args.merged_samplemap
        )
    if Path(args.outdir).is_dir() is False:
        raise NotADirectoryError(
            'The --outdir does not exist.',
            args.outdir
        )
    return


# MAIN ########################################################################

if __name__ == '__main__':
    VERSION = '0.0.11-alpha'

    if not sys.version_info >= (3, 10):
        raise OSError(
            'Python version must be 3.10 or greater.',
            sys.version_info
        )

    # Collect the command line arguments.

    args = collect_cli_arguments(version=VERSION)
    eval_cli_arguments(args=args)

    # Parse the input reagent files and create objects for downstream
    # processing.

    read_counts = mergefastq.ReadCountsSource(args=args)
    read_counts.tell_comp_differences()

# __END__
