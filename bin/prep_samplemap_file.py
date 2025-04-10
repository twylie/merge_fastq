#!/usr/bin/python3.10

# Project     : merge_fastq
# File Name   : prep_samplemap_file.py
# Description : Prepare a samplemap file for merging FASTQ files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Thu Apr 10 11:57:45 CDT 2025
# Copyright   : Copyright (C) 2024-2025 by T.N. Wylie. All rights reserved.

import argparse
import sys
from pathlib import Path
import pandas as pd  # type: ignore

"""Prep Samplemap File

This script helps to prepare a Samplemap file for merging FASTQ files
using the merge_fastq command. At the moment, this script's function is
to convert all spaces in the Library Name field to underscores and write
the revised file to disk. This script may be run prior to running the
merge_fastq code. The output file will be a comma-delimited
Samplemap.csv file.

See Also
--------
mergefastq.MergeFastq
"""


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
        description=('Prepare a samplemap file for merging FASTQ files '
                     'with merge_fastq.'),
        prog='prep_samplemap_file',
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
        '--samplemap',
        metavar='FILE',
        action='store',
        help='Samplemap file to be evaluated.',
        required=True
    )

    required_group.add_argument(
        '--out',
        metavar='FILE',
        action='store',
        help='File path to write the revised samplemap file.',
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
    FileExistsError
        The --rename-out file already exists.

    FileNotFoundError
        A --samplemap input file does not exist.

    Returns
    -------
    None
    """
    if Path(args.out).exists() is True:
        raise FileExistsError(
            'The --out file already exists.',
            args.out
        )
    if Path(args.samplemap).is_file() is False:
        raise FileNotFoundError(
            'A --samplemap input file does not exist.',
            args.samplemap
        )
    return


# MAIN ########################################################################

if __name__ == '__main__':
    VERSION = '1.0.0'

    if not sys.version_info >= (3, 10):
        raise OSError(
            'Python version must be 3.10 or greater.',
            sys.version_info
        )

    # Collect the command line arguments.

    args = collect_cli_arguments(version=VERSION)
    eval_cli_arguments(args=args)

    # We will be converting all of spaces to underscores for the values
    # in the Library Name column. The resultant samplemap file should
    # have no spaces in the source library names.

    df = pd.read_csv(args.samplemap)
    pre_names = list(df['Library Name'])
    df['Library Name'] = df['Library Name'].replace('\\s+', '_', regex=True)
    post_names = list(df['Library Name'])

    for i, pairs in enumerate(zip(pre_names, post_names)):
        pre, post = pairs
        if post != pre:
            print(f'{i}. "{pre}" --> "{post}"')

    df.to_csv(args.out, index=False)
    print()
    print(f'wrote: {Path(args.out).as_posix()}')

# __END__
