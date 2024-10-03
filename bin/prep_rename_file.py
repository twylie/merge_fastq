#!/usr/bin/python3.10

# Project     : merge_fastq
# File Name   : prep_rename_file.py
# Description : Prepare a rename file for merging FASTQ files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Thu Oct 03 10:41:46 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import sys
from pathlib import Path
import pandas as pd  # type: ignore

"""Prep Rename File

This script helps to prepare a rename file for merging FASTQ files using
the merge_fastq command. The rename file is tab-delimited and contains
the following fields:

    1. samplemap_sample_id
    2. revised_sample_id
    3. comments

This script will populate both columns 1 an 2 with the original sample
names found in the user supplied Samplemap.csv files. It is up to the
user to manually edit the rev_sample_name field, when appropriate.

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
        description=('Prepare a rename file for merging FASTQ files '
                     'with merge_fastq.'),
        prog='prep_rename_file',
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
        help='List of samplemap files to be evaluated.',
        required=True,
        nargs='+'
    )

    required_group.add_argument(
        '--rename-out',
        metavar='FILE',
        action='store',
        help='File path to write the rename file.',
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
    if Path(args.rename_out).exists() is True:
        raise FileExistsError(
            'The --rename-out file already exists.',
            args.rename_out
        )
    for smap_file in args.samplemap:
        if Path(smap_file).is_file() is False:
            raise FileNotFoundError(
                'A --samplemap input file does not exist.',
                smap_file
            )
    return


# MAIN ########################################################################

if __name__ == '__main__':
    VERSION = '1.0.0-alpha'

    if not sys.version_info >= (3, 10):
        raise OSError(
            'Python version must be 3.10 or greater.',
            sys.version_info
        )

    # Collect the command line arguments.

    args = collect_cli_arguments(version=VERSION)
    eval_cli_arguments(args=args)

    # Concatenate the list of Samplemap.csv files and get unique sample
    # names.

    dfs: list = list()
    for smap_path in args.samplemap:
        df = pd.read_csv(smap_path)
        dfs.append(df)
    df_smaps = pd.concat(dfs).reset_index(drop=True)
    sample_names = list(df_smaps['Library Name'].unique())

    # Write the prepared rename file.

    header_row = '\t'.join(
        ['samplemap_sample_id', 'revised_sample_id', 'comments']
    )
    with open(args.rename_out, 'w') as fho:
        fho.write(header_row + '\n')
        for sample_name in sorted(sample_names):
            row = '\t'.join([sample_name, sample_name, ''])
            fho.write(row + '\n')

# __END__
