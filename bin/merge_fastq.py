#!/usr/bin/python3.10

# Project     : merge_fastq
# File Name   : merge_fastq.py
# Description : Main script for merging FASTQ files by project.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Fri Sep 13 15:35:31 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

"""Merge FASTQ Files

This script will take samplemap and project-specific information and
merge multi-lane FASTQ files that were split during GTAC@MGI sequencing.
New FASTQ files will be created and written to disk, at the sample
level.
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
        description=('Merge FASTQ files at sample-level.'),
        prog='merge_fastq',
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

    parser.add_argument(
        '--lsf-image',
        metavar='STR',
        action='store',
        help='Docker image for LSF processing. [twylie/merge_fastq:latest]',
        required=False,
        default='twylie/merge_fastq:latest'
    )

    parser.add_argument(
        '--lsf-group',
        metavar='STR',
        action='store',
        help='Group for LSF processing. [compute-kwylie]',
        required=False,
        default='compute-kwylie'
    )

    parser.add_argument(
        '--lsf-queue',
        metavar='STR',
        action='store',
        help='Queue for LSF processing. [general]',
        required=False,
        default='general'
    )

    parser.add_argument(
        '--lsf-dry',
        action='store_true',
        help='Dry run for LSF processing. (default)',
        required=False
    )

    parser.add_argument(
        '--no-lsf-dry',
        action='store_false',
        dest='lsf_dry',
        help='Executes LSF processing.',
        required=False
    )

    parser.set_defaults(lsf_dry=True)

    # Required arguments.

    required_group = parser.add_argument_group('required')

    required_group.add_argument(
        '--samplemap',
        metavar='FILE',
        action='store',
        help='List of samplemap files to be merged.',
        required=True,
        nargs='+'
    )

    required_group.add_argument(
        '--outdir',
        metavar='DIR',
        action='store',
        help='Output directory path to write results.',
        required=True
    )

    required_group.add_argument(
        '--rename',
        metavar='FILE',
        action='store',
        help='Path to a sample rename file.',
        required=True
    )

    required_group.add_argument(
        '--lsf-vol',
        metavar='PATH',
        action='store',
        help='Top-level Docker volume path(s) for WashU LSF processing.',
        required=True,
        nargs='+'
    )

    project_tags = ('MIDAS', 'PLACENTA', 'PTLD')
    required_group.add_argument(
        '--project',
        action='store',
        help='Project tag/name.',
        choices=project_tags,
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
        The --rename input file does not exist.

    FileNotFoundError
        A --samplemap input file does not exist.

    Returns
    -------
    None
    """
    if Path(args.rename).is_file() is False:
        raise FileNotFoundError(
            'The --rename input file does not exist.',
            args.rename
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
    VERSION = '1.0.14'

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

    rename_samples = mergefastq.RenameSamples(args=args)
    samplemap = mergefastq.Samplemap(args=args, rename=rename_samples)
    merge_fastq = mergefastq.MergeFastq(
        args=args,
        rename=rename_samples,
        samplemap=samplemap
    )

    # Proceed creating FASTQ merge commands mediated by LSF bsub jobs,
    # one job per sample. We will also perform read count evaluations.

    merge_fastq.setup_output_dirs()
    rename_samples.copy_rename_file(outdir=args.outdir)
    samplemap.copy_samplemaps(outdir=args.outdir)
    merged_df = Path(args.outdir) / 'merged_samplemap.tsv'
    merge_fastq.write_df(file_path=str(merged_df.resolve()))
    read_counts = mergefastq.ReadCountsGtac(
        args=args,
        merged_tsv=str(merged_df.resolve())
    )
    read_counts.calc_gtac_read_coverage()
    gtac_read_counts = Path(args.outdir) / 'gtac_read_counts.tsv'
    read_counts.write_df(file_path=str(gtac_read_counts.resolve()))
    merge_fastq.prepare_lsf_cmds()
    merge_fastq.launch_lsf_jobs()

# __END__
