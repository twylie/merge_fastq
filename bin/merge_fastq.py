# Project     : merge_fastq
# File Name   : merge_fastq.py
# Description : Main script for merging FASTQ files by project.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Fri Sep 13 15:35:31 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

"""Merge FASTQ Files

This script will take samplemap, QC, and project-specific information
and merge multi-lane FASTQ files that were split during GTAC@MGI
sequencing. New FASTQ files will be created and written to disk, at the
sample level.
"""

import mergefastq  # type: ignore
import argparse
import sys
from pathlib import Path


def collect_cli_arguments(version: str) -> argparse.Namespace:
    """Collect the command line arguments."""
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
        '--target-rp-count',
        metavar='INT',
        type=int,
        action='store',
        help='Overwrites the target read pair count for project.',
        required=False
    )

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
        '--qc',
        metavar='FILE',
        action='store',
        help='List of QC files to be merged.',
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

    project_tags = ('BVI', 'COLOCARE', 'MIDAS', 'PLACENTA', 'PTLD', 'SHINE')
    required_group.add_argument(
        '--project',
        action='store',
        help='Project tag/name.',
        choices=project_tags,
        required=True
    )
    return parser.parse_args()


def populate_target_rp_counts(args: argparse.Namespace) -> None:
    """Set target read pair counts based on project tag.

    We will set the desired minimal read pair count for a sample based
    on the project. Values are in terms of read pair count per sample.
    These values will be used downstream during sequence throughput
    assessment and reporting.

    If the --target-rp-count argument was manually passed from the
    command line, that value will be used in place of the default
    project-specific values.

    Projects
    --------
    BVI      : Wylie/Stout maternal virome R01 project.
    COLOCARE : L. Ding project.
    MIDAS    : Wylie/Fritz MRSA related R01 project.
    PLACENTA : Wylie/Ernst placental inflammation R01 project.
    PTLD     : Dharnidharka PTLD R01 project.
    SHINE    : Subset of Fritz data.

    Raises
    ------
    ValueError : Project tag/name does not have associated minimal
                 target read pair count.
    """
    # TODO: I need to update the minimal read pair count values for each
    # project when Kristine provides the actual numbers for the #
    # projects.
    if not args.target_rp_count:
        match args.project:
            case 'BVI':
                args.target_rp_count = 50_000_000
            case 'COLOCARE':
                args.target_rp_count = 50_000_000
            case 'MIDAS':
                args.target_rp_count = 40_000_000
            case 'PLACENTA':
                args.target_rp_count = 40_000_000
            case 'PTLD':
                args.target_rp_count = 50_000_000
            case 'SHINE':
                args.target_rp_count = 50_000_000
            case _:
                raise ValueError(
                    ('Project tag/name does not have associated '
                     'minimal target read pair count.'),
                    args.project
                )
    return


def eval_cli_arguments(args: argparse.Namespace) -> None:
    """Evaluate the individual command line arguments.

    We can perform a few simple evaluations on reagent input files and
    directories prior to moving forward with downstream functions.

    Raises
    ------
    FileNotFoundError : The --rename input file does not exist.

    FileNotFoundError : A --qc input file does not exist.

    FileNotFoundError : A --samplemap input file does not exist.

    IsADirectoryError : The --outdir directory already exists.
    """
    if Path(args.rename).is_file() is False:
        raise FileNotFoundError(
            'The --rename input file does not exist.',
            args.rename
        )
    for qc_file in args.qc:
        if Path(qc_file).is_file() is False:
            raise FileNotFoundError(
                'A --qc input file does not exist.',
                qc_file
            )
    for smap_file in args.samplemap:
        if Path(smap_file).is_file() is False:
            raise FileNotFoundError(
                'A --samplemap input file does not exist.',
                smap_file
            )
    if Path(args.outdir).is_dir() is True:
        raise IsADirectoryError(
            'The --outdir directory already exists.',
            args.outdir
        )
    return


if __name__ == '__main__':
    VERSION = '0.0.8'

    if not sys.version_info >= (3, 10):
        raise OSError(
            'Python version must be 3.10 or greater.',
            sys.version_info
        )

    args = collect_cli_arguments(version=VERSION)
    populate_target_rp_counts(args=args)
    eval_cli_arguments(args=args)

    mergefastq.RenameSamples(args=args)

# __END__
