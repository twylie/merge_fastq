# Project     : merge_fastq
# File Name   : merge_fastq.py
# Description : Main script for merging FASTQ files by project.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Fri Sep 13 15:35:31 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

"""Merge FASTQ Files

This script will take samplemap information and merge multi-lane FASTQ
files that were split during GTAC@MGI sequencing. New FASTQ files will
be created and written to disk, at the sample level.
"""

# import mergefastq
import argparse
import sys


def collect_cli_arguments(version: int) -> argparse.Namespace:
    """Collect the command line arguments."""
    # TODO: Add argparse command line options here...
    return


if __name__ == '__main__':
    VERSION = '0.0.1'

    if not sys.version_info >= (3, 10):
        raise OSError(
            'Python version must be 3.10 or greater.',
            sys.version_info
        )

    args = collect_cli_arguments(version=VERSION)

# __END__
