# Project     : merge_fastq
# File Name   : rename_samples.py
# Description : Functions for renaming FASTQ samples.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Tue Sep 17 14:30:14 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from typing_extensions import Self
from pandas import DataFrame  # type: ignore


class RenameSamples:
    """A class for renaming FASTQ sample ids.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    args.rename : STR
        The path to the --rename file used to map samplemap sample names
        to revised/renamed sample names.
    """

    def __init__(self: Self, args: argparse.Namespace) -> None:
        """Construct the class."""
        self.args = args
        self.__populate_df()
        return

    def __populate_df(self: Self) -> None:
        """Load the dataframe and place in the object.

        We are loading the --rename input file into a dataframe within
        an object instance. The rename file is tab-delimited and
        contains the following fields:
            1. samplemap_sample_id
            2. revised_sample_id
            3. comments

        Raises
        ------
        ValueError : Input file column names are incorrect.

        ValueError : Input file contains no values.

        ValueError : The samplemap_sample_id column values are not
                     unique.

        ValueError : The revised_sample_id column values are not unique.
        """
        self.df = pd.read_csv(self.args.rename, sep='\t')

        fields = ['samplemap_sample_id', 'revised_sample_id', 'comments']
        if list(self.df.columns) != fields:
            raise ValueError(
                'Input file column names are incorrect.',
                self.df.columns
            )

        if len(self.df.values) == 0:
            raise ValueError(
                'Input file contains no values.',
                len(self.df.values)
            )

        samplemap_sample_id_count = len(self.df['samplemap_sample_id'].values)
        samplemap_sample_id_ucount = len(
            set(self.df['samplemap_sample_id'].values)
        )
        if samplemap_sample_id_count != samplemap_sample_id_ucount:
            raise ValueError('The samplemap_sample_id column values '
                             'are not unique.')

        revised_sample_id_count = len(self.df['revised_sample_id'].values)
        revised_sample_id_ucount = len(
            set(self.df['revised_sample_id'].values)
        )
        if revised_sample_id_count != revised_sample_id_ucount:
            raise ValueError('The revised_sample_id column values '
                             'are not unique.')
        return

    def copy_df(self: Self) -> DataFrame:
        """Return a copy of the dataframe from the RenameSamples class."""
        return self.df.copy()

# __END__
