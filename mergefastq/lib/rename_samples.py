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
from pathlib import Path


class RenameSamples:
    """A class for renaming FASTQ sample ids.

    This class is used to import and interact with a --rename file. The
    rename file is used to map old sample names to new ones, when the
    associated FASTQ files are merged--see Samplemap and MergeFastq
    classes for details on merging.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    Attributes
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    df : DataFrame
        A dataframe representation of an imported --rename file.

    Methods
    -------
    copy_df() -> DataFrame
        Return a copy of the dataframe from the RenameSamples class.

    copy_rename_file(outdir)
        Copy the rename file to another directory.

    Examples
    --------
    rename_samples = RenameSamples(args=args)
    df = rename_samples.copy_df()
    df.copy_rename_file(outdir=outdir)
    """

    def __init__(self: Self, args: argparse.Namespace) -> None:
        """Construct the class.

        Parameters
        ----------
        args : argparse.Namespace
            Arguments for merging FASTQ files as provided by argparse.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        self.args = args
        self.df: DataFrame = DataFrame()
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

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Input file column names are incorrect.

        ValueError
            Input file contains no values.

        ValueError
            The samplemap_sample_id column values are not unique.

        ValueError
            The revised_sample_id column values are not unique.

        ValueError
            Samplemap to revised id is not one-to-one.

        Returns
        -------
        None
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

        dfg = self.df.groupby('samplemap_sample_id')
        for sample_id in dfg.groups:
            sample_id_count = len(dfg.get_group(sample_id))
            if sample_id_count > 1:
                raise ValueError('Samplemap to revised id is not one-to-one.')

        dfg = self.df.groupby('revised_sample_id')
        for sample_id in dfg.groups:
            sample_id_count = len(dfg.get_group(sample_id))
            if sample_id_count > 1:
                raise ValueError('Samplemap to revised id is not one-to-one.')
        return

    def copy_df(self: Self) -> DataFrame:
        """Return a copy of the dataframe from the RenameSamples class.

        Parameters
        ----------
        None

        Raises
        ------
        None

        Returns
        -------
        DataFrame
            Returns a copy of the RenameSamples dataframe.
        """
        return self.df.copy()

    def copy_rename_file(self: Self, outdir: str) -> None:
        """Copy the rename file to another directory.

        Parameters
        ----------
        outdir : str
            A qualified path to copy the source --rename file.

        Raises
        ------
        IsADirectoryError
            Outdir directory does not exist.

        Returns
        -------
        None
        """
        if Path(outdir).is_dir() is False:
            raise IsADirectoryError(
                'Outdir directory does not exist.',
                outdir
            )
        else:
            src = Path(self.args.rename)
            dest = Path(outdir) / 'rename.tsv'
            dest.write_text(src.read_text())
        return

# __END__
