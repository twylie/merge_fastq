# Project     : merge_fastq
# File Name   : read_counts.py
# Description : Methods for sequence throughput evaluation.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Fri Sep 27 12:21:33 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from pandas import DataFrame
from typing_extensions import Self
from pathlib import Path
import hashlib


class ReadCountsGtac:
    """A class for read count evaluation.

    This class provides methods for calculating and evaluating FASTQ
    file read counts. We will handle read counts as provided by GTAC@MGI
    (see Samplemap class).

    Parameters
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    merged_tsv : str
        A qualified path to a merged dataframe file as written by the
        MergeFastq class.

    Attributes
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    df_gtac_seqcov : DataFrame
        A dataframe of evaluated sequence throughput for samples from
        the MergeFastq object. Calculations are based on the original
        FASTQ read counts provided by GTAC@MGI.

    df_merged : DataFrame
        A merged dataframe as provided by the MergeFastq class.

    merged_tsv : str
        A qualified path to a merged dataframe file as written by the
        MergeFastq class.

    target_counts : tuple
        A tuple of values use to evaluate the expected read count
        throughput of a given sample. Update this tuple to reduce or
        expand the number of evaluations.

    target_min_perct = int [80]
        The minimum percent of target sequence throughput for a sample
        to pass.

    Methods
    -------
    calc_gtac_read_coverage()
        Calculates the GTAC sample read count sequence throughput.

    write_df(file_path)
        Write the GTAC read count dataframe to a tab-delimited file.

    Examples
    --------
    read_counts = ReadCountsGtac(
        args=args,
        merged_tsv=merged_tsv
    )
    read_counts.calc_gtac_read_coverage()
    gtac_read_counts = Path(args.outdir) / 'gtac_read_counts.tsv'
    read_counts.write_df(file_path=str(gtac_read_counts.resolve()))
    """

    def __init__(self: Self, args: argparse.Namespace,
                 merged_tsv: str) -> None:
        """Construct the class.

        Parameters
        ----------
        args : argparse.Namespace
            Arguments for merging FASTQ files as provided by argparse.

        merged_tsv : str
            A qualified path to a merged dataframe file as written by
            the MergeFastq class.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        self.args = args
        self.merged_tsv = merged_tsv
        self.df_gtac_seqcov: DataFrame = DataFrame()
        self.df_merged: DataFrame = DataFrame()
        self.target_min_perct = 80
        self.target_counts: tuple = tuple()
        self.__set_target_coverages()
        self.__populate_df()
        return

    def __populate_df(self: Self) -> None:
        """Populate the merged samplemap dataframe.

        Parameters
        ----------
        None

        Raises
        ------
        FileNotFoundError
            Merged samplemap file does not exist.

        Returns
        -------
        None
        """
        merged_tsv_path = Path(self.merged_tsv)
        if merged_tsv_path.is_file() is False:
            raise FileNotFoundError(
                'Merged samplemap file does not exist.',
                self.merged_tsv
            )
        else:
            self.df_merged = pd.read_csv(self.merged_tsv, sep='\t')
        return

    def __set_target_coverages(self: Self) -> None:
        """Set the target coverage index lookup.

        We are populating a predefined set of target read counts values
        that will be used to evaluate actual read counts from FASTQ
        files.

        Parameters
        ----------
        None

        Raises
        ------
        None

        Returns
        -------
        None
        """
        target_read_counts = (
            100_000,
            200_000,
            300_000,
            400_000,
            500_000,
            1_000_000,
            1_500_000,
            2_000_000,
            2_500_000,
            3_000_000,
            3_500_000,
            4_000_000,
            4_500_000,
            5_000_000,
            10_000_000,
            20_000_000,
            30_000_000,
            40_000_000,
            50_000_000
        )
        self.target_counts = target_read_counts
        return

    def calc_gtac_read_coverage(self: Self) -> None:
        """Calculates the GTAC sample read count sequence throughput.

        Calculations are based on the original FASTQ read counts
        provided by GTAC@MGI in accompanying Samplemap.csv files. We
        calculate the percent of sequence throughput against a set of
        predefined target read throughput values--see
        __set_target_coverages() for details. Essentially:

            percent = (sample_counts / target_counts) * 100

        A minimum percent of coverage defines a pass/fail state at the
        sample read count level, set by the self.target_min_perct value.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Read counts differ for R1 and R2 columns.

        ValueError
            Read count columns vary in length.

        ValueError
            R1 and R2 read count sum differ froms sample count.

        Returns
        -------
        None
        """
        df = self.df_merged.copy()
        dfg = df.groupby('revised_sample_name')
        col_sample_name: list = list()
        col_r1_counts: list = list()
        col_r2_counts: list = list()
        col_sample_counts: list = list()
        col_target_min_perct: list = list()
        target_cols: dict = dict()

        for i, sample_name in enumerate(dfg.groups):
            col_sample_name.append(sample_name)
            col_target_min_perct.append(self.target_min_perct)
            dfi = dfg.get_group(sample_name)
            sample_counts = list(set(dfi['gtac_sample_reads']))[0]
            col_sample_counts.append(sample_counts)
            dfg_rn = dfi.groupby('read_number')
            r1_counts = list(set(
                dfg_rn.get_group(1)['gtac_end_pair_reads']
            ))[0]
            col_r1_counts.append(r1_counts)
            r2_counts = list(set(
                dfg_rn.get_group(2)['gtac_end_pair_reads']
            ))[0]
            col_r2_counts.append(r2_counts)
            target_cols[i] = {
                'col_target_counts': list(),
                'col_perct_of_target': list(),
                'col_is_pass_perct_target': list()
            }
            for target_counts in self.target_counts:
                target_cols[i]['col_target_counts'].append(target_counts)
                perct_of_target = round(
                    (sample_counts / target_counts) * 100, 2
                )
                target_cols[i]['col_perct_of_target'].append(perct_of_target)
                if perct_of_target < self.target_min_perct:
                    is_passed_perct_target = False
                else:
                    is_passed_perct_target = True
                target_cols[i]['col_is_pass_perct_target'].append(
                    is_passed_perct_target
                )

        if (col_r1_counts == col_r2_counts) is False:
            raise ValueError('Read counts differ for R1 and R2 columns.')

        if (
            len(col_sample_name) ==
            len(col_r1_counts) ==
            len(col_r2_counts) ==
            len(col_sample_counts) ==
            len(col_target_min_perct)
        ) is False:
            raise ValueError('Read count columns vary in length.')

        df_tmp = pd.DataFrame({
            'r1_counts': col_r1_counts,
            'r2_counts': col_r2_counts,
            'sample_counts': col_sample_counts
        })
        df_tmp['sum'] = df_tmp['r1_counts'] + df_tmp['r2_counts']
        if list(df_tmp['sample_counts']) == list(df_tmp['sum']) is False:
            raise ValueError(
                'R1 and R2 read count sum differ froms sample count.'
            )

        df_seqcov = pd.DataFrame()
        df_seqcov['sample_name'] = col_sample_name
        df_seqcov['R1_read_counts'] = col_r1_counts
        df_seqcov['R2_read_counts'] = col_r2_counts
        df_seqcov['sample_read_counts'] = col_sample_counts
        df_seqcov['min_target_perct_cov'] = col_target_min_perct

        collection: dict = dict()
        for i, target_count in enumerate(target_cols[0]['col_target_counts']):
            # We may use the first entry to setup the target labels for
            # all of the entries in the collection.
            collection[i] = {
                'perct_label': f'perct_of_{target_count}',
                'col_perct_of_target': list(),
                'is_passed_label': f'is_pased_{target_count}',
                'is_pass_perct_target': list()
            }

        for i in sorted(collection.keys()):
            for sample_i in sorted(target_cols):
                perct = target_cols[sample_i]['col_perct_of_target'][i]
                is_pass = target_cols[sample_i]['col_is_pass_perct_target'][i]
                collection[i]['col_perct_of_target'].append(perct)
                collection[i]['is_pass_perct_target'].append(is_pass)

        for i in sorted(collection):
            perct_label = collection[i]['perct_label']
            is_passed_label = collection[i]['is_passed_label']
            df_seqcov[perct_label] = collection[i]['col_perct_of_target']
            df_seqcov[is_passed_label] = collection[i]['is_pass_perct_target']

        self.df_gtac_seqcov = df_seqcov
        return

    def __calc_file_md5(self: Self, file_path: str) -> str:
        """Returns the MD5 hash for a specified file.

        Parameters
        ----------
        file_path : str
            A qualified file path for which to calculate the MD5 hash.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        with open(file_path, 'rb') as fh:
            data = fh.read()
            md5sum = hashlib.md5(data).hexdigest()
        return md5sum

    def write_df(self: Self, file_path: str) -> None:
        """Write the GTAC read count dataframe to a tab-delimited file.

        Parameters
        ----------
        file_path : str
            A qualified file path to write the read counts dataframe.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        self.df_gtac_seqcov.to_csv(file_path, sep='\t', index=False)
        md5 = self.__calc_file_md5(file_path=file_path)
        md5_path = file_path + '.MD5'
        with open(md5_path, 'w') as fh:
            fh.write(f'MD5 ({file_path}) = {md5}\n')
        return

# __END__
