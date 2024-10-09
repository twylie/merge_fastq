# Project     : merge_fastq
# File Name   : read_counts_source.py
# Description : Methods for sequence throughput evaluation.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Tue Oct 01 16:14:01 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from pandas import DataFrame
import numpy as np
from typing_extensions import Self
from pathlib import Path
import hashlib


class ReadCountsSource():
    """A class for FASTQ sequence throughput evaluation.

    This class provides methods for calculating and evaluating source
    FASTQ file read counts. We will handle read counts as provided by
    sample_name.fastq.gz.counts files, see MergeFastq class for details.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    Attributes
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    df_comp : DataFrame
        A dataframe comparing GTAC FASTQ read counts to source FASTQ
        read counts.

    df_gtac_counts : DataFrame
        A dataframe of GTAC FASTQ read counts.

    df_merged : DataFrame
        A dataframe of merged FASTQ file information.

    df_src_seqcov : DataFrame
        A dataframe of source FASTQ read counts and sequence coverage
        metrics.

    gtac_counts : str
        A qualified path to the GTAC read counts file.

    merged_samplemap : str
        A qualified path to the merged FASTQ samplemap file.

    target_counts : tuple
        A tuple of values use to evaluate the expected read count
        throughput of a given sample. Update this tuple to reduce or
        expand the number of evaluations.

    target_min_perct = int [80]
        The minimum percent of target sequence throughput for a sample
        to pass.

    report_path : str
        Path to the comparison read counts report file.

    Methods
    -------
    tell_comp_differences()
        If there are differences between GTAC and source FASTQ read
        counts, tell the user which differences are present per sample.

    Examples
    --------
    read_counts = ReadCountsSource(args=args)
    read_counts.tell_comp_differences()
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
        self.merged_samplemap = self.args.merged_samplemap
        self.gtac_counts = self.args.gtac_counts
        self.df_merged: DataFrame = DataFrame()
        self.df_src_seqcov: DataFrame = DataFrame()
        self.df_comp: DataFrame = DataFrame()
        self.df_gtac_counts: DataFrame = DataFrame()
        self.target_min_perct = 80
        self.target_counts: tuple = tuple()
        self.__set_target_coverages()
        self.__populate_gtac_df()
        self.__populate_merged_df()
        self.__update_merged_df_read_counts()
        self.__calc_src_read_coverage()
        src_tsv = Path(self.args.outdir) / 'src_read_counts.tsv'
        self.__write_src_counts_df(file_path=str(src_tsv.resolve()))
        self.__compare_gtac_to_src_counts()
        report_path = Path(self.args.outdir) / 'read_count_comparisons.tsv'
        self.report_path = str(report_path.resolve())
        self.__write_comp_counts_df(file_path=self.report_path)
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

    def __populate_merged_df(self: Self) -> None:
        """Populate the merged FASTQ dataframe.

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
        merged_tsv_path = Path(self.merged_samplemap)
        if merged_tsv_path.is_file() is False:
            raise FileNotFoundError(
                'Merged samplemap file does not exist.',
                self.merged_samplemap
            )
        else:
            self.df_merged = pd.read_csv(self.merged_samplemap, sep='\t')
        return

    def __col_src_end_pair_reads(self: Self) -> None:
        """Add the source FASTQ read counts column in dataframe.

        Parameters
        ----------
        None

        Raises
        ------
        FileNotFoundError
            FASTQ counts file not found.

        ValueError
            Read counts column length and dataframe length are not
            equal.

        Returns
        -------
        None
        """
        df = self.df_merged.copy()
        col_counts: list = list()
        for i in df.index:
            dfi = df.loc[i]
            fastq_counts = dfi['merged_fastq_path'] + '.counts'
            if Path(fastq_counts).is_file() is False:
                raise FileNotFoundError(
                    'FASTQ counts file not found.',
                    fastq_counts
                )
            else:
                with open(fastq_counts, 'r') as fhi:
                    for line in fhi:
                        read_counts = line.strip()
                        col_counts.append(int(read_counts))

        if len(col_counts) == len(df.index) is False:
            raise ValueError('Read counts column length and dataframe '
                             'length are not equal.')
        else:
            self.df_merged['src_end_pair_reads'] = col_counts
        return

    def __col_src_sample_reads(self: Self) -> None:
        """Update the source sample counts in the dataframe.

        We will update the merged, sample-level read counts based on the
        independent, manual FASTQ read count review. The sample read
        counts includes both the R1 and R2 FASTQ file read counts.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Odd length of reads pairs.

        ValueError
            Sample counts column length does not match dataframe.

        ValueError
            Source sample counts do not match for R1 & R2.

        ValueError
            Source end pair counts do not match for R1 & R2.

        Returns
        -------
        None
        """
        df = self.df_merged.copy()

        count_index: dict = dict()
        df = self.df_merged.copy()
        dfg = df.groupby(['sample_name', 'read_number'])
        for i in dfg.groups:
            sample_name, read_number = i
            dfi = dfg.get_group(i)
            if len(dfi['src_end_pair_reads']) == 1:
                sample_counts = dfi['src_end_pair_reads'].item() * 2
            elif len(dfi['src_end_pair_reads']) == 2:
                sample_counts = dfi.sum()['src_end_pair_reads']
            else:
                raise ValueError(
                    'Odd length of reads pairs.',
                    len(dfi['src_end_pair_reads'])
                )
            count_index[sample_name] = {'sample_counts': int(sample_counts)}

        col_sample_counts: list = list()
        for i in df.index:
            sample_name = df.loc[i]['sample_name']
            sample_counts = count_index[sample_name]['sample_counts']
            col_sample_counts.append(sample_counts)

        if len(col_sample_counts) != len(self.df_merged.index):
            raise ValueError(
                'Sample counts column length does not match dataframe.'
            )

        self.df_merged['src_sample_reads'] = col_sample_counts

        dfg = self.df_merged.groupby('revised_sample_name')
        for sample_name in dfg.groups:
            dfi = dfg.get_group(sample_name)
            dfg_r1_fastq = dfi.groupby('read_number').get_group(1)
            r1_sample_counts = list(set(dfg_r1_fastq['src_sample_reads']))[0]
            r1_ep_counts = list(set(dfg_r1_fastq['src_end_pair_reads']))[0]
            dfg_r2_fastq = dfi.groupby('read_number').get_group(2)
            r2_sample_counts = list(set(dfg_r2_fastq['src_sample_reads']))[0]
            r2_ep_counts = list(set(dfg_r2_fastq['src_end_pair_reads']))[0]
            if r1_sample_counts != r2_sample_counts:
                raise ValueError(
                    'Source sample counts do not match for R1 & R2.',
                    sample_name,
                    f'R1={r1_sample_counts}',
                    f'R2={r2_sample_counts}'
                )
            if r1_ep_counts != r2_ep_counts:
                raise ValueError(
                    'Source end pair counts do not match for R1 & R2.',
                    sample_name,
                    f'R1={r1_ep_counts}',
                    f'R2={r2_ep_counts}'
                )
        return

    def __update_merged_df_read_counts(self: Self) -> None:
        """Update the source read counts in the dataframe.

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
        self.__col_src_end_pair_reads()
        self.__col_src_sample_reads()
        return

    def __calc_src_read_coverage(self: Self) -> None:
        """Calculate the source sequencing read coverage values.

        Calculations are based on the independent, manual FASTQ read
        counts provided by the sample_name.fastq.gz.counts files
        generated during MergeFastq processing. We calculate the percent
        of sequence throughput against a set of predefined target read
        throughput values--see __set_target_coverages() for details.
        Essentially:

            percent = ((sample_counts / 2) / target_counts) * 100

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

        Return
        ------
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
            sample_counts = list(set(dfi['src_sample_reads']))[0]
            col_sample_counts.append(sample_counts)
            dfg_rn = dfi.groupby('read_number')
            r1_counts = list(set(
                dfg_rn.get_group(1)['src_end_pair_reads']
            ))[0]
            col_r1_counts.append(r1_counts)
            r2_counts = list(set(
                dfg_rn.get_group(2)['src_end_pair_reads']
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
                    ((sample_counts / 2) / target_counts) * 100, 2
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
                'is_passed_label': f'is_passed_{target_count}',
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

        self.df_src_seqcov = df_seqcov
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

    def __write_src_counts_df(self: Self, file_path: str) -> None:
        """Write the source read count dataframe to a tab-delimited file.

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
        self.df_src_seqcov.to_csv(file_path, sep='\t', index=False)
        md5 = self.__calc_file_md5(file_path=file_path)
        md5_path = file_path + '.MD5'
        with open(md5_path, 'w') as fh:
            fh.write(f'MD5 ({file_path}) = {md5}\n')
        return

    def __write_comp_counts_df(self: Self, file_path: str) -> None:
        """Write the count comparisons dataframe to a tab-delimited file.

        Parameters
        ----------
        file_path : str
            A qualified file path to write the read counts comparison
            dataframe.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        self.df_comp.to_csv(file_path, sep='\t', index=False)
        md5 = self.__calc_file_md5(file_path=file_path)
        md5_path = file_path + '.MD5'
        with open(md5_path, 'w') as fh:
            fh.write(f'MD5 ({file_path}) = {md5}\n')
        return

    def __populate_gtac_df(self: Self) -> None:
        """Populate the GTAC FASTQ dataframe.

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
        gtac_tsv_path = Path(self.gtac_counts)
        if gtac_tsv_path.is_file() is False:
            raise FileNotFoundError(
                'GTAC counts file does not exist.',
                self.gtac_counts
            )
        else:
            self.df_gtac_counts = pd.read_csv(self.gtac_counts, sep='\t')
        return

    def __compare_gtac_to_src_counts(self: Self) -> None:
        """Compare the read counts for source and GTAC reports.

        Given that we have the two read count throughput tables
        (gtac_read_counts.tsv and src_read_counts.tsv), we may compare
        the GTAC-based FASTQ read counts to those calculated by
        independent, manual review. If there are any differences in read
        counts (i.e. R1 counts, R2 counts, and sample-level counts)
        between the two tables will be highlighted in a subsequent
        comparison report file.

         The comparison report file's fields will be blank if there are
        differences between the read count tables. Differences will be
        displayed as GTAC:SRC count pairs. The is_no_difference column
        is a boolean field indicating if there are any differences
        across the count fields for the sample.

        REPORT FIELDS
        -------------
        1. sample_name        : STR
        2  R1_read_counts     : STR | NaN
        3. R2_read_counts     : STR | NaN
        4. sample_read_counts : STR | NaN
        5. is_no_difference   : BOOL

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Read count columns vary in length.

        Returns
        -------
        None
        """
        df_src = self.df_src_seqcov.copy()
        df_gtac = self.df_gtac_counts.copy()
        df_eval = df_src == df_gtac
        col_sample_names: list = list()
        col_r1_read_counts: list = list()
        col_r2_read_counts: list = list()
        col_sample_read_counts: list = list()

        for i in df_eval.index:
            sample_name = df_gtac.loc[i]['sample_name']
            col_sample_names.append(sample_name)
            is_r1_same = bool(df_eval.loc[i]['R1_read_counts'])
            is_r2_same = bool(df_eval.loc[i]['R2_read_counts'])
            is_sample_same = bool(df_eval.loc[i]['sample_read_counts'])
            if is_r1_same is False:
                gtac_counts = df_gtac.loc[i]['R1_read_counts']
                src_counts = df_src.loc[i]['R1_read_counts']
                r1_comp = f'{gtac_counts}:{src_counts}'
                col_r1_read_counts.append(r1_comp)
            elif is_r1_same is True:
                col_r1_read_counts.append(np.nan)
            if is_r2_same is False:
                gtac_counts = df_gtac.loc[i]['R2_read_counts']
                src_counts = df_src.loc[i]['R2_read_counts']
                r2_comp = f'{gtac_counts}:{src_counts}'
                col_r2_read_counts.append(r2_comp)
            elif is_r2_same is True:
                col_r2_read_counts.append(np.nan)
            if is_sample_same is False:
                gtac_counts = df_gtac.loc[i]['sample_read_counts']
                src_counts = df_src.loc[i]['sample_read_counts']
                r2_comp = f'{gtac_counts}:{src_counts}'
                col_sample_read_counts.append(r2_comp)
            elif is_r2_same is True:
                col_sample_read_counts.append(np.nan)

        self.df_comp['sample_name'] = col_sample_names
        self.df_comp['R1_read_counts'] = col_r1_read_counts
        self.df_comp['R2_read_counts'] = col_r2_read_counts
        self.df_comp['sample_read_counts'] = col_sample_read_counts

        col_bool: list = list()
        for i in df_eval.index:
            dfi = df_eval.loc[i][[
                'R1_read_counts',
                'R2_read_counts',
                'sample_read_counts'
            ]]
            is_all_true = bool(dfi.all())
            col_bool.append(is_all_true)

        if (
                len(col_sample_names) ==
                len(col_r1_read_counts) ==
                len(col_r2_read_counts) ==
                len(col_sample_read_counts) ==
                len(col_bool)
        ) is False:
            raise ValueError('Read count columns vary in length.')
        else:
            self.df_comp['is_no_difference'] = col_bool
        return

    def tell_comp_differences(self: Self) -> None:
        """Tell the user the comparison differences.

        If there are differences between GTAC and source FASTQ read
        counts, tell the user which differences are present per sample
        using STDOUT.

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
        df = self.df_comp.copy()
        if bool(df['is_no_difference'].all()) is True:
            print('\nFINISHED:')
            print('There are no differences between the GTAC and '
                  'source read counts.')
            print('\nOUTPUT:')
            print(f'{self.report_path}\n')
        else:
            dfg = df.groupby('is_no_difference')
            dfg_false = dfg.get_group(False)
            print('\nFINISHED:')
            print('Differences detected between the GTAC and source read '
                  'counts.\n')
            print('Investigate the following sample names:\n')
            for i in dfg_false.index:
                sample_name = dfg_false.loc[i]['sample_name']
                print(sample_name)
            print('\nOUTPUT:')
            print(f'{self.report_path}\n')
        return

# __END__
