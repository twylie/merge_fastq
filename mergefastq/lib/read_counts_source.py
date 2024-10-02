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
from typing_extensions import Self
from pathlib import Path
import hashlib


class ReadCountsSource():

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
        """Update the source sample counts in the dataframe."""
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
        """Update the source read counts in the dataframe."""
        self.__col_src_end_pair_reads()
        self.__col_src_sample_reads()
        return

    def __calc_src_read_coverage(self: Self) -> None:
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
        """Compare the read counts for source and gtac reports."""
        df_src = self.df_src_seqcov.copy()
        df_gtac = self.df_gtac_counts.copy()
        df_eval = df_src == df_gtac

        r1_false: list = list()
        false_i = df_eval.loc[df_eval['R1_read_counts']].index
        if len(false_i) >= 1:
            r1_false += list(false_i)

        r2_false: list = list()
        false_i = df_eval.loc[~df_eval['R2_read_counts']].index
        if len(false_i) >= 1:
            r2_false += list(false_i)

        sample_false: list = list()
        false_i = df_eval.loc[~df_eval['sample_read_counts']].index
        if len(false_i) >= 1:
            sample_false += list(false_i)

        src_sample_names = df_src.loc[r1_false]['sample_name']
        src_r1_counts = df_src.loc[r1_false]['R1_read_counts']
        gtac_sample_names = df_gtac.loc[r1_false]['sample_name']
        gtac_r1_counts = df_gtac.loc[r1_false]['R1_read_counts']
        for i in zip(
                src_sample_names,
                src_r1_counts,
                gtac_sample_names,
                gtac_r1_counts
        ):
            aa, bb, cc, dd = i
            print(f'R1_sample_counts: {aa} {bb} != {dd}')
        return

# __END__
