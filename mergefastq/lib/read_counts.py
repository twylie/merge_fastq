# Project     : merge_fastq
# File Name   : read_counts.py
# Description : Methods for sequence throughput evaluation.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Fri Sep 27 12:21:33 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from typing_extensions import Self
from pathlib import Path


class ReadCounts:
    """A class for read count evaluation."""

    def __init__(self: Self, args: argparse.Namespace,
                 merged_tsv: str) -> None:
        """Construct the class."""
        self.args = args
        self.merged_tsv = merged_tsv
        self.__set_target_coverages()
        self.__populate_df()
        self.target_min_perct = 0.05  # FIXME
        return

    def __populate_df(self: Self) -> None:
        """Populate the merged samplemap dataframe.

        Raises
        ------
        FileNotFoundError : Merged samplemap file does not exist.
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
        """
        target_read_counts = (
            200_000,
            500_000,
            1_000_000,
            2_000_000,
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
        """Calculates the GTAC sample read count sequence throughput."""
        # FIXME: Needs evaluations throughout.
        df = self.df_merged.copy()
        dfg = df.groupby('revised_sample_name')
        col_sample_name: list = list()
        col_r1_counts: list = list()
        col_r2_counts: list = list()
        col_sample_counts: list = list()
        col_target_min_perct: list = list()
        target_cols: dict = dict()
        for i, sample_name in enumerate(dfg.groups):
            # FIXME: We need labels for the target counts values.
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
                    (sample_counts / target_counts) * 100, 1
                )
                target_cols[i]['col_perct_of_target'].append(perct_of_target)
                if perct_of_target < self.target_min_perct:
                    is_passed_perct_target = False
                else:
                    is_passed_perct_target = True
                target_cols[i]['col_is_pass_perct_target'].append(
                    is_passed_perct_target
                )
        df_seqcov = pd.DataFrame()
        df_seqcov['sample_name'] = col_sample_name
        df_seqcov['R1_read_counts'] = col_r1_counts
        df_seqcov['R2_read_counts'] = col_r2_counts
        df_seqcov['sample_read_counts'] = col_sample_counts
        df_seqcov['min_target_perct_cov'] = col_target_min_perct
        collection: dict = dict()
        for i, target_count in enumerate(target_cols[0]['col_target_counts']):
            collection[i] = {
                'perct_label': f'perct_of_{target_count}',
                'col_perct_of_target': list(),
                'is_passed_label': f'is_pased_{target_count}',
                'is_pass_perct_target': list()
            }
        for i in sorted(collection.keys()):
            for sample in sorted(target_cols):
                perct = target_cols[sample]['col_perct_of_target'][i]
                is_pass = target_cols[sample]['col_is_pass_perct_target'][i]
                collection[i]['col_perct_of_target'].append(perct)
                collection[i]['is_pass_perct_target'].append(is_pass)
        for i in sorted(collection):
            perct_label = collection[i]['perct_label']
            is_passed_label = collection[i]['is_passed_label']
            df_seqcov[perct_label] = collection[i]['col_perct_of_target']
            df_seqcov[is_passed_label] = collection[i]['is_pass_perct_target']
        self.df_gtac_seqcov = df_seqcov
        return

# __END__
