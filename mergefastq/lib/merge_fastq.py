# Project     : merge_fastq
# File Name   : merge_fastq.py
# Description : Functions for merging FASTQ files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Thu Sep 19 17:01:30 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from pandas import DataFrame  # type: ignore
from . rename_samples import RenameSamples  # type: ignore
from . samplemap import Samplemap  # type: ignore
from typing_extensions import Self
from pathlib import Path
import re
import hashlib
import gzip


class MergeFastq:

    def __init__(self: Self, args: argparse.Namespace, rename: RenameSamples,
                 samplemap: Samplemap) -> None:
        """Construct the class."""
        self.args = args
        self.rename = rename
        self.samplemap = samplemap
        self.copy_cmds: dict = dict()
        self.__parse_fastq_copy_types()
        self.__setup_copy_cmds()
        # self.__setup_merge_cmds()
        return

    def __parse_fastq_copy_types(self: Self) -> None:
        """Parse the Samplemap dataframe and determine FASTQ copy types.

        We may determine the copy types for samples based on the
        cardinality of sample ids to flowcells. For example, if a unique
        sample has 2 FASTQ (R1 & R2) and 1 associated flowcell, we don't
        need to merge the FASTQ, but rather simply copy them using the
        revised sample name. Samples with more than this will require
        merging of FASTQ files.

        Raises
        ------
        ValueError : FASTQ copy type counts differ.
        """
        df = self.samplemap.copy_df()
        dfg = df.groupby('revised_sample_name')
        single_copy_ids: set = set()
        merge_copy_ids: set = set()
        for i, sample_name in enumerate(dfg.groups, 1):
            df_sample_name = dfg.get_group(sample_name)
            flow_cell_count = len(df_sample_name['flow_cell_id'].unique())
            fastq_count = len(df_sample_name)
            if flow_cell_count == 1 and fastq_count == 2:
                single_copy_ids.add(sample_name)
            else:
                merge_copy_ids.add(sample_name)
        if (len(single_copy_ids) + len(merge_copy_ids)) != i:
            raise ValueError('FASTQ copy type counts differ.')
        else:
            self.single_copy_ids = single_copy_ids
            self.merge_copy_ids = merge_copy_ids
        return

    def __setup_copy_cmds(self: Self) -> None:
        """Setup commands for FASTQ not split across flow cells or lanes.

        We will handle sample ids that are of single paired-end copying
        type---i.e. they do not require merging of FASTQ file. We will
        be making a copy of these files for downstream analysis, using a
        new file name for the FASTQ files. Commands are being collected
        and associated on a per-sample level.

        We will handle source FASTQ files that are potentially (1)
        compressed or (2) decompressed or (3) a mixture of compressed
        and decompressed. Compression is handled using gzip.

        In addition to copying the source to destination FASTQ, we will
        also count the number of reads in the FASTQ file and write a
        corresponding "counts" file.

        Raises
        ------
        ValueError : Copy type should only have 2 associated FASTQ files.

        ValueError : R1 and R2 index sequences do not match.

        ValueError : R1 and R2 flow cell ids do not match.

        ValueError : R1 and R2 lane numbers do not match.
        """
        df_smaps = self.samplemap.copy_df()
        dfg_ids = df_smaps.groupby('revised_sample_name')
        for sample_name in self.single_copy_ids:
            df = dfg_ids.get_group(sample_name)

            fastq_count = len(df.index)
            if fastq_count != 2:
                raise ValueError(
                    'Copy type should only have 2 associated FASTQ files.',
                    sample_name,
                    fastq_count
                )

            dfg_read_num = df.groupby('read_number')
            r1_index = dfg_read_num.get_group(1)['index_sequence'].item()
            r2_index = dfg_read_num.get_group(2)['index_sequence'].item()
            if r1_index != r2_index:
                raise ValueError(
                    'R1 and R2 index sequences do not match.',
                    sample_name,
                    r1_index,
                    r2_index
                )

            r1_fc_id = dfg_read_num.get_group(1)['flow_cell_id'].item()
            r2_fc_id = dfg_read_num.get_group(2)['flow_cell_id'].item()
            if r1_fc_id != r2_fc_id:
                raise ValueError(
                    'R1 and R2 flow cell ids do not match.',
                    sample_name,
                    r1_fc_id,
                    r2_fc_id
                )

            r1_lane_num = dfg_read_num.get_group(1)['lane_number'].item()
            r2_lane_num = dfg_read_num.get_group(2)['lane_number'].item()
            if r1_lane_num != r2_lane_num:
                raise ValueError(
                    'R1 and R2 lane numbers do not match.',
                    sample_name,
                    r1_lane_num,
                    r2_lane_num
                )

            src_r1_fq = dfg_read_num.get_group(1)['fastq_path'].item()
            src_r2_fq = dfg_read_num.get_group(2)['fastq_path'].item()
            name_r1 = dfg_read_num.get_group(1)['revised_sample_name'].item()
            name_r2 = dfg_read_num.get_group(2)['revised_sample_name'].item()
            dest_name_r1 = f'{name_r1}.R1.fastq.gz'
            dest_name_r2 = f'{name_r2}.R2.fastq.gz'
            sample_dir = Path(self.args.outdir) / Path(sample_name)
            dest_r1_fq = sample_dir / Path(dest_name_r1)
            dest_r2_fq = sample_dir / Path(dest_name_r2)

            # FIXME:
            # src_r1_fq = '/Users/toddwylie/Desktop/DEVELOPMENT/virosearchTutorial/virosearch/testing/sample1.100k.r1.fastq.gz'
            # src_r2_fq = '/Users/toddwylie/Desktop/DEVELOPMENT/virosearchTutorial/virosearch/testing/sample1.100k.r2.fastq.gz'

            src_r1_fq = '/tmp/sample1.100k.r2.fastq'
            src_r2_fq = '/tmp/sample1.100k.r2.fastq'

            with gzip.open(src_r1_fq, 'r') as fhi:
                try:
                    fhi.read(1)
                    is_src_r1_fq_gzip = True
                except gzip.BadGzipFile:
                    is_src_r1_fq_gzip = False

            with gzip.open(src_r2_fq, 'r') as fhi:
                try:
                    fhi.read(1)
                    is_src_r2_fq_gzip = True
                except gzip.BadGzipFile:
                    is_src_r2_fq_gzip = False

            # The read count commands look a little busy, but this
            # method should be a quick way to count the lines and
            # calculate the read count per FASTQ file.

            r1_copy_cmds: list = list()
            if is_src_r1_fq_gzip is True:
                r1_copy_cmds.append(f'cp {src_r1_fq} {dest_r1_fq}')
            else:
                dest_r1_fq = dest_r1_fq.as_posix()[:-3]
                r1_copy_cmds.append(f'cp {src_r1_fq} {dest_r1_fq}')
                r1_copy_cmds.append(f'gzip {dest_r1_fq}')
            r1_copy_cmds.append(
                f'md5sum {dest_r1_fq}.gz > {dest_r1_fq}.gz.MD5'
            )
            cmd = (f'echo $(zgrep -c ^ {dest_r1_fq}.gz / 4 | bc '
                   f'> {dest_r1_fq}.gz.counts')
            r1_copy_cmds.append(cmd)

            r2_copy_cmds: list = list()
            if is_src_r2_fq_gzip is True:
                r2_copy_cmds.append(f'cp {src_r2_fq} {dest_r2_fq}')
            else:
                dest_r2_fq = dest_r2_fq.as_posix()[:-3]
                r2_copy_cmds.append(f'cp {src_r2_fq} {dest_r2_fq}')
                r2_copy_cmds.append(f'gzip {dest_r2_fq}')
            r2_copy_cmds.append(
                f'md5sum {dest_r2_fq}.gz > {dest_r2_fq}.gz.MD5'
            )
            cmd = (f'echo $(zgrep -c ^ {dest_r2_fq}.gz / 4 | bc '
                   f'> {dest_r2_fq}.gz.counts')
            r2_copy_cmds.append(cmd)

            self.copy_cmds.update({
                sample_name: zip(r1_copy_cmds, r2_copy_cmds)
            })
        return

# __END__
