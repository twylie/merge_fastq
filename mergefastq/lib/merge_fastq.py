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
        self.merge_cmds: dict = dict()
        self.__parse_fastq_copy_types()
        self.__setup_copy_cmds()
        self.__setup_merge_cmds()
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
        """Setup FASTQ commands not split across flow cells or lanes.

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
        # TODO: Is revised sample name the appropriate key?
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
                r1_copy_cmds.append(
                    f'md5sum {dest_r1_fq} > {dest_r1_fq}.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {dest_r1_fq} / 4 | bc '
                       f'> {dest_r1_fq}.counts')
                r1_copy_cmds.append(cmd)
            elif is_src_r1_fq_gzip is False:
                dest_r1_fq_stem = dest_r1_fq.as_posix()[:-3]
                r1_copy_cmds.append(f'cp {src_r1_fq} {dest_r1_fq_stem}')
                r1_copy_cmds.append(f'gzip {dest_r1_fq_stem}')
                r1_copy_cmds.append(
                    f'md5sum {dest_r1_fq_stem}.gz > {dest_r1_fq_stem}.gz.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {dest_r1_fq_stem}.gz / 4 | bc '
                       f'> {dest_r1_fq_stem}.gz.counts')
                r1_copy_cmds.append(cmd)

            r2_copy_cmds: list = list()
            if is_src_r2_fq_gzip is True:
                r2_copy_cmds.append(f'cp {src_r2_fq} {dest_r2_fq}')
                r2_copy_cmds.append(
                    f'md5sum {dest_r2_fq} > {dest_r2_fq}.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {dest_r2_fq} / 4 | bc '
                       f'> {dest_r2_fq}.counts')
                r2_copy_cmds.append(cmd)
            elif is_src_r2_fq_gzip is False:
                dest_r2_fq_stem = dest_r2_fq.as_posix()[:-3]
                r2_copy_cmds.append(f'cp {src_r2_fq} {dest_r2_fq_stem}')
                r2_copy_cmds.append(f'gzip {dest_r2_fq_stem}')
                r2_copy_cmds.append(
                    f'md5sum {dest_r2_fq_stem}.gz > {dest_r2_fq_stem}.gz.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {dest_r2_fq_stem}.gz / 4 | bc '
                       f'> {dest_r2_fq_stem}.gz.counts')
                r2_copy_cmds.append(cmd)

            self.copy_cmds.update({sample_name: (r1_copy_cmds, r2_copy_cmds)})
        return

    def __setup_merge_cmds(self: Self) -> None:
        """Setup FASTQ commands merging across flow cells and lanes.

        Any sample handled intheiin this routine will be split across
        multiple FASTQ files, and will require merging into a single,
        paired-end FASTQ file. The order of the first-in FASTQ in the
        merge is arbitrary; however, this needs to be consistent between
        the R1 and R2 ordering. That is, R1.a + R1.b should be the same
        order as R2.a + R2.b entries. Sorting the FASTQ files (like-end)
        by flowcell id should be able to keep the ordering the same
        between R1 and R2 files.

        Raises
        ------
        ValueError : Merge-samples should all have the same index
                     sequence.

        ValueError : Merge-samples sort order was not maintained.

        ValueError : Merge-samples revised sample names are not equal.

        ValueError : Merge-samples revised names must be unique.
        """
        df_smaps = self.samplemap.copy_df()
        # TODO: Is revised sample name the appropriate key?
        dfg_ids = df_smaps.groupby('revised_sample_name')
        for sample_name in self.merge_copy_ids:
            df = dfg_ids.get_group(sample_name)
            dfg_read_num = df.groupby('read_number')

            if len(set(df['index_sequence'])) != 1:
                raise ValueError(
                    'Merge-samples should all have the same index sequence.',
                    sample_name,
                    set(df['index_sequence'])
                )

            df_read_num_1 = dfg_read_num.get_group(1).sort_values(
                ['flow_cell_id', 'lane_number']
            )

            df_read_num_2 = dfg_read_num.get_group(2).sort_values(
                ['flow_cell_id', 'lane_number']
            )

            sort_order_1 = list(df_read_num_1[['flow_cell_id', 'lane_number']])
            sort_order_2 = list(df_read_num_2[['flow_cell_id', 'lane_number']])
            if sort_order_1 != sort_order_2:
                raise ValueError(
                    'Merge-samples sort order was not maintained.',
                    sort_order_1,
                    sort_order_2
                )

            src_r1_fq = list(df_read_num_1['fastq_path'].values)
            src_r2_fq = list(df_read_num_2['fastq_path'].values)
            name_r1 = list(df_read_num_1['revised_sample_name'].values)
            name_r2 = list(df_read_num_2['revised_sample_name'].values)

            if set(name_r1) != set(name_r2):
                raise ValueError(
                    'Merge-samples revised sample names are not equal.',
                    name_r1,
                    name_r2
                )

            name_r1_set = set(name_r1)
            name_r2_set = set(name_r2)

            if len(name_r1_set) > 1 or len(name_r2_set) > 1:
                raise ValueError(
                    'Merge-samples revised names must be unique.',
                    name_r1_set,
                    name_r2_set
                )

            name_r1_str = list(name_r1)[0]
            name_r2_str = list(name_r2)[0]
            dest_name_r1 = f'{name_r1_str}.R1.fastq.gz'
            dest_name_r2 = f'{name_r2_str}.R2.fastq.gz'
            sample_dir = Path(self.args.outdir) / Path(sample_name)
            dest_r1_fq = sample_dir / Path(dest_name_r1)
            dest_r2_fq = sample_dir / Path(dest_name_r2)

            r1_is_gzip: list = list()
            for r1_fq in src_r1_fq:
                r1_fq = '/tmp/sample1.100k.r2.fastq'  # FIXME:
                with gzip.open(r1_fq, 'r') as fhi:
                    try:
                        fhi.read(1)
                        is_src_r1_fq_gzip = True
                    except gzip.BadGzipFile:
                        is_src_r1_fq_gzip = False
                r1_is_gzip.append(is_src_r1_fq_gzip)

            r2_is_gzip: list = list()
            for r2_fq in src_r2_fq:
                r2_fq = '/tmp/sample1.100k.r2.fastq'  # FIXME:
                with gzip.open(r2_fq, 'r') as fhi:
                    try:
                        fhi.read(1)
                        is_src_r2_fq_gzip = True
                    except gzip.BadGzipFile:
                        is_src_r2_fq_gzip = False
                r2_is_gzip.append(is_src_r2_fq_gzip)

            # We may have a mixture of compressed and uncompressed
            # source FASTQ files, which complicates merging things. We
            # may have to copy and compress some of the source FASTQ
            # files prior to merging.

            r1_merge_cmds: list = list()
            tmp_r1: set = set()
            r1_cat_order: list = list()
            for i, is_gzip in enumerate(r1_is_gzip):
                if is_gzip is True:
                    r1_cat_order.append(src_r1_fq[i])
                elif is_gzip is False:
                    from_r1_fq = src_r1_fq[i][:-3]
                    to_r1_fq = sample_dir / Path(src_r1_fq[i]).stem
                    tmp_r1.add(to_r1_fq)
                    r1_merge_cmds.append(f'cp {from_r1_fq} {to_r1_fq}')
                    r1_merge_cmds.append(f'gzip {to_r1_fq}')
                    to_r1_gzip_fq = str(to_r1_fq) + '.gz'
                    r1_cat_order.append(to_r1_gzip_fq)
            cat_files = ' '.join(r1_cat_order)
            r1_merge_cmds.append(f'cat {cat_files} > {dest_r1_fq}')
            r1_merge_cmds.append(f'md5sum {dest_r1_fq} > {dest_r1_fq}.MD5')
            cmd = (f'echo $(zgrep -c ^ {dest_r1_fq} / 4 | bc '
                   f'> {dest_r1_fq}.counts')
            r1_merge_cmds.append(cmd)
            if tmp_r1:
                for tmp_file in tmp_r1:
                    r1_merge_cmds.append(f'rm {tmp_file}')

            r2_merge_cmds: list = list()
            tmp_r2: set = set()
            r2_cat_order: list = list()
            for i, is_gzip in enumerate(r2_is_gzip):
                if is_gzip is True:
                    r2_cat_order.append(src_r2_fq[i])
                elif is_gzip is False:
                    from_r2_fq = src_r2_fq[i][:-3]
                    to_r2_fq = sample_dir / Path(src_r2_fq[i]).stem
                    tmp_r2.add(to_r2_fq)
                    r2_merge_cmds.append(f'cp {from_r2_fq} {to_r2_fq}')
                    r2_merge_cmds.append(f'gzip {to_r2_fq}')
                    to_r2_gzip_fq = str(to_r2_fq) + '.gz'
                    r2_cat_order.append(to_r2_gzip_fq)
            cat_files = ' '.join(r2_cat_order)
            r2_merge_cmds.append(f'cat {cat_files} > {dest_r2_fq}')
            r2_merge_cmds.append(f'md5sum {dest_r2_fq} > {dest_r2_fq}.MD5')
            cmd = (f'echo $(zgrep -c ^ {dest_r2_fq} / 4 | bc '
                   f'> {dest_r2_fq}.counts')
            r2_merge_cmds.append(cmd)
            if tmp_r2:
                for tmp_file in tmp_r2:
                    r2_merge_cmds.append(f'rm {tmp_file}')

            self.merge_cmds.update({
                sample_name: (r1_merge_cmds, r2_merge_cmds)
            })
        return

# __END__
