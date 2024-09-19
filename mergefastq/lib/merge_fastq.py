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


class MergeFastq:

    def __init__(self: Self, args: argparse.Namespace, rename: RenameSamples,
                 samplemap: Samplemap) -> None:
        """Construct the class."""
        self.rename = rename
        self.samplemap = samplemap
        self.__parse_fastq_copy_types()
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

# __END__
