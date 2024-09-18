# Project     : merge_fastq
# File Name   : samplemap.py
# Description : Functions to handle samplemap files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Tue Sep 17 15:40:25 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from typing_extensions import Self


class Samplemap:

    def __init__(self: Self, args: argparse.Namespace) -> None:
        """Construct the class."""
        self.args = args
        self.__parse_samplemaps()
        return

    def __type_samplemap_format(self: Self, smap: str) -> str:
        """Return the Samplemap file's format type.

        GTAC@MGI provides a Samplemap.csv file with FASTQ files. This
        file provides basic information related to the sequencing batch.
        Unfortunately, this file has been inconsistent across projects
        and sequencing batches, and the fields provided are not always
        in the same order. As we are collecting a core subset of values
        needed for FASTQ merging, we will type the Samplemap file and
        deal with the columns based on the Samplemap format type.

        Parameters
        ----------
        smap : str
            Fully qualified path to a Samplemap file.

        Returns
        -------
        smap_type : str
            A defined format type for a given Samplemap file.
        """
        smap_type: str = str()
        df = pd.read_csv(smap)
        cols = df.columns
        return smap_type

    def __parse_samplemaps(self: Self) -> None:
        """Parse all of the input Samplemap files."""
        for smap in self.args.samplemap:
            smap_type = self.__type_samplemap_format(smap=smap)
            print(smap_type)
        return

# __END__
