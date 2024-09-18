# Project     : merge_fastq
# File Name   : samplemap.py
# Description : Functions to handle samplemap files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Tue Sep 17 15:40:25 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import argparse
import pandas as pd  # type: ignore
from pandas import DataFrame  # type: ignore
from typing_extensions import Self
from pathlib import Path
import re


class Samplemap:

    def __init__(self: Self, args: argparse.Namespace) -> None:
        """Construct the class."""
        self.args = args
        self.smaps: dict = dict()
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

        NOTE: There are 2 samplemap files provided by GTAC@MGI
        (Samplemap.csv and Samplemap2.csv). We will always expect the
        Samplemap.csv file as input.

        Parameters
        ----------
        smap : str
            Fully qualified path to a Samplemap file.

        Format Types
        ------------
        smap_mid_2024_format:
             1. FASTQ
             2. Flowcell ID
             3. Index Sequence
             4. Flowcell Lane
             5. ESP ID
             6. Pool Name
             7. Species
             8. Illumina Sample Type
             9. Library Type
            10. Library Name
            11. Date Complete
            12. Total Reads
            13. Total Bases
            14. PhiX Error Rate
            15. % Pass Filter Clusters
            16. % >Q30
            17. Avg Q Score

        Returns
        -------
        smap_type : str
            A defined format type for a given Samplemap file.

        Raises
        ------
        ValueError : File name does not match Samplemap.csv format.

        TypeError : Unknown samplemap format type encountered.
        """
        smap_file_name = Path(smap).stem + Path(smap).suffix
        if smap_file_name != 'Samplemap.csv':
            raise ValueError(
                'File name does not match Samplemap.csv format.',
                smap_file_name
            )

        # The order of the columns is irrelevant as long as all fields
        # are present in the Samplemap.csv file.

        df = pd.read_csv(smap)
        cols = set(df.columns)

        smap_mid_2024_format = {'FASTQ', 'Flowcell ID', 'Index Sequence',
                                'Flowcell Lane', 'ESP ID', 'Pool Name',
                                'Species', 'Illumina Sample Type',
                                'Library Type', 'Library Name',
                                'Date Complete', 'Total Reads', 'Total Bases',
                                'PhiX Error Rate', '% Pass Filter Clusters',
                                '% >Q30', 'Avg Q Score'}

        if cols == smap_mid_2024_format:
            smap_type = 'smap_mid_2024_format'
        else:
            raise TypeError(
                'Unknown samplemap format type encountered.',
                cols
            )
        return smap_type

    def __smap_mid_2024_to_df(self: Self, smap: str) -> DataFrame:
        """Return dataframe of Samplemap information.

        We only require a uniform, subset of Samplemap values for FASTQ
        file merging and sequence throughput assessment. We will parse
        the desired fields from the Samplemap file and place them within
        a dataframe for downstream processing.

        We are making the following assumptions regarding the new
        2024 Samplemap formats:

        1. Library Name == Sample Name.
        2. Library Name is a unique identifier.
        3. The Library Name cardinality is one-to-one for sample names.

        Samplemap Format Type
        ---------------------
         1. FASTQ
         2. Flowcell ID
         3. Index Sequence
         4. Flowcell Lane
         5. ESP ID
         6. Pool Name
         7. Species
         8. Illumina Sample Type
         9. Library Type
        10. Library Name
        11. Date Complete
        12. Total Reads
        13. Total Bases
        14. PhiX Error Rate
        15. % Pass Filter Clusters
        16. % >Q30
        17. Avg Q Score

        Example Samplemap Input Entry
        -----------------------------
        {'% >Q30': 94.0,
         '% Pass Filter Clusters': 67.26,
         'Avg Q Score': 38.81,
         'Date Complete': '2024-03-18 01:45:24.403828+00:00',
         'ESP ID': 'LIB028751-DIL01',
         'FASTQ': 'LIB028751-DIL01_2235JKLT4_S210_L006_R1_001.fastq.gz',
         'Flowcell ID': '2235JKLT4',
         'Flowcell Lane': 6,
         'Illumina Sample Type': nan,
         'Index Sequence': 'CCTCTATAGA-GTACTTCCAA',
         'Library Name': 'H141_1_3_2',
         'Library Type': 'WGS',
         'PhiX Error Rate': 0.28,
         'Pool Name': nan,
         'Species': 'human',
         'Total Bases': '22,707,956,820',
         'Total Reads': '75,191,910'}

        Example Dataframe Output Entry
        ------------------------------
        {'esp_id': 'LIB028751-DIL01',
         'fastq': 'LIB028751-DIL01_2235JKLT4_S210_L006_R1_001.fastq.gz',
         'flow_cell_id': '2235JKLT4',
         'index_sequence': 'CCTCTATAGA-GTACTTCCAA',
         'lane_number': 6,
         'library_type': 'WGS',
         'pool_name': nan,
         'read_number': 1,
         'sample_name': 'H141_1_3_2',
         'samplemap_path': '/pwd/Samplemap.csv',
         'total_bases': 22707956820,
         'total_reads': 75191910}

        Raises
        ------
        ValueError : Read-pair number tag is unknown.

        ValueError : Dataframe field row lengths not of equal length.
        """
        (file_name_col,
         flow_cell_id_col,
         index_sequence_col,
         lane_number_col,
         read_number_col,
         sample_name_col,
         library_type_col,
         total_bases_col,
         samplemap_path_col,
         total_reads_col,
         esp_id_col,
         pool_name_col) = (list() for i in range(12))

        df = pd.read_csv(smap, thousands=',')
        for i in df.index:
            df_i = df.loc[i]
            fastq = df_i['FASTQ']
            flow_cell_id = df_i['Flowcell ID']
            index_sequence = df_i['Index Sequence']
            lane_number = df_i['Flowcell Lane']
            sample_name = df_i['Library Name']
            library_type = df_i['Library Type']
            total_bases = df_i['Total Bases']
            samplemap_path = smap
            total_reads = df_i['Total Reads']
            esp_id = df_i['ESP ID']
            pool_name = df_i['Pool Name']
            read_num_tag = re.search('_R[12]_', df_i['FASTQ']).group()
            match read_num_tag:
                case '_R1_':
                    read_number = 1
                case '_R2_':
                    read_number = 2
                case _:
                    raise ValueError(
                        'Read-pair number tag is unknown.',
                        read_num_tag
                    )

            file_name_col.append(fastq)
            flow_cell_id_col.append(flow_cell_id)
            index_sequence_col.append(index_sequence)
            lane_number_col.append(lane_number)
            read_number_col.append(read_number)
            sample_name_col.append(sample_name)
            library_type_col.append(library_type)
            total_bases_col.append(total_bases)
            samplemap_path_col.append(samplemap_path)
            total_reads_col.append(total_reads)
            esp_id_col.append(esp_id)
            pool_name_col.append(pool_name)

        field_lengths: set = set()
        field_lengths.add(len(file_name_col))
        field_lengths.add(len(flow_cell_id_col))
        field_lengths.add(len(index_sequence_col))
        field_lengths.add(len(lane_number_col))
        field_lengths.add(len(read_number_col))
        field_lengths.add(len(sample_name_col))
        field_lengths.add(len(library_type_col))
        field_lengths.add(len(total_bases_col))
        field_lengths.add(len(samplemap_path_col))
        field_lengths.add(len(total_reads_col))
        field_lengths.add(len(esp_id_col))
        field_lengths.add(len(pool_name_col))
        if len(field_lengths) > 1:
            raise ValueError(
                'Dataframe field row lengths not of equal length.',
                field_lengths
            )

        df_subset = DataFrame()
        df_subset['fastq'] = file_name_col
        df_subset['flow_cell_id'] = flow_cell_id_col
        df_subset['index_sequence'] = index_sequence_col
        df_subset['lane_number'] = lane_number_col
        df_subset['read_number'] = read_number_col
        df_subset['sample_name'] = sample_name_col
        df_subset['library_type'] = library_type_col
        df_subset['total_bases'] = total_bases_col
        df_subset['samplemap_path'] = samplemap_path_col
        df_subset['total_reads'] = total_reads_col
        df_subset['esp_id'] = esp_id_col
        df_subset['pool_name'] = pool_name_col
        return df_subset

    def __parse_samplemap(self: Self, i: int, smap: str,
                          smap_type: str) -> None:
        """Parse a Samplemap file and add to the object instance.

        We will parse a given Samplemap.csv file and convert the
        information to a uniform, subset of values in dataframe format.
        The dataframe will be added to a running list of batch
        dataframes in the object instance. Order will be retained from
        the order of samplemaps called from the command line interface.

        Parameters
        ----------
        i: int
            Samplemap "batch" number/order from the command line
            arguments. Order will be retained.

        smap : str
            A path to a Samplemap.csv file.

        smap_type : str
            The identified Samplemap format type, as determined by the
            __type_samplemap_format() method.

        Raises
        ------
        TypeError : Unknown samplemap format type encountered.
        """
        if smap_type == 'smap_mid_2024_format':
            df_subset = self.__smap_mid_2024_to_df(smap=smap)
            self.smaps.update({
                i: {
                    'samplemap_type': smap_type,
                    'samplemap_path': str(Path(smap).as_posix()),
                    'fastq_dir': str(Path(smap).parent),
                    'df': df_subset,
                }
            })
        else:
            raise TypeError(
                'Unknown samplemap format type encountered.',
                smap_type
            )
        return

    def __parse_samplemaps(self: Self) -> None:
        """Parse all of the input Samplemap files."""
        for i, smap in enumerate(self.args.samplemap, 1):
            smap_type = self.__type_samplemap_format(smap=smap)
            self.__parse_samplemap(i=i, smap=smap, smap_type=smap_type)
            # self.__eval_cross_batch_sample_ids()
            # self.__concatenate_samplemaps()
        return

# __END__
