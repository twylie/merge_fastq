# Project     : merge_fastq
# File Name   : merge_fastq.py
# Description : Functions for merging FASTQ files.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Thu Sep 19 17:01:30 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

import mergefastq  # type: ignore
import argparse
from . rename_samples import RenameSamples  # type: ignore
from . samplemap import Samplemap  # type: ignore
from typing_extensions import Self
from pathlib import Path
import hashlib
import gzip
from datetime import datetime


class MergeFastq:
    """A class for merging split FASTQ provided by GTAC@MGI.

    This class handles merging split FASTQ files as provided provided by
    GTAC@MGI. We will also (potentially) rename the samples once merged.
    Sample mappings are based on the original Samplemap.csv files
    provided alongside the FASTQ files. Every unique sample will resolve
    to a pair (R1 & R2) of FASTQ files. We will also manually perform an
    independent read count review of each FASTQ file. Merging/counting
    commands will be executed in LSF jobs using compute1 at WashU.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments for merging FASTQ files as provided by argparse.

    rename : RenameSamples
        A RenameSamples object from the the mergefastq package.

    samplemap : Samplemap
        A post-concatenation Samplemap object from the mergefastq
        package.

    Attributes
    ----------
    copy_cmds : dict
        A dictionary of sample names and associated (lists) of shell
        commands for "copy" type FASTQ.

    dest_fq_index : dict
        A dictionary of sample names and associated pairs of R1 & R2
        destination file name paths.

    log_dir_path : Path
        A pathlib.Path object for the path to write the LSF shell
        commands, etc.

    merge_cmds : dict
        A dictionary of sample names and associated (lists) of shell
        commands for "merge" type FASTQ.

    merge_copy_ids : set
        The full, unique set of sample names that are of "merge" FASTQ
        type.

    all_jobs : list
        A list of Bsub objects for running the commands for "merge" and
        "copy" type FASTQ commands.

    sample_dir : dict
        A dictionary of sample names and associated output directories
        to write merged FASTQ files. Parent directories are "old" sample
        names whereas the merged files will use "new" sample names.

    samplemap_merged : DataFrame
        A post-concatenation samplemap object as provide by the
        Samplemap class.

    single_copy_ids : set
        The full, unique set of sample names that are of "copy" FASTQ
        type.

    Methods
    -------
    launch_lsf_jobs()
        Launch the copy and merge LSF jobs.

    prepare_lsf_cmds()
        Prepare the LSF jobs for merging FASTQ commands.

    setup_output_dirs()
        Setup all of the sample output directories.

    write_df()
        Write the merged FASTQ dataframe to a tab-delimited file.

    Examples
    --------
    merge_fastq = mergefastq.MergeFastq(
        args=args,
        rename=rename_samples,
        samplemap=samplemap
    )
    merge_fastq.setup_output_dirs()
    merge_fastq.prepare_lsf_cmds()
    merge_fastq.launch_lsf_jobs()
    """

    def __init__(self: Self, args: argparse.Namespace, rename: RenameSamples,
                 samplemap: Samplemap) -> None:
        """Construct the class.

        Parameters
        ----------
        args : argparse.Namespace
            Arguments for merging FASTQ files as provided by argparse.

        rename : RenameSamples
            A RenameSamples object from the the mergefastq package.

        samplemap : Samplemap
            A post-concatenation Samplemap object from the mergefastq
            package.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        self.args = args
        self.rename = rename
        self.samplemap = samplemap
        self.samplemap_merged = self.samplemap.copy_df()
        self.dest_fq_index: dict = dict()
        self.copy_cmds: dict = dict()
        self.merge_cmds: dict = dict()
        self.log_dir_path: Path = Path()
        self.sample_dir: dict = dict()
        self.all_jobs: list = list()
        self.single_copy_ids: set = set()
        self.merge_copy_ids: set = set()
        self.__parse_fastq_copy_types()
        self.__setup_copy_cmds()
        self.__setup_merge_cmds()
        self.__update_df_cmds()
        self.__update_df_dest_fq()
        self.__update_df_read_counts()
        return

    def __parse_fastq_copy_types(self: Self) -> None:
        """Parse the Samplemap dataframe and determine FASTQ copy types.

        We may determine the copy types for samples based on the
        cardinality of sample ids to flowcells. For example, if a unique
        sample has 2 FASTQ (R1 & R2) and 1 associated flowcell, we don't
        need to merge the FASTQ, but rather simply copy them using the
        original sample name. Samples with more these attributes will
        require merging of FASTQ files. Therefore, types fall into
        "copy" or "merge" categories.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            FASTQ copy type counts differ.

        Returns
        -------
        None
        """
        df = self.samplemap.copy_df()
        dfg = df.groupby('sample_name')
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
        # i = total unique sample names
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

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Copy type should only have 2 associated FASTQ files.

        ValueError
            R1 and R2 index sequences do not match.

        ValueError
            R1 and R2 flow cell ids do not match.

        ValueError
            R1 and R2 lane numbers do not match.

        ValueError
            Sample name should not be duplicated in FASTQ index.

        ValueError
            FASTQ R1 comp or decomp file not found.

        ValueError
            FASTQ R2 comp or decomp file not found.

        Returns
        -------
        None
        """
        df_smaps = self.samplemap.copy_df()
        dfg_ids = df_smaps.groupby('sample_name')
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

            if sample_name not in self.dest_fq_index.keys():
                self.dest_fq_index[sample_name] = {
                    'R1': str(dest_r1_fq.resolve()),
                    'R2': str(dest_r2_fq.resolve())
                }
            else:
                raise ValueError(
                    'Sample name should not be duplicated in FASTQ index.',
                    sample_name
                )

            if Path(src_r1_fq).is_file() is True:
                src_r1_fq_eval = src_r1_fq
            elif Path(src_r1_fq).is_file() is False:
                if Path(src_r1_fq[:-3]).is_file() is True:
                    src_r1_fq_eval = src_r1_fq[:-3]
                elif Path(src_r1_fq[:-3]).is_file() is False:
                    raise FileNotFoundError(
                        'FASTQ R1 comp or decomp file not found.',
                        src_r1_fq
                    )

            if Path(src_r2_fq).is_file() is True:
                src_r2_fq_eval = src_r2_fq
            elif Path(src_r2_fq).is_file() is False:
                if Path(src_r2_fq[:-3]).is_file() is True:
                    src_r2_fq_eval = src_r2_fq[:-3]
                elif Path(src_r2_fq[:-3]).is_file() is False:
                    raise FileNotFoundError(
                        'FASTQ R2 comp or decomp file not found.',
                        src_r2_fq
                    )

            with gzip.open(src_r1_fq_eval, 'r') as fhi:
                try:
                    fhi.read(1)
                    is_src_r1_fq_gzip = True
                except gzip.BadGzipFile:
                    is_src_r1_fq_gzip = False

            with gzip.open(src_r2_fq_eval, 'r') as fhi:
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
                r1_copy_cmds.append(
                    f'cp {src_r1_fq} {str(dest_r1_fq.resolve())}'
                )
                r1_copy_cmds.append(
                    f'md5sum {str(dest_r1_fq.resolve())} > '
                    f'{str(dest_r1_fq.resolve())}.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {str(dest_r1_fq.resolve())}) / 4 '
                       f'| bc > {str(dest_r1_fq.resolve())}.counts')
                r1_copy_cmds.append(cmd)
            elif is_src_r1_fq_gzip is False:
                src_r1_fq_stem = src_r1_fq[:-3]
                dest_r1_fq_stem = str(dest_r1_fq.resolve())[:-3]
                r1_copy_cmds.append(f'cp {src_r1_fq_stem} {dest_r1_fq_stem}')
                r1_copy_cmds.append(f'gzip {dest_r1_fq_stem}')
                r1_copy_cmds.append(
                    f'md5sum {dest_r1_fq_stem}.gz > {dest_r1_fq_stem}.gz.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {dest_r1_fq_stem}.gz) / 4 | bc '
                       f'> {dest_r1_fq_stem}.gz.counts')
                r1_copy_cmds.append(cmd)

            r2_copy_cmds: list = list()
            if is_src_r2_fq_gzip is True:
                r2_copy_cmds.append(
                    f'cp {src_r2_fq} {str(dest_r2_fq.resolve())}'
                )
                r2_copy_cmds.append(
                    f'md5sum {str(dest_r2_fq.resolve())} > '
                    f'{str(dest_r2_fq.resolve())}.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {str(dest_r2_fq.resolve())}) / 4 '
                       f'| bc > {str(dest_r2_fq.resolve())}.counts')
                r2_copy_cmds.append(cmd)
            elif is_src_r2_fq_gzip is False:
                src_r2_fq_stem = src_r2_fq[:-3]
                dest_r2_fq_stem = str(dest_r2_fq.resolve())[:-3]
                r2_copy_cmds.append(f'cp {src_r2_fq_stem} {dest_r2_fq_stem}')
                r2_copy_cmds.append(f'gzip {dest_r2_fq_stem}')
                r2_copy_cmds.append(
                    f'md5sum {dest_r2_fq_stem}.gz > {dest_r2_fq_stem}.gz.MD5'
                )
                cmd = (f'echo $(zgrep -c ^ {dest_r2_fq_stem}.gz) / 4 | bc '
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

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Merge-samples should all have the same index sequence.

        ValueError
            Merge-samples sort order was not maintained.

        ValueError
            Merge-samples revised sample names are not equal.

        ValueError
            Merge-samples revised names must be unique.

        ValueError
            FASTQ R1 comp or decomp file not found.

        ValueError
            FASTQ R2 comp or decomp file not found.
        """
        df_smaps = self.samplemap.copy_df()
        dfg_ids = df_smaps.groupby('sample_name')
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

            if sample_name not in self.dest_fq_index.keys():
                self.dest_fq_index[sample_name] = {
                    'R1': str(dest_r1_fq.resolve()),
                    'R2': str(dest_r2_fq.resolve())
                }
            else:
                raise ValueError(
                    'Sample name should not be duplicated in FASTQ index.',
                    sample_name
                )

            r1_is_gzip: list = list()
            for r1_fq in src_r1_fq:
                if Path(r1_fq).is_file() is True:
                    r1_fq_eval = r1_fq
                elif Path(r1_fq).is_file() is False:
                    if Path(r1_fq[:-3]).is_file() is True:
                        r1_fq_eval = r1_fq[:-3]
                    elif Path(r1_fq[:-3]).is_file() is False:
                        raise FileNotFoundError(
                            'FASTQ R1 comp or decomp file not found.',
                            r1_fq
                        )
                with gzip.open(r1_fq_eval, 'r') as fhi:
                    try:
                        fhi.read(1)
                        is_src_r1_fq_gzip = True
                    except gzip.BadGzipFile:
                        is_src_r1_fq_gzip = False
                r1_is_gzip.append(is_src_r1_fq_gzip)

            r2_is_gzip: list = list()
            for r2_fq in src_r2_fq:
                if Path(r2_fq).is_file() is True:
                    r2_fq_eval = r2_fq
                elif Path(r2_fq).is_file() is False:
                    if Path(r2_fq[:-3]).is_file() is True:
                        r2_fq_eval = r2_fq[:-3]
                    elif Path(r2_fq[:-3]).is_file() is False:
                        raise FileNotFoundError(
                            'FASTQ R2 comp or decomp file not found.',
                            r2_fq
                        )
                with gzip.open(r2_fq_eval, 'r') as fhi:
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
                    r1_merge_cmds.append(
                        f'cp {from_r1_fq} {str(to_r1_fq.resolve())}'
                    )
                    r1_merge_cmds.append(
                        f'gzip {str(to_r1_fq.resolve())}'
                    )
                    to_r1_gzip_fq = str(to_r1_fq) + '.gz'
                    r1_cat_order.append(to_r1_gzip_fq)
            cat_files = ' '.join(r1_cat_order)
            r1_merge_cmds.append(
                f'cat {cat_files} > {str(dest_r1_fq.resolve())}'
            )
            r1_merge_cmds.append(
                f'md5sum {str(dest_r1_fq.resolve())} > '
                f'{str(dest_r1_fq.resolve())}.MD5'
            )
            cmd = (f'echo $(zgrep -c ^ {str(dest_r1_fq.resolve())}) / 4 '
                   f'| bc > {str(dest_r1_fq.resolve())}.counts')
            r1_merge_cmds.append(cmd)
            if tmp_r1:
                for tmp_file in tmp_r1:
                    rm_file = str(tmp_file.resolve()) + '.gz'
                    r1_merge_cmds.append(f'rm {rm_file}')

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
                    r2_merge_cmds.append(
                        f'cp {from_r2_fq} {str(to_r2_fq.resolve())}'
                    )
                    r2_merge_cmds.append(
                        f'gzip {str(to_r2_fq.resolve())}'
                    )
                    to_r2_gzip_fq = str(to_r2_fq) + '.gz'
                    r2_cat_order.append(to_r2_gzip_fq)
            cat_files = ' '.join(r2_cat_order)
            r2_merge_cmds.append(
                f'cat {cat_files} > {str(dest_r2_fq.resolve())}'
            )
            r2_merge_cmds.append(
                f'md5sum {str(dest_r2_fq.resolve())} > '
                f'{str(dest_r2_fq.resolve())}.MD5'
            )
            cmd = (f'echo $(zgrep -c ^ {str(dest_r2_fq.resolve())}) / 4 '
                   f'| bc > {str(dest_r2_fq.resolve())}.counts')
            r2_merge_cmds.append(cmd)
            if tmp_r2:
                for tmp_file in tmp_r2:
                    rm_file = str(tmp_file.resolve()) + '.gz'
                    r2_merge_cmds.append(f'rm {rm_file}')

            self.merge_cmds.update({
                sample_name: (r1_merge_cmds, r2_merge_cmds)
            })
        return

    def setup_output_dirs(self: Self) -> None:
        """Setup all of the sample output directories.

        We are setting up the output directory, the bsub log directory,
        and all of the per-sample write-directories. None of these
        should exist at this point, so any preexisting directories will
        throw an exception. We will update the MergeFastq object with
        the directories that we create.

        Parameters
        ----------
        None

        Raises
        ------
        IsADirectoryError
            The outdir directory already exists.

        Returns
        -------
        None
        """
        if Path(self.args.outdir).is_dir() is True:
            raise IsADirectoryError(
                'The outdir directory already exists.',
                self.args.outdir
            )
        else:
            Path(self.args.outdir).mkdir(parents=False, exist_ok=False)

        log_dir = Path(self.args.outdir) / '__bsub'
        self.log_dir_path = log_dir

        for sample_name in self.copy_cmds:
            sample_dir = Path(self.args.outdir) / Path(sample_name)
            sample_dir.mkdir(parents=False, exist_ok=False)
            self.sample_dir.update({
                sample_name: str(sample_dir.resolve())
            })

        for sample_name in self.merge_cmds:
            sample_dir = Path(self.args.outdir) / Path(sample_name)
            sample_dir.mkdir(parents=False, exist_ok=False)
            self.sample_dir.update({
                sample_name: str(sample_dir.resolve())
            })
        return

    def prepare_lsf_cmds(self: Self) -> None:
        """Prepare the LSF jobs for merging FASTQ commands.

        Each unique sample gets its own set of shell commands and
        associated LSF jobs information. We will be using the Bsub class
        to formulate LSF jobs to run at WashU.

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
        log_dir = self.log_dir_path
        lsf_vols: dict = dict()
        for lsf_vol in self.args.lsf_vol:
            if lsf_vol.endswith('/') is True:
                lsf_vols.update({lsf_vol[:-1]: lsf_vol[:-1]})
            else:
                lsf_vols.update({lsf_vol: lsf_vol})
        all_cmds = self.copy_cmds | self.merge_cmds
        for i, sample_name in enumerate(all_cmds, 1):
            r1_cmds, r2_cmds = all_cmds[sample_name]
            cmds = r1_cmds + r2_cmds
            error_log = f'{i}_merge_fastq_bsub.err'
            cmd_name = f'{i}_merge_fastq.sh'
            out_log = f'{i}_merge_fastq_bsub.out'
            job_name = f'{i}_merge_fastq'
            config = f'{i}_merge_fastq_bsub.yaml'
            bsub_cmd = f'{i}_merge_fastq_bsub.sh'
            bsub_log_dir = str(log_dir.resolve())
            lsf_job = mergefastq.Bsub(
                log_dir=bsub_log_dir,
                docker_volumes=lsf_vols,
                docker_image=self.args.lsf_image,
                group=self.args.lsf_group,
                queue=self.args.lsf_queue,
                command=cmds,
                error_log=error_log,
                output_log=out_log,
                command_name=cmd_name,
                config=config,
                bsub_command_name=bsub_cmd,
                job_name=job_name
            )
            self.all_jobs.append(lsf_job)
        return

    def launch_lsf_jobs(self: Self) -> None:
        """Launch the copy and merge LSF jobs.

        This method will launch the LSF bsub jobs that are defined in
        the MergeFastq object. Both the copy and merge jobs will be
        launched, in succession. If the --lsf-dry argument is passed,
        all of the bsub jobs will be written to output, but they will
        not be launched--i.e. a "dry run" of the pipeline; else, all
        jobs will be launched and submitted to the LSF queue.

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
        for job in self.all_jobs:
            job.execute(dry=self.args.lsf_dry)
        return

    def __update_df_dest_fq(self: Self) -> None:
        """Update the destination FASTQ column in dataframe.

        This method will add the destination FASTQ paths to the existing
        merged dataframe. These paths are to the merged FASTQ files (R1
        & R2) for a unique sample name. Destination sample names reflect
        the RenameSamples id mappings.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Sample name not in the destination FASTQ index.

        ValueError
            Destination FASTQ indexes lengths differ.

        Returns
        -------
        None
        """
        col_dest_fq_path: list = list()
        for i in self.samplemap_merged.index:
            df_i = self.samplemap_merged.loc[i]
            sample_name = df_i['sample_name']
            read_num = df_i['read_number']
            if sample_name in self.dest_fq_index.keys():
                if read_num == 1:
                    col_dest_fq_path.append(
                        self.dest_fq_index[sample_name]['R1']
                    )
                elif read_num == 2:
                    col_dest_fq_path.append(
                        self.dest_fq_index[sample_name]['R2']
                    )
            else:
                raise ValueError(
                    'Sample name not in the destination FASTQ index.',
                    sample_name
                )
        if len(col_dest_fq_path) != len(self.samplemap_merged.index):
            raise ValueError('Destination FASTQ indexes lengths differ.')
        else:
            self.samplemap_merged['merged_fastq_path'] = col_dest_fq_path
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
        str
            Returns a MD5 hash value for the supplied file.
        """
        with open(file_path, 'rb') as fh:
            data = fh.read()
            md5sum = hashlib.md5(data).hexdigest()
        return md5sum

    def write_df(self: Self, file_path: str) -> None:
        """Write the merged FASTQ dataframe to a tab-delimited file.

        For reference, we will write the dataframe to disk as both plain
        text and as a binary (pickle) file. The pickle file retains data
        objects in the dataframe, whereas the text version does not.

        Output Files
        ------------
        file.tsv
        file.tsv.MD5
        file.tsv.pickle

        Parameters
        ----------
        file_path : str
            A qualified file path to write the merged dataframe.

        Raises
        ------
        None

        Returns
        -------
        None
        """
        header = self.__format_tsv_header(file_path=file_path)
        with open(file_path, 'w') as fho:
            fho.write(header)
            self.samplemap_merged.to_csv(fho, sep='\t', index=False)
            fho.write('# End file text.\n')
        md5 = self.__calc_file_md5(file_path=file_path)
        md5_path = file_path + '.MD5'
        with open(md5_path, 'w') as fh:
            fh.write(f'MD5 ({file_path}) = {md5}\n')
        pickle = file_path + '.pickle'
        self.samplemap_merged.to_pickle(pickle)
        self.samplemap_merged_pkl = pickle
        return

    def __format_tsv_header(self: Self, file_path) -> str:
        """Format a header for the dataframe tsv file.

        Parameters
        ----------
        file_path : str
            A qualified file path to write the merged dataframe.

        Raises
        ------
        None

        Returns
        -------
        str
            Returns the header information block as a formatted string.
        """
        file_path_name = Path(file_path).name
        date = datetime.now().strftime("%a %d. %b %H:%M:%S %Z %Y")
        from textwrap import dedent
        project = self.args.project
        header = dedent(f"""\
        # Begin file text.
        #
        # Project     : {project}
        # File Name   : {file_path_name}
        # Description : A dataframe of merged FASTQ file sequencing
        #               information based on Samplemap.csv files.
        # Created     : {date}
        #
        # Fields
        # ------
        # fastq : STR
        #     The name of the FASTQ file. Will be either R1 or R2 type.
        #
        # flow_cell_id : STR
        #     Sequencing flow cell id from which the sample's FASTQ file was
        #     derived.
        #
        # index_sequence : STR
        #     Molecular barcode used to identify the sequenced sample.
        #
        # lane_number : INT
        #     Flow cell lane in which the sample was sequenced.
        #
        # read_number : INT
        #     Read-end, R1 or R2, for the sequenced read-pair.
        #
        # sample_name : STR
        #     Original sample name as provided by the sequencing core.
        #
        # library_type : STR
        #     Library type as provided by the sequencing core.
        #
        # total_bases : INT
        #     Summed base count using both read-pairs in the sample, as
        #     provided by the sequencing core.
        #
        # samplemap_path : STR
        #     The Samplemap.csv file from which the sample was taken.
        #
        # gtac_fastq_reads : INT
        #     Read count for the FASTQ file provided by the sequencing core.
        #
        # esp_id : STR
        #     ESP ID is an auto-generated ID from the sequencing core's LIMS
        #     system; each entity type has one.
        #
        # pool_name : STR
        #     Sequencing pool identifier.
        #
        # batch_id : INT
        #     The sequencing batch as determined by individual Samplemap.csv
        #     files.
        #
        # fastq_path : STR
        #     File path to the original, source FASTQ file used in merging.
        #
        # project : STR
        #     The project name/tag associated with the FASTQ file.
        #
        # revised_sample_name : STR
        #     A revised sample name used in merging and downstream analysis.
        #
        # merged_commands : STR
        #     The shell commands used to create the associated merged
        #     FASTQ file.
        #
        # merged_fastq_path : STR
        #     File path to the merged FASTQ file.
        #
        # gtac_end_pair_reads : INT
        #     Summed read count using all end-pairs in a FASTQ merge, as
        #     provided by the sequencing core.
        #
        # gtac_sample_reads : INT
        #     Summed read count using both read-pairs in the sample, as
        #     provided by the sequencing core.
        #
        """)
        return header

    def __update_df_read_counts(self: Self) -> None:
        """Update the GTAC supplied read counts in the dataframe.

        We will be tracking the read counts as supplied by the GTAC@MGI
        sequencing core. Of interest are the per-FASTQ read counts, the
        merged end pair (R1 & R2) counts, and the total number of reads
        per sample. All of these values are calculated based on the
        original Total Reads values supplied by GTAC@MGI in the
        Samplemap files.

        IMPORTANT: The read counts provided by this method are not the
        manual read counts performed in the MergeFastq LSF jobs.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            End pair counts column length does not match dataframe.

        ValueError
            Sample counts column length does not match dataframe.

        ValueError
            GTAC sample counts do not match for R1 & R2.

        ValueError
            GTAC end pair counts do not match for R1 & R2.

        Returns
        -------
        None
        """
        count_index: dict = dict()
        df = self.samplemap_merged.copy()
        dfg = df.groupby(['sample_name', 'read_number'])
        for i in dfg.groups:
            sample_name, read_number = i
            dfi = dfg.get_group(i)
            merged_count = dfi.sum()['gtac_fastq_reads']
            count_index[sample_name] = {'merged_count': merged_count}

        col_end_pair_counts: list = list()
        col_sample_counts: list = list()
        for i in df.index:
            sample_name = df.loc[i]['sample_name']
            merged_reads = count_index[sample_name]['merged_count']
            all_reads = merged_reads * 2
            col_end_pair_counts.append(merged_reads)
            col_sample_counts.append(all_reads)

        if len(col_end_pair_counts) != len(self.samplemap_merged.index):
            raise ValueError(
                'End pair counts column length does not match dataframe.'
            )

        if len(col_sample_counts) != len(self.samplemap_merged.index):
            raise ValueError(
                'Sample counts column length does not match dataframe.'
            )

        self.samplemap_merged['gtac_end_pair_reads'] = col_end_pair_counts
        self.samplemap_merged['gtac_sample_reads'] = col_sample_counts

        dfg = self.samplemap_merged.groupby('revised_sample_name')
        for sample_name in dfg.groups:
            dfi = dfg.get_group(sample_name)
            dfg_r1_fastq = dfi.groupby('read_number').get_group(1)
            r1_sample_counts = list(set(dfg_r1_fastq['gtac_sample_reads']))[0]
            r1_ep_counts = list(set(dfg_r1_fastq['gtac_end_pair_reads']))[0]
            dfg_r2_fastq = dfi.groupby('read_number').get_group(2)
            r2_sample_counts = list(set(dfg_r2_fastq['gtac_sample_reads']))[0]
            r2_ep_counts = list(set(dfg_r2_fastq['gtac_end_pair_reads']))[0]
            if r1_sample_counts != r2_sample_counts:
                raise ValueError(
                    'GTAC sample counts do not match for R1 & R2.',
                    sample_name,
                    f'R1={r1_sample_counts}',
                    f'R2={r2_sample_counts}'
                )
            if r1_ep_counts != r2_ep_counts:
                raise ValueError(
                    'GTAC end pair counts do not match for R1 & R2.',
                    sample_name,
                    f'R1={r1_ep_counts}',
                    f'R2={r2_ep_counts}'
                )
        return

    def __update_df_cmds(self: Self) -> None:
        """Update the merged dataframe with merging commands.

        We are storing the per-FASTQ merging commands as executed by the
        Bsub class, for reference. These entries are useful for sanity
        checks and archiving how the merged FASTQ files were created.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            Unknown read number.

        KeyError
            Sample name missing in copy commands keys.

        Returns
        -------
        None
        """
        df_merged = self.samplemap_merged.copy()
        col_merge_cmds: list = list()
        for i in df_merged.index:
            df = df_merged.loc[i]
            sample_name = df['sample_name']
            read_number = df['read_number']
            if sample_name in self.copy_cmds.keys():
                r1_cmd, r2_cmd = self.copy_cmds[sample_name]
                if read_number == 1:
                    col_merge_cmds.append(r1_cmd)
                elif read_number == 2:
                    col_merge_cmds.append(r2_cmd)
                else:
                    raise ValueError(
                        'Unknown read number.',
                        read_number
                    )
            elif sample_name in self.merge_cmds.keys():
                r1_cmd, r2_cmd = self.merge_cmds[sample_name]
                if read_number == 1:
                    col_merge_cmds.append(r1_cmd)
                elif read_number == 2:
                    col_merge_cmds.append(r2_cmd)
                else:
                    raise ValueError(
                        'Unknown read number.',
                        read_number
                    )
            else:
                raise KeyError(
                    'Sample name missing in copy commands keys.',
                    sample_name
                )
        self.samplemap_merged['merged_commands'] = col_merge_cmds
        return

# __END__
