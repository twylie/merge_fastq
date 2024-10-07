<!-- Begin document text. -->
<!-- T.N. Wylie -->
<!-- twylie@wustl.edu -->
<!-- version 1.0.0  -->
<!-- Sun Oct  6 20:10:26 CDT 2024 -->

# merge_fastq

## Overview

### What does it do?

The `merge_fastq` code base and associated commands may be used to evaluate and merge split FASTQ files. Sequencing FASTQ files generate by GTAC@MGI Sequencing Core can potentially be split across multiple lanes for a single, unique sample. We need to merge the split FASTQ files into a single, unified read pair (R1 and R2 FASTQ files) for a unique sample. This code also allows for renaming GTAC sample ids to any arbitrary naming convention during the FASTQ merging process. Finally, sample-level read count assessment is performed across a range of predefined target read count requirements—e.g. >= 80% of 30,000,000 read pairs for passing samples.

### How does it do it?

A user can merge FASTQ files using a few simple commands from a terminal session on `compute1` at WashU. The merging command automatically launches parallel `bsub` jobs for all samples being evaluated using the IBM LSF system—i.e. one sample equals one `bsub` job. There is no direct manipulation of the source FASTQ files other than copying or reading them; no source FASTQ files are ever altered or deleted by `merge_fastq` commands. The merging process also creates consolidated dataframe files which store important information related to all the FASTQ files being evaluated. A Docker image supplies all of the required code and dependencies needed to merge FASTQ files at WashU.

### What prerequisites are required?

There are several dependencies required in order to effectively run `merge_fastq` commands. They are described below.

+ **Docker Image** — You will need access to the Docker image for the latest `merge_fastq` code. The commands must be executed using the supported Docker image. During the FASTQ merging step, each sample will receive its own `bsub` job on WashU's `compute1` LSF queue, which requires execution using a Docker container.
+ **Samplemap.csv Files** — MGI@GTAC provides a metadata file for each batch of sequencing FASTQ files it provides, named `Samplemap.csv`. You will need at least one `Samplemap.csv` and associated FASTQ files to perform a merge.
+ **FASTQ Files** — MGI@GTAC provides compressed (gzipped) FASTQ files, one for each end-pair (R1 and R2). You will need these FASTQ files to perform a merge.
+ **Resources** — You will need an account and access to `compute1` at WashU to run this code. Adequate disk space is required, commensurate with the size of the original sequencing batch being evaluated.

## Steps for Merging FASTQ

>[!note] Conceptual Overview
> This section outlines the conceptual steps in merging a FASTQ batch more so than command line instructions. I recommend reading over this section to understand the steps required before attempting to run related commands. For detailed instructions on command usage, see the **Commands** and **Test Data Set and Demonstration** sections below. Sitting down and running the commands in **Running Test Data Demo** provides the best step-by-step instructions for command line formation.

### 1. Locate Sequencing Batches

You will need, minimally, a batch of sequencing FASTQ files from GTAC@MGI, accompanied by a `Samplemap.csv` metadata file. All of the FASTQ files should be at the same level in a directory and the `Samplemap.csv` must be at the same level. As of early 2024, GTAC@MGI provides two samplemap files: 1) `Samplemap.csv`; 2) `Samplemap2.csv`. We will only be using the `Samplemap.csv` formatted file.

The FASTQ files will be listed in the `Samplemap.csv` as being compressed; however, the `merge_fastq` code handles compressed (gzipped), uncompressed (plaintext), or a mixture. For example, `sample2.A_R1_001.fastq` and `sample2.A_R2_001.fastq.gz` are both viable.

The concept of a *batch* here is arbitrary. The only significance is that a batch must be a single directory with one set of FASTQ files and single, representative `Samplemap.csv` file, all at the same level. The following example shows two batches.

```plaintext
batch1
├── Samplemap.csv
├── sample1.A_R1_001.fastq.gz
├── sample1.A_R2_001.fastq.gz
├── sample1.B_R1_001.fastq.gz
└── sample1.B_R2_001.fastq.gz
batch2
├── Samplemap.csv
├── sample2.A_R1_001.fastq
├── sample2.A_R2_001.fastq.gz
├── sample2.B_R1_001.fastq.gz
├── sample2.B_R2_001.fastq
├── sample3_R1_001.fastq
└── sample3_R2_001.fastq

2 directories, 12 files
```

We may process multiple batches at one time, provided that the batches follow the rules listed above.

As of writing (2024-10-04) the GTAC@MGI `Samplemap.csv` format looks like the following.

```plaintext
Samplemap.csv Fields
--------------------
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
```

```python
## Example dictionary for a FASTQ entry in Samplemap.csv file.
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
```

If the format of the `Samplemap.csv` file changes in the future, the `merge_fastq` code base will require updates to work properly.

> [!IMPORTANT] Summary
> + A batch is a set of sequencing FASTQ files accompanied by a `Samplemap.csv` file.
> + All of the FASTQ files and `Samplemap.csv` should be at the same level in a directory.
> + Use `Samplemap.csv` files and not `Samplemap2.csv` files.
> + Any `Samplemap.csv` file format that does not match the early 2024 version will break the code base.

###  2. Prepare a Rename File

It is not uncommon to need to rename the FASTQ files in a batch. GTAC@MGI often names the samples within the `Samplemap.csv` in a manner different than what may be desired for a sample's final, canonical name. Also, a sample may have been sent to sequencing with a name that later needs revision. For example, the following samples need conversion from the samplemap naming convention to a final, revised naming convention.

```plaintext
Samplemap.csv Sample Names
--------------------------
1. WLAB-H141_10_12_2-1
2. WLAB-H147_3_13_2-1-1
3. WLAB-H147_3_13_2-2-1
4. WLAB-H141_10_13_2_Swift-1
5. WLAB-H149_1_13_2_Swift-2-1
6. WLAB-H152_11_15_2-H152_11_15_2

Merged FASTQ Sample Names
-------------------------
1. H141_10_12_2
2. H147_3_13_2
3. H147_3_13_2
4. H141_10_13_2
5. H149_1_13_2
6. H152_11_15_2
```

We need a mechanism to map on-file sequence sample names to revised sample names, with a one-to-one mapping. Running the `prep_rename_file` command provides this functionality. You will need to log into WashU's `compute1` network to run the command. Running the command produces a `rename.tsv` file as a starting point, based on all of the `Samplemap.csv` files used as input. A user can manually edit this file as desired to change any of the sample renaming mappings, before supplying the `rename.tsv` file downstream to the `merge_fastq` command.

If no sample renaming is required, the `samplemap_sample_id` and `revised_sample_id` column values should remain equal. The free-form `comments` field is for the user's reference and is never evaluated by any of the `merge_fastq` commands.

Technically, this step may be omitted provided a user manually creates a `rename.tsv` which adheres to the format specification—i.e. tab-delimited, provides the required column types, and all sequencing sample names in the supplied `Samplemap.csv` files are unique and present. The `rename.tsv` file name is arbitrary and any name may be used for the resultant file produced by the `prep_rename_file` command.

Example tab-delimited `rename.tsv` file:

```plaintext
samplemap_sample_id             revised_sample_id  comments
WLAB-H141_10_12_2-1             H141_10_12_2
WLAB-H147_3_13_2-1-1            H147_3_13_2
WLAB-H147_3_13_2-2-1            H147_3_13_2
WLAB-H141_10_13_2_Swift-1       H141_10_13_2       Swift kit test sample.
WLAB-H149_1_13_2_Swift-2-1      H149_1_13_2        Swift kit test sample.
WLAB-H152_11_15_2-H152_11_15_2  H152_11_15_2
```

>[!IMPORTANT] Summary
> + It is not uncommon to need to rename the FASTQ files in a batch.
> + Running the `prep_rename_file` command provides a initial file for renaming samples.
> + A user can manually edit a `rename.tsv` file before supplying it to the `merge_fastq` command.
> + The `samplemap_sample_id` and `revised_sample_id` values should be the same if there is no sample name change.

### 3. Merge the FASTQ Files

Once you have located the FASTQ batches you wish to merge (Step 1) and have created an associated `rename.tsv` file (Step 2), you are ready to run the `merge_fastq` command. This command will launch all of the parallel jobs required to appropriately merge the FASTQ files and write related metadata. You will need to log into WashU's `compute1` network to run the command.

Running `merge_fastq` command with the following input files produced the output listed below.

**INPUT**

**BATCH 1:** `Samplemap.csv`

```plaintext
FASTQ,Flowcell ID,Index Sequence,Flowcell Lane,ESP ID,Pool Name,Species,Illumina Sample Type,Library Type,Library Name,Date Complete,Total Reads,Total Ba
ses,% >Q30,% Pass Filter Clusters,PhiX Error Rate,Avg Q Score
sample1.A_R1_001.fastq.gz,227NFFLT4,GTACAGCGGA-CAATCTGTGT,4,LIB043956-DIL01,,human,,WGS,H214_1_3_2,2024-06-27 14:22:09.194473+00:00,"100","30,200",91.0,5
6.98,0.39,38.13
sample1.A_R2_001.fastq.gz,227NFFLT4,GTACAGCGGA-CAATCTGTGT,4,LIB043956-DIL01,,human,,WGS,H214_1_3_2,2024-06-27 14:22:09.194473+00:00,"100","30,200",87.0,5
6.98,0.48,37.41
sample1.B_R1_001.fastq.gz,22CGWFLT4,GTACAGCGGA-CAATCTGTGT,6,LIB043956-DIL02,,human,,WGS,H214_1_3_2,2024-07-07 22:10:37.349467+00:00,"100","30,200",93.0,6
6.67,0.47,38.48
sample1.B_R2_001.fastq.gz,22CGWFLT4,GTACAGCGGA-CAATCTGTGT,6,LIB043956-DIL02,,human,,WGS,H214_1_3_2,2024-07-07 22:10:37.349467+00:00,"100","30,200",89.0,6
6.67,0.52,37.65
```

**BATCH 2:** `Samplemap.tsv`

```plaintext
FASTQ,Flowcell ID,Index Sequence,Flowcell Lane,ESP ID,Pool Name,Species,Illumina Sample Type,Library Type,Library Name,Date Complete,Total Reads,Total Ba
ses,% >Q30,% Pass Filter Clusters,PhiX Error Rate,Avg Q Score
sample3_R1_001.fastq.gz,227NFFLT4,CTCATCATCT-TAAGCCTCTA,4,LIB043957-DIL01,,human,,WGS,H214_3_3_2,2024-06-27 14:22:09.194473+00:00,"100","30,200",76.0,56.
98,0.39,35.02
sample3_R2_001.fastq.gz,227NFFLT4,CTCATCATCT-TAAGCCTCTA,4,LIB043957-DIL01,,human,,WGS,H214_3_3_2,2024-06-27 14:22:09.194473+00:00,"100","30,200",78.0,56.
98,0.48,35.36
sample2.A_R1_001.fastq.gz,227NFFLT4,TTCGTTCACA-CTGATTAGTG,4,LIB043958-DIL01,,human,,WGS,H214_2_6_2,2024-06-27 14:22:09.194473+00:00,"100","30,200",91.0,5
6.98,0.39,38.24
sample2.A_R2_001.fastq.gz,227NFFLT4,TTCGTTCACA-CTGATTAGTG,4,LIB043958-DIL01,,human,,WGS,H214_2_6_2,2024-06-27 14:22:09.194473+00:00,"100","30,200",87.0,5
6.98,0.48,37.42
sample2.B_R1_001.fastq.gz,22CGWFLT4,TTCGTTCACA-CTGATTAGTG,6,LIB043958-DIL02,,human,,WGS,H214_2_6_2,2024-07-07 22:10:37.349467+00:00,"100","30,200",93.0,6
6.67,0.47,38.53
sample2.B_R2_001.fastq.gz,22CGWFLT4,TTCGTTCACA-CTGATTAGTG,6,LIB043958-DIL02,,human,,WGS,H214_2_6_2,2024-07-07 22:10:37.349467+00:00,"100","30,200",89.0,6
6.67,0.52,37.69
```

`rename.tsv`

```plaintext
samplemap_sample_id	revised_sample_id	comments
H214_1_3_2	H214_1	Test samples.
H214_2_6_2	H214_2	Test samples.
H214_3_3_2	H214_3_3_2	Test samples.
```

**OUTPUT (merge_fastq)**

```plaintext
results
|-- H214_1_3_2
|   |-- H214_1_3_2.R1.fastq.gz
|   |-- H214_1_3_2.R1.fastq.gz.MD5
|   |-- H214_1_3_2.R1.fastq.gz.counts
|   |-- H214_1_3_2.R2.fastq.gz
|   |-- H214_1_3_2.R2.fastq.gz.MD5
|   `-- H214_1_3_2.R2.fastq.gz.counts
|-- H214_2_6_2
|   |-- H214_2_6_2.R1.fastq.gz
|   |-- H214_2_6_2.R1.fastq.gz.MD5
|   |-- H214_2_6_2.R1.fastq.gz.counts
|   |-- H214_2_6_2.R2.fastq.gz
|   |-- H214_2_6_2.R2.fastq.gz.MD5
|   `-- H214_2_6_2.R2.fastq.gz.counts
|-- H214_3_3_2
|   |-- H214_3_3_2.R1.fastq.gz
|   |-- H214_3_3_2.R1.fastq.gz.MD5
|   |-- H214_3_3_2.R1.fastq.gz.counts
|   |-- H214_3_3_2.R2.fastq.gz
|   |-- H214_3_3_2.R2.fastq.gz.MD5
|   `-- H214_3_3_2.R2.fastq.gz.counts
|-- __bsub
|   |-- 1_merge_fastq.sh
|   |-- 1_merge_fastq_bsub.err
|   |-- 1_merge_fastq_bsub.out
|   |-- 1_merge_fastq_bsub.sh
|   |-- 1_merge_fastq_bsub.yaml
|   |-- 2_merge_fastq.sh
|   |-- 2_merge_fastq_bsub.err
|   |-- 2_merge_fastq_bsub.out
|   |-- 2_merge_fastq_bsub.sh
|   |-- 2_merge_fastq_bsub.yaml
|   |-- 3_merge_fastq.sh
|   |-- 3_merge_fastq_bsub.err
|   |-- 3_merge_fastq_bsub.out
|   |-- 3_merge_fastq_bsub.sh
|   `-- 3_merge_fastq_bsub.yaml
|-- gtac_read_counts.tsv
|-- gtac_read_counts.tsv.MD5
|-- merged_samplemap.tsv
|-- merged_samplemap.tsv.MD5
|-- merged_samplemap.tsv.pickle
|-- rename.tsv
`-- src_samplemaps
    |-- all_samplemaps.tsv
    |-- all_samplemaps.tsv.MD5
    |-- batch_1
    |   `-- Samplemap.csv
    `-- batch_2
        `-- Samplemap.csv

7 directories, 43 files
```

A breakdown of output files and directories follows.

**FASTQ Merge Directories**

Each sample gets its own FASTQ merge directory.

```plaintext
|-- H214_1_3_2
|   |-- H214_1.R1.fastq.gz
|   |-- H214_1.R1.fastq.gz.MD5
|   |-- H214_1.R1.fastq.gz.counts
|   |-- H214_1.R2.fastq.gz
|   |-- H214_1.R2.fastq.gz.MD5
|   `-- H214_1.R2.fastq.gz.counts
```

The name of the directory is taken from the `Samplemap.csv` file—i.e. the `samplemap_sample_id` field in the `rename.tsv` file. However, the merged FASTQ file name is taken from the `revised_sample_id` column from the `rename.tsv` file. This is useful for a sanity check, as the file path will include both the original and revised sample name for the FASTQ files. Example:

```zsh
H214_1_3_2/H214_1.R2.fastq.gz
H214_1_3_2/H214_1.R1.fastq.gz
```

There will be 1 pair of compressed FASTQ files (R1 and R2) per sample, with corresponding MD5 checksum files. Each FASTQ will also have a `.counts` file which has an independently calculated read count value. The R1 and R2 read counts should be equal.

**The bsub Directory**

The bsub (batch submission) directory contains all of the individual sample merging commands executed across the LSF system at WashU. Except for the LSF log files, the user likely will never need to review these files unless process troubleshooting is required.

Every unique sample gets a set of commands and connected metadata. These commands are automatically executed by the `merge_fastq` command when the `--no-lsf-dry` argument is passed. The `*_merge_fastq.sh` files are bash commands that do the actual file merging. The `*_merge_fastq_bsub.sh` files are responsible for executing the `bsub` commands that run the `*_merge_fastq.sh` commands. The `*_merge_fastq_bsub.yaml` files contain metadata related to the parallel LSF jobs. The `*._merge_fastq_bsub.err` and `*_merge_fastq_bsub.out` are log files generated during LSF job execution.

```plaintext
|-- __bsub
|   |-- 1_merge_fastq.sh
|   |-- 1_merge_fastq_bsub.err
|   |-- 1_merge_fastq_bsub.out
|   |-- 1_merge_fastq_bsub.sh
|   |-- 1_merge_fastq_bsub.yaml
|   |-- 2_merge_fastq.sh
|   |-- 2_merge_fastq_bsub.err
|   |-- 2_merge_fastq_bsub.out
|   |-- 2_merge_fastq_bsub.sh
|   |-- 2_merge_fastq_bsub.yaml
|   |-- 3_merge_fastq.sh
|   |-- 3_merge_fastq_bsub.err
|   |-- 3_merge_fastq_bsub.out
|   |-- 3_merge_fastq_bsub.sh
|   `-- 3_merge_fastq_bsub.yaml
```

**The merged_samplemap.tsv File**

This is the principal dataframe which tracks all of the samples and related FASTQ file information. This file may be used to access the merged FASTQ file paths, sample renaming history, constituent original FASTQ files used for merging, merging commands, etc.

```plaintext
Fields
------
 1. fastq
 2. flow_cell_id
 3. index_sequence
 4. lane_number
 5. read_number
 6. sample_name
 7. library_type
 8. total_bases
 9. samplemap_path
10. gtac_fastq_reads
11. esp_id
12. pool_name
13. batch_id
14. fastq_path
15. project
16. revised_sample_name
17. merged_commands
18. merged_fastq_path
19. gtac_end_pair_reads
20. gtac_sample_reads
```

The dataframe makes asking questions of the merged samples easy. For example, say we want to see the original sample name, the revised sample name, and the final merged FASTQ path for the batch using Python and Pandas.

```python
## Example of interrogating the merged FASTQ dataframe.
import pandas as pd
df = pd.read_pickle('merged_samplemap.tsv.pickle')
df_paths = df[
    ['sample_name', 'revised_sample_name', 'merged_fastq_path']
].drop_duplicates().reset_index(drop=True)
print(df_paths)

##   sample_name revised_sample_name                            merged_fastq_path
## 0  H214_1_3_2              H214_1      /tmp/test/H214_1_3_2/H214_1.R1.fastq.gz
## 1  H214_1_3_2              H214_1      /tmp/test/H214_1_3_2/H214_1.R2.fastq.gz
## 2  H214_3_3_2          H214_3_3_2  /tmp/test/H214_3_3_2/H214_3_3_2.R1.fastq.gz
## 3  H214_3_3_2          H214_3_3_2  /tmp/test/H214_3_3_2/H214_3_3_2.R2.fastq.gz
## 4  H214_2_6_2              H214_2      /tmp/test/H214_2_6_2/H214_2.R1.fastq.gz
## 5  H214_2_6_2              H214_2      /tmp/test/H214_2_6_2/H214_2.R2.fastq.gz
```

**The merged_samplemap.tsv.pickle File**

This is the same information in the `merged_samplemap.tsv` file, but in binary form. The Python data structures used to create the dataframe are retained in the binary version. For example, the `merged_commands` field retains the `list` type for the commands.

**The rename.tsv File**

This is just a copy of the original `rename.tsv` supplied to the `merge_fastq` command, for reference.

**The gtac_read_counts.tsv File**

GTAC@MGI provides some read count per FASTQ file information within their `Samplemap.csv` files. We have collected this information and carried forward read counts for the merged FASTQ. Based on these values, the `gtac_read_counts.tsv` file provides a matrix of samples and their sequence "throughput" based on (1) a set of predefined read count target values and (2) a minimum percent of target read count required to "pass" a sample. The set of predefined targets provides a lookup table for potential targets, as minimum target read counts will vary by project. Read counts are calculated by read pair counts—i.e. R1 + R2 as pairs. Currently, samples are termed "passed" when they meet or exceed 80% of a given target read count. Target read counts being evaluated are as follows: 100000; 200000; 300000; 400000; 500000; 1000000; 1500000; 2000000; 2500000; 3000000; 3500000; 4000000; 4500000; 5000000; 10000000; 20000000; 30000000; 40000000; 50000000.

**The src_samplemaps Directory**

This is just a copy of the original input `Samplemap.csv` files as supplied to the `merge_fastq` commands, for reference. The `all_samplemaps.tsv`  file is a concatenated dataframe of all of the samplemap files prior to any FASTQ merging steps.

```plaintext
-- src_samplemaps
   |-- all_samplemaps.tsv
   |-- all_samplemaps.tsv.MD5
   |-- batch_1
   |   `-- Samplemap.csv
   `-- batch_2
       `-- Samplemap.csv
```

>[!IMPORTANT] Summary
> + The `merge_fastq` command launches parallel jobs to merge FASTQ files and write metadata.
> + You will need to log into WashU's `compute1` to run `merge_fastq`.
> + The `__bsub` directory contains individual sample merging commands executed across the LSF system.
> + The `merged_samplemap.tsv` file is a dataframe which tracks all samples and merging information.
> + The dataframe makes asking questions of the merged samples easy.
> + The `gtac_read_counts.tsv` file provides a matrix of samples and their calculated sequence throughput.

### 4. Review Parallel Processing Logs

Each unique sample gets its own LSF `bsub` job once the `merge_fastq` command has been successfully run. You may check the process state (PEND, RUN, etc.) of your submitted jobs using the `bjobs` command.

When all jobs in the queue have completed processing, `bjobs` will no longer list any `merge_fastq` jobs; however, this does not necessarily mean that all of the jobs have completed successfully. You will still need to manually validate that all of the jobs fully completed before exiting. You can validate that there were no complications by reviewing the LSF logs in the results `__bsub/` directory.

A few simple shell commands can help assess if the jobs properly completed.

First, determine how many total jobs should have run in the LSF queue. This will be equal to the total number of samples in the `rename.tsv` file. The total number of samples should also be equal to the number of shell commands for merging FASTQ run in the `__bsub/` directory, as well as the connected LSF error and output files.

```zsh
## Run commands from the merge_fastq --outdir directory.
ls __bsub/*fastq.sh | wc -l  # total job count
ls __bsub/*err | wc -l  # total job count
ls __bsub/*out | wc -l  # total job count
```

If processing was complete, you may proceed to check for other errors.

We can review the LSF processing logs by searching for potential error keywords This is not an exhaustive search but catches the majority of issues.

```zsh
## Run commands from the merge_fastq --outdir directory.
grep -i "exit" __bsub/*out | wc -l  # should be 0
grep -i "error" __bsub/*err | wc -l  # should be 0
grep -i "success" __bsub/*out | wc -l  # should equal total job count
```

Finally, each sample should produce one pair (R1 and R2) of merged FASTQ files.

```zsh
## Run commands from the merge_fastq --outdir directory.
find * -type f | grep "R1.fastq.gz$" | wc -l  # should equal total job count
find * -type f | grep "R2.fastq.gz$" | wc -l  # should equal total job count
```

If the counts above differ, processing was incomplete and you will need to investigate, fix, and reprocess the FASTQ batch. Provided that you are satisfied all jobs completed correctly, you may now proceed to running the `eval_fastq_counts` command.

>[!IMPORTANT] Summary
> + Each unique sample automatically gets its own LSF bsub job for merging commands.
> + Check the process state of submitted jobs using the `bjobs` command.
> + LSF logs were written to the results `__bsub/` directory.
> + Shell commands can help assess if LSF jobs completed correctly.

### 5. Evaluate the Merged  FASTQ Read Counts

Newly merged FASTQ files were independently evaluated for read counts during the `merge_fastq` processing. Each merged FASTQ file received a corresponding `.counts` file. Using the counts files we may produce a `src_read_counts.tsv` file to review target sequence coverage, much like the previously discussed `gtac_read_counts.tsv` file. Running the `eval_fastq_counts` command will produce `src_read_counts.tsv` and `read_count_comparisons.tsv` files.

The `src_read_counts.tsv` provides a matrix of samples and their sequence "throughput" based on (1) a set of predefined read count target values and (2) a minimum percent of target read count required to "pass" a sample. The set of predefined targets provides a lookup table for potential targets, as minimum target read counts will vary by project. Read counts are calculated by read pair counts—i.e. R1 + R2 as pairs. Currently, samples are termed "passed" when they meet or exceed 80% of a given target read count. Target read counts being evaluated are as follows: 100000; 200000; 300000; 400000; 500000; 1000000; 1500000; 2000000; 2500000; 3000000; 3500000; 4000000; 4500000; 5000000; 10000000; 20000000; 30000000; 40000000; 50000000.

The `read_count_comparisons.tsv` compares the values in `gtac_read_counts.tsv` (GTAC@MGI reported counts) to `src_read_counts.tsv` (`merged_fastq` counts) and reports any differences. This file highlights any differences in sample read counts between the sequencing core and the final, merged FASTQ files.

>[!IMPORTANT] Summary
> + FASTQ files were independently evaluated for read counts during merging.
> + Each merged FASTQ file receives a `.counts` file.
> + The `src_read_counts.tsv` file provides a matrix of merged samples and their calculated sequence throughput.
> + The `read_count_comparisons.tsv` compares the counts from GTAC@MGI to `merged_fastq` counts and reports any differences.

## Commands

The following commands are supported in the latest version of the `merge_fastq` Docker image. See [https://hub.docker.com/r/twylie/merge_fastq/tags](https://hub.docker.com/r/twylie/merge_fastq/tags) for details.

To pull the latest image and start an interactive shell, do the following.

```zsh
## Example is for Mac running Docker Desktop.
docker pull twylie/merge_fastq:latest
docker container run --platform linux/amd64 -it twylie/merge_fastq:latest zsh
merge_fastq

## usage: merge_fastq [-h] [--version] [--lsf-image STR] [--lsf-group STR] [--lsf-queue STR]
##         [--lsf-dry] [--no-lsf-dry] --samplemap FILE [FILE ...] --outdir DIR --rename FILE
##         --lsf-vol PATH [PATH ...] --project {MIDAS,PLACENTA,PTLD}
## merge_fastq: error: the following arguments are required: --samplemap, --outdir, --rename,
##                     --lsf-vol, --project
```

### 1. Running `prep_rename_file` Command

#### Usage

```plaintext
usage: prep_rename_file [-h] [--version] --samplemap FILE [FILE ...]
       --rename-out FILE

Prepare a rename file for merging FASTQ files with merge_fastq.

options:
  -h, --help            Display the extended usage statement.
  --version             Display the software version number.

required:
  --samplemap FILE [FILE ...]
                        List of samplemap files to be evaluated.
  --rename-out FILE     File path to write the rename file.
```

#### Example Commands

```zsh
prep_rename_file \
    --samplemap "batch1/Samplemap.csv" "batch2/Samplemap.csv" \
    --rename-out "prep-rename.tsv"
```

### 2. Running `merge_fastq` Command

#### Usage

```plaintext
usage: merge_fastq [-h] [--version] [--lsf-image STR] [--lsf-group STR] [--lsf-queue STR]
       [--lsf-dry] [--no-lsf-dry] --samplemap FILE [FILE ...] --outdir DIR --rename FILE
       --lsf-vol PATH [PATH ...] --project
                   {MIDAS,PLACENTA,PTLD}

Merge FASTQ files at sample-level.

options:
  -h, --help            Display the extended usage statement.
  --version             Display the software version number.
  --lsf-image STR       Docker image for LSF processing. [twylie/merge_fastq:latest]
  --lsf-group STR       Group for LSF processing. [compute-kwylie]
  --lsf-queue STR       Queue for LSF processing. [general]
  --lsf-dry             Dry run for LSF processing. (default)
  --no-lsf-dry          Executes LSF processing.

required:
  --samplemap FILE [FILE ...]
                        List of samplemap files to be merged.
  --outdir DIR          Output directory path to write results.
  --rename FILE         Path to a sample rename file.
  --lsf-vol PATH [PATH ...]
                        Top-level Docker volume path(s) for WashU LSF processing.
  --project {MIDAS,PLACENTA,PTLD}
                        Project tag/name.
```

#### Example Commands

```zsh
merge_fastq \
    --samplemap "batch1/Samplemap.csv" "batch2/Samplemap.csv" \
    --outdir "/tmp/test" \
    --rename "rename.tsv" \
    --lsf-vol "/storage1/fs1/PTB" "/home/twylie" \
    --project "PTLD"
```

### 3. Running `eval_fastq_counts` Command

#### Usage

```plaintext
usage: eval_fastq_counts [-h] [--version] --merged-samplemap FILE
       --gtac-counts FILE --outdir DIR

Evaluate and compare read counts for merged FASTQ files.

options:
  -h, --help            Display the extended usage statement.
  --version             Display the software version number.

required:
  --merged-samplemap FILE
                        A merged_samplemap.tsv file.
  --gtac-counts FILE    A gtac_read_counts.tsv file.
  --outdir DIR          Output directory path to write results.
```

#### Example Commands

```zsh
eval_fastq_counts \
    --merged-samplemap "/tmp/test/merged_samplemap.tsv" \
    --gtac-counts "/tmp/test/gtac_read_counts.tsv" \
    --outdir "/tmp/test/"
```

## Test Data Set and Demonstration

###  Data

A set of contrived test data are included with this distribution for demonstration purposes. You may download the test data by cloning this GitHub repository and locating the `test_data/` directory. Use these data for the demo outlined below.

```plaintext
test_data
├── batch1
│   ├── Samplemap.csv
│   ├── sample1.A_R1_001.fastq.gz
│   ├── sample1.A_R2_001.fastq.gz
│   ├── sample1.B_R1_001.fastq.gz
│   └── sample1.B_R2_001.fastq.gz
├── batch2
│   ├── Samplemap.csv
│   ├── sample2.A_R1_001.fastq
│   ├── sample2.A_R2_001.fastq.gz
│   ├── sample2.B_R1_001.fastq.gz
│   ├── sample2.B_R2_001.fastq
│   ├── sample3_R1_001.fastq
│   └── sample3_R2_001.fastq
├── cmds.sh
└── rename.tsv

3 directories, 14 files
```

### Running Test Data Demo

The following demo assumes you are an employee in the Wylie Lab at WashU with proper credentials and access to `compute1` and `storage1`. Do not use `twylie` as the user when running these commands, but rather use your own user name, home directory, space, etc. The demo is using `scratch1` space; however, use larger `storage1` space when running real FASTQ batches.

==TK== Add instructions for cloning the GitHub repository here.

```zsh
## Log into WashU.

ssh compute1-client-1.ris.wustl.edu

## Change directory to your scratch area and create directories for
## processing.

cd /scratch1/fs1/twylie/
mkdir merge_fastq_demo
cd merge_fastq_demo/

## Clone the merge_fastq GitHub repository and change directory to the
## test data set.

## Clone the repository step here.

cd test_data/

## Get an interactive Docker session by calling the merge_fastq:latest
## Docker image. You should still be in the test_data/ directory.

LSF_DOCKER_VOLUMES='/scratch1/fs1/twylie/merge_fastq_demo/test_data:/scratch1/fs1/twylie/merge_fastq_demo/test_data' \
                  bsub -Is -R 'select[mem>16GB] rusage[mem=16GB]' \
                  -G compute-kwylie \
                  -q general-interactive \
                  -a 'docker(twylie/merge_fastq:latest)' \
                  zsh

## Running the prep_rename_file command is optional since the test_data/
## directory already has a rename.tsv file. However, running the
## prep_rename_file command looks like the following.

prep_rename_file \
    --samplemap "batch1/Samplemap.csv" "batch2/Samplemap.csv" \
    --rename-out "-prep_rename.tsv"

## Now run the merge_fastq command using the test_data/rename.tsv file as
## input. You should include your home directory and the working
## directory in the --lsf-vol argument. We are passiong two Samplemap.csv
## because there are two sequencing batches in the test_data/ set.

merge_fastq \
    --samplemap "batch1/Samplemap.csv" "batch2/Samplemap.csv" \
    --outdir "results" \
    --rename "rename.tsv" \
    --lsf-vol "/home/twylie" "/scratch1/fs1/twylie" \
    --project "MIDAS" \
    --no-lsf-dry

## Monitor the LSF job progress using the bsub command.

bjobs

## When all of the LSF jobs have completed, you can validate that there
## were no complications by reviewing the LSF logs in the __bsub/
## directory. Also, half the number of merged FASTQ files (R1 + R2 are a
## sample pair) should equal the total number samples in the rename.tsv
## file.

grep -i "exit" __bsub/*out | wc -l  # should be 0
grep -i "error" __bsub/*err | wc -l  # should be 0
ls __bsub/*fastq.sh | wc -l  # total job count
grep -i "success" __bsub/*out | wc -l  # should equal total job count
find * -type f | grep "R1.fastq.gz$" | wc -l  # should equal total job count
find * -type f | grep "R2.fastq.gz$" | wc -l  # should equal total job count

## Provided that you are satisfied all jobs completed correctly, you may
## now run the eval_fastq_counts command.

eval_fastq_counts \
    --merged-samplemap "results/merged_samplemap.tsv" \
    --gtac-counts "results/gtac_read_counts.tsv" \
    --outdir .
```

## Contact

**Todd N. Wylie, Associate Professor** \
Department of Pediatrics \
Division of Gastroenterology, Hepatology, & Nutrition \
Washington University School of Medicine \
660 S. Euclid Avenue, MSC 8208-0043-10 \
St. Louis, MO 63110 \
email — twylie@wustl.edu

## License

*Begin license text.*

**MIT License**
https://opensource.org/license/mit

Copyright 2024 T.N. Wylie. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*End license text.*

<!-- End document text. -->
