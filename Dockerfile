# Project     : merge_fastq
# File Name   : Dockerfile
# Description : Create a Docker image for merging FASTQ.
# Author      : Todd N. Wylie
# Email       : twylie@wustl.edu
# Created     : Mon Sep 30 22:40:51 CDT 2024
# Copyright   : Copyright (C) 2024 by T.N. Wylie. All rights reserved.

FROM ubuntu:22.04
MAINTAINER T.N. Wylie <twylie@wustl.edu>
LABEL description = "merge_fastq: Merge split GTAC@MGI FASTQ files."

# Updates and Upgrades ########################################################

RUN apt-get update -y && apt-get upgrade -y 

RUN apt-get install -y \
    build-essential

RUN apt-get install -y \
    autoconf \
    automake

RUN apt-get install -y \
    zsh

# Essential Tools #############################################################

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    bzip2 \
    cmake \
    git \
    unzip \
    wget \
    curl \
    libbz2-dev \
    zile \
    tree \
    bc

# Python ######################################################################

# Python 3.10

RUN apt-get install -y python3-pip
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade setuptools
RUN python3 -m pip install pandas
RUN pip3 install typing-extensions
RUN pip3 install biopython
RUN pip3 install pyyaml

# Clean-Up ####################################################################

RUN apt-get clean

# Local Python Installations ##################################################

COPY ./bin/prep_rename_file.py /usr/bin/prep_rename_file
COPY ./bin/prep_samplemap_file.py /usr/bin/prep_samplemap_file
COPY ./bin/merge_fastq.py /usr/bin/merge_fastq
COPY ./bin/eval_fastq_counts.py /usr/bin/eval_fastq_counts
COPY ./mergefastq /usr/lib/python3.10/mergefastq
COPY ./test_data /test_data

# __END__
