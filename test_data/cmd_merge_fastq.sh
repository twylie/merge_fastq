#!/usr/bin/sh
merge_fastq \
    --samplemap batch1/Samplemap.csv batch2/Samplemap.csv \
    --outdir /tmp/test \
    --rename rename.tsv \
    --lsf-vol "/foo" "/bar" \
    --project PTLD

sh /tmp/test/__bsub/1_merge_fastq.sh
sh /tmp/test/__bsub/2_merge_fastq.sh
sh /tmp/test/__bsub/3_merge_fastq.sh

tree /tmp/test
