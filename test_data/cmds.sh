mkdir /tmp/renaming

prep_rename_file \
    --samplemap batch1/Samplemap.csv batch2/Samplemap.csv \
    --rename-out /tmp/renaming/prep-rename.tsv

merge_fastq \
    --samplemap batch1/Samplemap.csv batch2/Samplemap.csv \
    --outdir /tmp/test \
    --rename rename.tsv \
    --lsf-vol "/home/twylie" "/scratch1/fs1/twylie" \
    --project PTLD

sh /tmp/test/__bsub/1_merge_fastq.sh
sh /tmp/test/__bsub/2_merge_fastq.sh
sh /tmp/test/__bsub/3_merge_fastq.sh

tree /tmp/test

eval_fastq_counts \
    --merged-samplemap /tmp/test/merged_samplemap.tsv \
    --gtac-counts /tmp/test/gtac_read_counts.tsv \
    --outdir /tmp/test/
