#!/bin/bash
# Usage: sbatch --array=1-N mutation_calling_from_bam.sh <CFDNA_SAMPLE_FILE> <NORMAL_SAMPLE_FILE> <CFDNA_DIR> <NORMAL_DIR> [OUTPUT_DIR] [BED_DIR]
sbatch \
    --job-name=B23_2 \
    --time=2-00:00:00 \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem-per-cpu=10G \
    --mail-type=ALL \
    --mail-user=ssu42 \
    --partition=emoding \
    --array=1-6 \
    /oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/mutation_calling_from_bam.sh \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB23/NovaSeqB23_2/sample2barcodeB23_2.txt" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB25/sample2barcodeB25.txt" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB23/NovaSeqB23_2/demultiplexed/barcode-deduped/cfdna" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB25/demultiplexed/barcode-deduped/normal" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB23/NovaSeqB23_2/results_somatic_calling"