#!/bin/bash
# Usage: sbatch --array=1-N mutation_calling_from_bam.sh <CFDNA_SAMPLE_FILE> <NORMAL_SAMPLE_FILE> <CFDNA_DIR> <NORMAL_DIR> [OUTPUT_DIR] [BED_DIR]
sbatch \
    --job-name=B24_1 \
    --time=1-00:00:00 \
    --ntasks=1 \
    --cpus-per-task=6 \
    --mem-per-cpu=6G \
    --mail-type=ALL \
    --mail-user=ssu42 \
    --partition=emoding \
    --array=1-5 \
    /oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/mutation_calling_from_bam.sh \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB24/NovaSeqB24_1/sample2barcodeB24_1.txt" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB25/sample2barcodeB25.txt" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB24/NovaSeqB24_1/demultiplexed/barcode-deduped/cfdna" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB25/demultiplexed/barcode-deduped/normal" \
    "/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB24/NovaSeqB24_1/results_somatic_calling"