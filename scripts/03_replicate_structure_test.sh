#!/bin/bash
# Quantify the first 5 M reads of each of the 18 runs separately and plot a run-level PCA per time point.
# Needs the SRA caches in data/raw (step 02) and the Salmon index (step 01).
set -euo pipefail
mkdir -p results/quant_test results/deseq2
while IFS=$'\t' read -r NAME SRR; do
  [ -z "$NAME" ] && continue
  (cd data/raw && fastq-dump -X 5000000 -Z "$SRR" </dev/null) > tmp_test.fastq
  salmon quant -i data/reference/salmon_index -l A -r tmp_test.fastq \
    -p 8 --validateMappings -o "results/quant_test/$NAME" </dev/null
done < annotation/run_accessions.tsv
rm -f tmp_test.fastq
Rscript scripts/replicate_structure_test.R
