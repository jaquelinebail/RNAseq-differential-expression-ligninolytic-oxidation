#!/bin/bash
# MultiQC report and Salmon quantification: the 3 technical runs of each barcode are one biological sample.
set -euo pipefail
multiqc results/qc/ -o results/qc/multiqc/
mkdir -p results/quant
for T in 40h 96h; do for B in bc1 bc2 bc3; do
  salmon quant -i data/reference/salmon_index -l A \
    -r data/trimmed/${T}_rep1_${B}.trim.fastq.gz data/trimmed/${T}_rep2_${B}.trim.fastq.gz data/trimmed/${T}_rep3_${B}.trim.fastq.gz \
    -p 8 --validateMappings -o results/quant/${T}_${B}
done; done
for S in 40h_bc1 40h_bc2 40h_bc3 96h_bc1 96h_bc2 96h_bc3; do
  echo -e "$S\t$(grep -i 'mapping rate' results/quant/$S/logs/salmon_quant.log | sed 's/.*: //')"
done > results/quant/mapping_rates.tsv
