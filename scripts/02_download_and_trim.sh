#!/bin/bash
# Download each run, convert to FASTQ, trim with fastp, delete the raw FASTQ. Run from the project root.
set -euo pipefail
mkdir -p data/raw data/trimmed results/qc tmp_fq
while IFS=$'\t' read -r NAME SRR; do
  [ -z "$NAME" ] && continue
  if [ -s "data/trimmed/$NAME.trim.fastq.gz" ]; then echo "skip $NAME"; continue; fi
  (cd data/raw && prefetch "$SRR" </dev/null)
  (cd data/raw && fasterq-dump "$SRR" -O ../../tmp_fq/ -e 8 </dev/null)
  fastp -w 8 -i "tmp_fq/$SRR.fastq" -o "data/trimmed/$NAME.trim.fastq.gz" \
    --json "results/qc/${NAME}_fastp.json" --html "results/qc/${NAME}_fastp.html" </dev/null
  rm "tmp_fq/$SRR.fastq"
done < annotation/run_accessions.tsv
