#!/bin/bash
# Differential expression (96 h vs 40 h), peroxidase validation and figures.
set -euo pipefail
mkdir -p results/deseq2
Rscript scripts/deseq2_analysis.R
Rscript scripts/check_peroxidases.R
Rscript scripts/plot_peroxidases.R
