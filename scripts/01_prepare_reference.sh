#!/bin/bash
# Reference genome, annotation and Salmon index. Run from the project root.
set -euo pipefail
mkdir -p data/reference
cd data/reference
BASE=https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-63
wget -O genome.fasta.gz "$BASE/fasta/phanerochaete_chrysosporium/dna/Phanerochaete_chrysosporium.GCA000167175v1.dna.toplevel.fa.gz"
wget -O annotation.gff.gz "$BASE/gff3/phanerochaete_chrysosporium/Phanerochaete_chrysosporium.GCA000167175v1.63.gff3.gz"
gunzip -f genome.fasta.gz annotation.gff.gz
gffread -w transcripts.fa -g genome.fasta annotation.gff
salmon index -t transcripts.fa -i salmon_index -k 31
