#!/bin/bash/

# Load Packages needed
module load python/3.9.12/


# Constructing Reference Index

## Extract splice sites:
extract_splice_sites.py ./genes/A_Equina.gff3 >chrX.ss

## Extract Exons:
extract_exons.py genes/A_Equina.gff3 >A_Equina.exon

## Assemble Index:
$ hisat2-build --ss A_Equina.ss --exon A_Equina.exon RefGenome/A_Equina_RefGenome.fa A_Equina_tran

## Move into specific index folder:
mkdir ./indexes
mv *_tran* ./indexes/
