#!/bin/bash/

### Alignment of RNA reads to Reference genome

## map each sample using 1 threads

# 1. Generate new directory for alignments
mkdir map

# 2. perform alignment
# NOTE: you may need to change the paths to the data depending on how you have set
# up your directory!!
hisat2 -p 1 --rg-id=ERR188044 --rg SM:ERR188044  --rg PL:ILLUMINA --rg PU:aaaaaaa -x chrX_data/indexes/chrX_tran --dta --rna-strandness RF -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz -2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S map/ERR188044_chrX.sam
hisat2 -p 1 --rg-id=ERR188104 --rg SM:ERR188104  --rg PL:ILLUMINA --rg PU:aaaaaaa -x chrX_data/indexes/chrX_tran --dta --rna-strandness RF -1 chrX_data/samples/ERR188104_chrX_1.fastq.gz -2 chrX_data/samples/ERR188104_chrX_2.fastq.gz -S map/ERR188104_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188234_chrX_1.fastq.gz -2 chrX_data/samples/ERR188234_chrX_2.fastq.gz -S map/ERR188234_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran   -1 chrX_data/samples/ERR188245_chrX_1.fastq.gz -2 chrX_data/samples/ERR188245_chrX_2.fastq.gz -S map/ERR188245_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188257_chrX_1.fastq.gz -2 chrX_data/samples/ERR188257_chrX_2.fastq.gz -S map/ERR188257_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran-1 chrX_data/samples/ERR188273_chrX_1.fastq.gz -2 chrX_data/samples/ERR188273_chrX_2.fastq.gz -S map/ERR188273_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188337_chrX_1.fastq.gz -2 chrX_data/samples/ERR188337_chrX_2.fastq.gz -S map/ERR188337_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188383_chrX_1.fastq.gz -2 chrX_data/samples/ERR188383_chrX_2.fastq.gz -S map/ERR188383_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188401_chrX_1.fastq.gz -2 chrX_data/samples/ERR188401_chrX_2.fastq.gz -S map/ERR188401_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz -2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S map/ERR188428_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188454_chrX_1.fastq.gz -2 chrX_data/samples/ERR188454_chrX_2.fastq.gz -S map/ERR188454_chrX.sam
hisat2 -p 1 -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR204916_chrX_1.fastq.gz -2 chrX_data/samples/ERR204916_chrX_2.fastq.gz -S map/ERR204916_chrX.sam

echo "If an error has occured with reference to 'no file with name, etc...', check the pathways."

