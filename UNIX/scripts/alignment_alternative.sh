#!/bin/bash/

### Alternative method for aliging the reads to reference genome using HISAT2

## bash script for hisat2; align all .fastq.gz files to indexed reference genome to generate .sam files

# Name each sample specifically
SAMPLES="ERR188044 ERR188104 ERR188234 ERR188245 ERR188257 ERR188273 ERR188337 ERR188383 ERR188401 ERR188428 ERR188454 ERR204916"

# Loop using the 'SAMPLES' vector to align all samples to refernce genome
for SAMPLE in $SAMPLES; do
    hisat2 -p 11 --dta -x ~/chrX_data/indexes/chrX_tran -1 ~/chrX_data/samples/${SAMPLE}_chrX_1.fastq.gz -2 ~/chrX_data/samples/${SAMPLE}_chrX_2.fastq.gz -S ${SAMPLE}_chrX.sam
done

## Finally, move all of the files into a single folder:

# Again, we can use a loop for this:
for SAMPLE in $SAMPLES; do
	mv ${SAMPLE}_chrX.sam
done

