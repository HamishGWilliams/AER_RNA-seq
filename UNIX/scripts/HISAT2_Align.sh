#!/bin/bash/

#SBATCH --job-name=alignment
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk

SAMPLES="1-1 1-19 2-20 3-21 4-2 4-22 5-3 6-4 6-23 7-5 7-24 8-6 9-7 9-25 10-8 10-26 11-9 11-27 12-10 12-28 14-11 15-12 16-13 18-14 19-15 20-16 21-17 24-18"

for SAMPLE in $SAMPLES; do
        hisat2 -p 8 -x indexes/A_Equina_tran -1 samples/${SAMPLE}_1.fastq.gz -2 samples/${SAMPLE}_2.fastq.gz -S ${SAMPLE}.sam 
	samtools sort -@ 8 -o ${SAMPLE}.bam ${SAMPLE}.sam
	rm ${SAMPLE}.sam
done