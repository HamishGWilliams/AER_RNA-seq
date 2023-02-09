#!/bin/bash

#SBATCH --job-name=FeatureCounts
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk

echo "running FeatureCounts will all .bam files"

module load subread/2.0.2

featureCounts -p -g 'ID' -t exon -O -a genes/A_Equina.gff3 -o A_Equina_featureCounts_output.txt *.bam

echo "Counts read!"

### Options used:
# -p: If specified, libraries are assumed to contain paired-end
                     # reads. For any library that contains paired-end reads, the
                      #'countReadPairs' parameter controls if read pairs or reads
                      # should be counted.

# -g:  Specify attribute type in GTF annotation. 'gene_id' by
                      ##default. Meta-features used for read counting will be
                      # extracted from annotation using the provided value.

			## 'Parent' should be specified.


# -t: Specify feature type(s) in a GTF annotation. If multiple
                     # types are provided, they should be separated by ',' with
                     # no space in between. 'exon' by default. Rows in the
                     # annotation with a matched feature will be extracted and
                     # used for read mapping.

# -O: Assign reads to all their overlapping meta-features (or
                     # features if -f is specified).