# Overview

## Description
This folder will contain all of the relevant UNIX scripts which are generated for UNIX for all purposes. Some of these purposes include:
- Loading and building commands to be used
- Shell scripts of components of analyses
- A complete master script which will run all of the individual components consecutively

## Order of Components
The script components should be ran in the following order:
1. **Index_Assembly**
  - Generate an index of the reference genome, marking the *Splice Sites* and *Exons*.
2. **HISAT2_align**
  - Align the sample data to a reference genome (.gff3 file). Generates .sam files
  - Converts .sam files -> .bam files
  - sorts .bam files
  - removes redundant .sam files
3. **FeatureCounts**
  - Counts the read data in reference to a paired-end reference library (.gff3 file)
