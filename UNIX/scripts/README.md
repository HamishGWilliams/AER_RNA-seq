# Overview

## Description
This folder will contain all of the relevant UNIX scripts which are generated for UNIX for all purposes. Some of these purposes include:
- Loading and building commands to be used
- Shell scripts of components of analyses
- A complete master script which will run all of the individual components consequtively

## Order of Components
The script components should be ran in the following order:
1. **alignment**
  - *Align* the sample reads (1 + 2) to a **reference genome**, generating a .SAM file
  - NOTE: HISAT2 works with .zip files, so there is no need to unpack the read files
2. **sam2bam**
  - *Sorts* and *Converts* .sam files to .bam files
3. **index**
  - Indexes BAM files using SAMtools
4. **assembleANDquantify**
  - *Assembles** the reads and *quantifies* them using **StringTie**
5. **merge**
  - *Merges* the assembled  

Will update the above section when I have completed a successful run of the example data
----------------------------------------------------------------------------------------------

## How to use scripts idea:
I'm not sure how well this will work, but I could try to pull the scripts at all times from the GitHub so that I can either edit on a seperate interface, or just within GitHub to test UNIX code. 
