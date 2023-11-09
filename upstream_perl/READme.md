These are the perl/bash codes for the RNA-seq upstream analysis:
1. combine_fastq.pl: combine reads from 4 seperate lanes for illumina sequencing reads.
2. fastqc.pl: conduct reads QC using fastqc program.
3. runSTAR.pl: run STAR for reads mapping.
4. runFeatureCounts.sh: run featureCounts for gene expression quantification. This is a bash script.
5. counts.txt.cleanup.sh: a bash script to clean up the count matrix, i.e., shorten the sample name and delete the first line.
6. submit_job.sh: Bash script to run the perl in linux.
