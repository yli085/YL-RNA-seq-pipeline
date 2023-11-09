## usage
## A perl script to repetitively run STAR on fastq files
## STAR parameters: --twopassMode Basic
## --clip5pNbases 1 1
## change: fastq file directory
## run: perl runSTAR.pl fastq_names.txt

use warnings;
use strict;

## input fastq file directory
my $file_dir = "/home/yli/work/RNAseq/project_id/data/combined_fastq/";

# load the genome and exit
system "STAR --genomeDir /home/yli/work/programs/genomes/GRCm39.vM27/GRCm39_genome_index/ --genomeLoad LoadAndExit";

# loop
while (<>) {
    chomp;
    print "Processing alignment for $_ file...\n";

    my $file_prefix = $file_dir.$_;
    system "STAR --runThreadN 8 \\
--genomeDir /home/yli/work/programs/genomes/GRCm39.vM27/GRCm39_genome_index \\
--twopassMode Basic \\
--readFilesCommand zcat \\
--outSAMunmapped Within \\
--outSAMtype BAM Unsorted \\
--readFilesPrefix $file_prefix \\
--readFilesIn _R1.fastq.gz _R2.fastq.gz \\
--clip5pNbases 1 1 \\
--outFileNamePrefix $_.";

}

# remove
system "STAR --genomeDir /home/yli/work/programs/genomes/GRCm39.vM27/GRCm39_genome_index/ --genomeLoad Remove";

## end


