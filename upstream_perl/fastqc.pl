## usage:
## A perl script to repetitively run fastqc on fastq files
## run for both of the ends of the PE reads
## run: perl fastqc.pl fastq_names.txt 

use warnings;
use strict;

# input fastq file directory here
my $filedir = "/home/yli/work/RNAseq/project_ID/data/combined_fastq";

while (<>) {
        chomp;
        my $fastq_R1 = $filedir."/$_"."_R1.fastq.gz";
        system "fastqc -t 32 --outdir=./ $fastq_R1";
        my $name_zip1 = $_."_R1_fastqc.zip";
        system "unzip -q $name_zip1";
        system "rm $name_zip1";	

	my $fastq_R2 = $filedir."/$_"."_R2.fastq.gz";
        system "fastqc -t 32 --outdir=./ $fastq_R2";
        my $name_zip2 = $_."_R2_fastqc.zip";
        system "unzip -q $name_zip2";
        system "rm $name_zip2";
}
