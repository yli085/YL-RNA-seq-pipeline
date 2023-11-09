## usage:
## A perl script to repetitively combine fastq reads from different lanes of each one of samples into one fastq file ## Here is for PE reads, do R1 and R2 separately.
## change: change the raw fastq file directory
## run: combine_fastq.pl fastq_names.txt

use strict;
use warnings;

## change the fastq file directory here, add / at the end
my $filedir = "/path/to/Fastq/";

while (<>) {
	chomp;
	my $current = $_;
	system "cat ${filedir}${current}_L00*_R1_001.fastq.gz > ${current}_R1.fastq.gz";
	system "cat ${filedir}${current}_L00*_R2_001.fastq.gz > ${current}_R2.fastq.gz";
}


