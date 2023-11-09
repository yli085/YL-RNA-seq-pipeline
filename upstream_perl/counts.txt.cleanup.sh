## usage
## A bash script to clean up the featureCount output table.
## Delete the long path name/bam suffix in the sample name, and delete the 2-6 columns
## change: change the path name (within the count table sample name) for the target dataset
## run: sh counts.txt.cleanup.sh

cp counts.txt counts2.txt

## change the project id in the path
less counts2.txt | sed '1d'| sed 's/\/home\/yli\/work\/RNAseq\/project_id\/data\/STAR\///g' | sed 's/\.Aligned\.out\.bam//g' > counts_renamed.txt

cut --complement -f2-6 counts_renamed.txt > counts_Rmatrix.txt

rm counts2.txt
