## Usage
## A bash script to run featureCounts on a set of BAM files
## Use protein coding only gtf as the reference for gene counting
## Run: sh runFeatureCounts.sh
## change: need to change the project ID in the file path.

featureCounts -p --countReadPairs -T 8 -s 2 -M --primary \
-a ~/work/programs/genomes/GRCm39.vM27/mouse.gencode.vM27.primary_assembly.annotation.protein_coding_only.gtf \
-o counts.txt \
~/work/RNAseq/project_id/data/STAR/*.bam

