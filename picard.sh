#!/bin/bash

#SBATCH -c 3                # number of core to be used
#SBATCH -t 0-06:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=50000         # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}

picard=/path/to/picard-2.9.2/picard.jar
genome=/path/to/Genomes/gencode/26/gencode.v26.refFlat.txt
fasta=/path/to/Genomes/gencode/26/GRCh38.p10.genome.fa
cd ../../bam


date
echo "Remove duplicates via picard..."
java -XX:ParallelGCThreads=2 -Xmx20g -jar $picard MarkDuplicates I=$sample.s.bam O=$sample.rmdup.bam METRICS_FILE=$sample.metrics ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE
echo "Done."
date

rm $sample.s.bam
samtools index $sample.rmdup.bam

date
echo "Calculate rna metrix via picard..."
java -XX:ParallelGCThreads=2 -Xmx20g -jar $picard CollectRnaSeqMetrics I=$sample.rmdup.bam O=$sample.rmdup.RNA_Metrics REF_FLAT=$genome STRAND=SECOND_READ_TRANSCRIPTION_STRAND
java -XX:ParallelGCThreads=2 -Xmx20g -jar $picard CollectAlignmentSummaryMetrics R=$fasta I=$sample.rmdup.bam O=$sample.output.txt
echo "Done."
date
