#!/bin/bash

#SBATCH -c 3                # number of core to be used
#SBATCH -t 2-00:00          # estimated run-time in D-HH:MM
#SBATCH -p mid              # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=100000        # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# variant calling via gatk: http://gatkforums.broadinstitute.org/gatk/discussion/6800/known-sites-for-indel-realignment-and-bqsr-in-hg38-bundle; https://www.broadinstitute.org/gatk/guide/article?id=1247; ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle
# samtools conversion, mpileup, snp calling, change flow cell id
# haplotypecaller troubles with downsampling and depth calculation - discussion on: http://gatkforums.broadinstitute.org/discussion/3356/haplotypecaller-dp-reports-low-values
# https://software.broadinstitute.org/gatk/documentation/article?id=3891


# Get sample name
sample=${PWD##*/}
# variables
gatk=/path/to/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar
fasta=/path/to/Genomes/gatk/hg38/Homo_sapiens_assembly38.fasta
picard=/path/to/picard-2.9.2/picard.jar


cd ../../bam


date
echo "Get same chromosome names as in gatk hg38 fasta file chr1-22,chrXYM..."
samtools view -@ 3 -b -o $sample.f.bam $sample.rmdup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
samtools index $sample.f.bam
echo "Done."
echo ""


date
echo "Assign read groups..."
flowcell=xxx
java -XX:ParallelGCThreads=2 -jar $picard AddOrReplaceReadGroups I=$sample.f.bam O=$sample.rmdup.groups.bam RGID=$sample RGLB=gatc_kit RGPL=ILLUMINA RGSM=$sample RGPU=$flowcell
samtools index $sample.rmdup.groups.bam
rm $sample.f.ba*
echo "Done."
echo ""


date
echo "RNA-Seq: Split'N'Trim and reassign mapping qualities..."
java -jar -Xmx6g -XX:ParallelGCThreads=2 $gatk -T SplitNCigarReads -R $fasta -I $sample.rmdup.groups.bam -o $sample.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
rm $sample.rmdup.groups.ba*
echo "Done."
echo ""


date
echo "Variant calling for RNA-Seq data..."
java -jar -Xmx6g -XX:ParallelGCThreads=2 $gatk -T HaplotypeCaller -R $fasta -I $sample.split.bam -dontUseSoftClippedBases -o ../snp/$sample.vcf -stand_call_conf 20.0
echo "Done."
echo ""



rm $sample.split.ba*



date
