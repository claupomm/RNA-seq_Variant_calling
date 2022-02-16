#!/bin/bash

#SBATCH -c 2                # number of core to be used
#SBATCH -t 1-12:00          # estimated run-time in D-HH:MM
#SBATCH -p mid              # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=10000         # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}
# variables
fasta=/path/to/Genomes/gatk/hg38/Homo_sapiens_assembly38.fasta
varscan=/path/to/VarScan.v2.4.3.jar


cd ../../snp


date
echo "Examine mutations via VarScan, single nucleotide mutations..."
samtools mpileup -q 10 -Q 20 -f $fasta ../bam/$sample.rmdup.bam | java -jar $varscan mpileup2snp -min-coverage 10 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.20 --p-value 0.01 --output-vcf 1 > $sample.varscan.vcf
echo "Done."
echo ""

date
echo "Examine mutations via VarScan, indels..."
samtools mpileup -q 10 -Q 20 -f $fasta ../bam/$sample.rmdup.bam | java -jar $varscan mpileup2indel -min-coverage 10 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.20 --p-value 0.01 --output-vcf 1 > $sample.varscan.indel.vcf
echo "Done."
echo ""

date
echo "Prepare varscan results for vep, merge indels and snps..."
grep -v "^#" $sample.varscan.indel.vcf > $sample.tmp
cat $sample.varscan.vcf $sample.tmp > $sample.tmp2
mv $sample.tmp2 $sample.varscan.vcf
sed 's/^chrM/chrMT/' $sample.varscan.vcf > $sample.tmp
sed 's/^chr//' $sample.tmp > $sample.ohneChr.varscan.vcf
rm $sample.tmp $sample.tmp2 $sample.varscan.indel.vcf
echo "Done."
echo ""

date
echo "Filter chromosomes 1-22, x,y,mt only..."
grep "^#" $sample.ohneChr.varscan.vcf > $sample.tmp
grep "^[1-9MYX]" $sample.ohneChr.varscan.vcf > $sample.tmp2
cat $sample.tmp $sample.tmp2 > $sample.ohneChr.varscan.vcf
rm $sample.tmp $sample.tmp2
echo "Done."
echo ""




date
