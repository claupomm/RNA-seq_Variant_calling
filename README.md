# RNA-seq_Variant-calling
Mutations can also be seen in RNA-seq data although correct allele frequency might not be captured by this. Nevertheless, this workflow describes how to extract mutations from RNA-seq data, which might be useful, if no genome sequencing data is not available. It includes download of public data (GEO/SRA), trimming, mapping, removing duplicate reads, mutation calling, and variant effect prediction/annotation.


## Prepare the project folder and fastq files
Set the Project folder:
```
DIR=/path/to/Project_folder
mkdir $DIR
cd $DIR
```
Download public data:
```
# geo: GSE181130
# Instrument: Illumina HiSeq 2500
# Strategy: RNA-Seq
# Source: TRANSCRIPTOMIC
# Selection: cDNA
# Layout: SINGLE
# Construction protocol: QIAGEN RNAeasy Kappa mRNA Hyperprep kit
# Samples:
# SRR15209007	CD4
# SRR15209010	MyLa
# SRR15209008	HH
# SRR15209011	HuT78
fastq=/path/to/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump
cd $DIR
mkdir GSE181130
cd $DIR/GSE181130
for geo_id in SRR15209007 SRR15209010 SRR15209008 SRR15209011; do
sbatch -c 2 -p mid -J $geo_id -o $geo_id.out -e $geo_id.err <<EOF
#!/bin/sh
$fastq -F --gzip $geo_id
EOF
done
```
Create samples subfolders:
```
mkdir raw trim fastqc star bam plots tables snp
cd raw
mkdir CD4 MyLa HH HuT78
# link fastq files to raw folder
ln -s ../../GSE181130/SRR15209007.fastq.gz CD4/.
ln -s ../../GSE181130/SRR15209010.fastq.gz MyLa/.
ln -s ../../GSE181130/SRR15209008.fastq.gz HH/.
ln -s ../../GSE181130/SRR15209011.fastq.gz HuT78/.
```


## Starting analysis pipeline via SLURM
```
# Get sample names
DIR=/path/to/Project_folder
cd $DIR
SAMPLES=$(find $DIR/raw/* -maxdepth 1 -type d)
# samples to be iterated need to be in single folder
for SAMPLE in $SAMPLES
do
cd $SAMPLE
sample=${PWD##*/}
# trim sequences and quality control via fastQC, plots, 30-50min
RES=$(sbatch -J $sample.1 -o $sample.1.out -e $sample.1.err ../../trim.sh)
# alignment without joining sequences, 6-30min
RES2=$(sbatch --dependency=afterok:${RES##* } -J $sample.2 -o $sample.2.out -e $sample.2.err ../../align.sh)
# Remove duplicates via picard + coverage, rna metrix, insert size, 1.5h
RES3=$(sbatch --dependency=afterok:${RES2##* } -J $sample.3 -o $sample.3.out -e $sample.3.err ../../picard.sh)
# variant calling via gatk, 5-10h
RES4=$(sbatch --dependency=afterok:${RES3##* } -J $sample.4 -o $sample.4.out -e $sample.4.err ../../gatk_rna.sh)
# variant calling via varscan, 10-15h
RES5=$(sbatch --dependency=afterok:${RES3##* } -J $sample.5 -o $sample.5.out -e $sample.5.err ../../varscan.sh)
done
```


## Monitoring and deletion of transient files
Delete transient trimmed sequence files, after alignment took place.
```
DIR=/path/to/Project_folder
cd $DIR
ls -lh $DIR/trim/*
# rm $DIR/trim/*
ls -lh $DIR/raw/*/*.1.err
grep error $DIR/raw/*/*.1.err
ls -lh $DIR/raw/*/*.2.err
ls -lh $DIR/raw/*/*.3.err
# check error files and delete, if empty/non-informative
# rm $DIR/raw/*/*.1.err
# rm $DIR/raw/*/*.2.err
# rm $DIR/raw/*/*.3.err

# trimming statistics
cd $DIR/raw
grep "an average of" */*.1.o*
grep "Clipped 'end' reads" */*.1.o*
grep "Total reads" */*.1.o*

# duplicates
grep "records as duplicates." */*.4.err

# alignment statistics
cd $DIR/star
grep "mapped reads %" *.Log.final*
grep "% of reads mapped to multiple loci" *.Log.final*
grep "mapped reads n" *.Log.final*
grep "Number of input reads" *.Log.final*

# info on sequencing + mapping, star statistics
R --file="stats_seq.R"
```


## Functional annotation of snps
The Variant effect preditor (VEP, ensembl) adds annotation to the mutations:
```
DIR=/path/to/Project_folder
cd $DIR/
# functional annotation of snps via VEP (ensembl)
cd $DIR/snp
ens=/path/to/Genomes/ens/90/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# vep parameter
# --pick_allele: chooses one line of consequence data per variant allele
for vcf in $(ls *ohneChr.varscan.vcf)
do
	date; echo $vcf
	vep --species homo_sapiens --fork 4 -i $vcf --output_file $vcf.txt --cache --cache_version 90 --assembly GRCh38 --merged --offline --pick_allele -sift b --polyphen b --force --stats_file $vcf.html --symbol --af --fasta $ens --quiet --gene_phenotype --fields Uploaded_variation,Location,STRAND,Allele,Gene,SYMBOL,Existing_variation,FREQS,GMAF,Consequence,cDNA_position,Protein_position,Amino_acids,Codons,SIFT,PolyPhen,GENE_PHENO,ClinVar
done
date
```


## Summarise results
Variant location and annotation of all samples are combined in one table:
```
DIR=/path/to/Project_folder
cd $DIR
R --file="gatk_varsc_single.R" &> R_mut_single.out &
```
