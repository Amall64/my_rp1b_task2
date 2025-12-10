#!/bin/bash

# Exit if any command fails
set -e

echo " Genome Mutation and Variant Calling Validation "
echo ""


# Conda Environment

echo "Using already-activated conda environment: pipeline"
echo "If not active, run:  conda activate pipeline"
echo ""


# Part 1: Python validation / simulation

echo "Running Part 1: Validation Code"

python3 Part1.py

echo "Part 1 Complete"


# E. coli (NC_000913.3) variant calling

echo " Processing E. coli (NC_000913.3) "

# Index reference
samtools faidx EcoliK12-MG1655.fasta
bwa index EcoliK12-MG1655.fasta

# Align reads
bwa mem EcoliK12-MG1655.fasta NC_000913.3.simulated_reads.fastq > Ecoli_aligned.sam

# Convert to BAM, sort, and index
samtools view -bS Ecoli_aligned.sam > Ecoli_aligned.bam
samtools sort Ecoli_aligned.bam -o Ecoli_aligned_sorted.bam
samtools index Ecoli_aligned_sorted.bam

# Call variants
bcftools mpileup -f EcoliK12-MG1655.fasta Ecoli_aligned_sorted.bam \
    | bcftools call -mv -Ov -o Ecoli_variants.vcf

# Count variants
echo "E. coli variants found:"
grep -v "^#" Ecoli_variants.vcf | wc -l
echo ""

# NC_037282.1 variant calling

echo " NC_037282.1 "

# Index reference
samtools faidx NC_037282.1.fasta
bwa index NC_037282.1.fasta

# Align reads
bwa mem NC_037282.1.fasta NC_037282.1.simulated_reads.fastq > plasmid_aligned.sam

# Convert to BAM, sorting and indexing
samtools view -bS plasmid_aligned.sam > plasmid_aligned.bam
samtools sort plasmid_aligned.bam -o plasmid_aligned_sorted.bam
samtools index plasmid_aligned_sorted.bam

# Call variants
bcftools mpileup -f NC_037282.1.fasta plasmid_aligned_sorted.bam \
    | bcftools call -mv -Ov -o plasmid_variants.vcf

# Count variants
echo "Plasmid variants found:"
grep -v "^#" plasmid_variants.vcf | wc -l



# Finished

echo " Variant calling complete "
echo " Output files created:"
echo "   - Ecoli_aligned_sorted.bam"
echo "   - Ecoli_variants.vcf"
echo "   - plasmid_aligned_sorted.bam"
echo "   - plasmid_variants.vcf"

echo "Now run part2.py"

