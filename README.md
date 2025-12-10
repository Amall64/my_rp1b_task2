# my_rp1b_task2
# RP1B- TASK2
This repository contains code and a pipeline for:

## Part 1 – Validation Code
- Code that makes SNPs and indels into reference genomes.
- Simulating 30 x depth of perfect 100bp reads from the mutated genomes.
- Calling variants using BWA + bcftools.
- Calculating precision and recall of the SNP/indel caller.

## Part 2 – Variant Calling Pipeline
- Mapping reads and running two variant callers (bcftools and snippy).
- Combining the results into one VCF.
- Evaluating on simulated data and running once on real E. coli data.
- Assigning a confidence score to each VCF record and checking support with tview.

### Data Location
Part 1 
- **NC_037282.1.fasta**   
  Path: `/home/jovyan/ shared-team/people/amal/Task_2/Part_1/NC_037282.1.fasta`
- **EcoliK12-MG1655.fasta**   
  Path: `/home/jovyan/ shared-team/people/amal/Task_2/Part_1/EcoliK12-MG1655.fasta`
Part 2  
- **Reference:** Same reference genome location: `/home/jovyan/ shared-team/people/amal/Task_2/Part_1/EcoliK12-MG1655.fasta`
- **Reads:** `shared-team/people/amal/Task_2`
  - `SRR25083113_1.fastq.gz`
  - `SRR25083113_2.fastq.gz`

### Install the following environment:
Make sure the current environment is : cd /shared/team/people/amal/Task_2/Part_1
- conda env create -f environment.yml
- conda activate pipeline
To run the scripts in a terminal:

- chmod +x Validation_script_done.sh (to make .sh executable)
- bash -x Validation_script_done.sh (to run)
- chmod +x Pipeline_script_done.sh (executable)
- bash -x Pipeline_script_done.sh (to run)

For snippy to work make sure you have the following versions: 
- samtools –version: samtools 1.20
- snippy –version: snippy 4.6.0 

### Required tools:
- BWA-MEM (for read mapping)
- samtools (for BAM processing)
- bcftools (variant calling)
- snippy (variant calling)

## Part 1 – Validation Code: Mutating the Genome

Running Part 1
python Part1_complete.py

What this does:
This script codes to simulate mutations to 2 reference genomes and validates their variant calling accuracy.
1.	First it loads the reference genomes (E. coli NC_000913.3 and Plasmid NC_037282.1)
2.	Introduces mutations: 
•	300 SNPs per genome
•	20 indels per genome (1-10 bp each)
Outputs: 
-	NC_037282.1.mutated
-	NC_000913.3.mutated
3.	Simulated reads: 
o	100bp read length
o	30x coverage
Outputs:
-	NC_037282.1.simulated_reads.fastq
-	NC_000913.3.simulated_reads.fastq
4.	Mapping with BWA and bcftools

```bash
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
``` 
5.	Lastly, I evaluated the variant calling performance:
- Compares truth variants against called variants
- Calculates precision, recall, and F1 scores

Output:
- part1_results.csv - Summary table with precision, recall, and F1 scores

This showed the total for both snps and indels

6.	For snps specifically, the output: part1_snp_results.csv 

## Part 2 – Variant Calling Pipeline
Running Part 2
python Part2_complete.py

What this does:
This script runs a complete variant calling pipeline using two independent variant callers and combines their results: 
Pipeline Steps
1.	Maps reads to reference using (BWA-MEM)
Outputs: Ecoli_aligned_sorted.BAM 
2.	Variant Calling using bcftools
3.	Variant Calling with snippy
4.	Merging VCF Files
- Compresses using bgzip
- Merges using bcftools merge 
5.	Trust Score Assignment
   
o	Assigns trust scores to each variant: 
- TRUST=1.0: If the variant is called by both bcftools AND snippy
- TRUST=0.5: Variant called by only one caller 

Assigning a trust score explanation 
After combining snippy and bcftools VCFs, the pipeline inspects what type of tool is in column 9 and 10 of the vcf. If both callers are present, it will assign 1.0, if only one is present it will assign 0.5. 

Outputs for part 2:
- test_pipeline_aligned_sorted.bam - Mapped reads
- test_pipeline_bcftools.vcf - bcftools variants
- test_pipeline_snippy/snps.vcf – directory for snippy and variants
- test_pipeline_merged.vcf - Combined variants
- test_pipeline_scored.vcf - Final VCF with trust scores
- part2_validation_results.csv – validation calcs from comparing pipeline against ground truth


For real E. coli data (Ecoli_real_pair):
- Ecoli_real_pair_aligned_sorted.bam - Mapped reads
- Ecoli_real_pair_bcftools.vcf - bcftools variants
- Ecoli_real_pair_snippy/snps.vcf – snippy directory and variants
- Ecoli_real_pair_merged.vcf - Combined variants
- Ecoli_real_pair_scored.vcf - Final VCF with trust scores
- real_scored_variants.csv – from ecoli variants 


  
