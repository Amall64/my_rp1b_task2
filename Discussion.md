# Discussion
## Part 1 results 

| Genome | Truth Variants | Called Variants | True Positives | False Positives | False Negatives | Precision | Recall | F1 Score |
|--------|----------------|------------------|----------------|------------------|------------------|-----------|--------|----------|
| **E. coli (NC_000913.3)** | 320 | 313 | 310 | 3 | 10 | 0.990 | 0.969 | 0.979 |
| **Plasmid (NC_037282.1)** | 320 | 319 | 314 | 5 | 6 | 0.984 | 0.981 | 0.983 |

- Showing mutation metrics from indels and snps
- Across both genomes, precision, recall and f1 scores were high
- Some variants may have been missed due to sensitivity issues

  

| Genome | Truth SNPs | Called SNPs | True Positives | False Positives | False Negatives | Precision | Recall | F1 Score |
|--------|------------|--------------|----------------|------------------|------------------|-----------|--------|----------|
| **E. coli (NC_000913.3)** | 300 | 295 | 295 | 0 | 5 | 1.000 | 0.983 | 0.992 |
| **Plasmid (NC_037282.1)** | 300 | 299 | 299 | 0 | 1 | 1.000 | 0.997 | 0.998 |

- Precision is perfect for both genomes scoring 1.0 for both, and high recall for both
- These metrics suggest snp calling performace was accurate and sensitive

## Part 2 results 
### Variant calls on simulated data 

| Pipeline | Precision | Recall | F1 Score | Interpretation |
|----------|-----------|--------|----------|----------------|
| **Part 1: bcftools only** | 100% | 92.2% | 95.9% | Very accurate but misses variants |
| **Part 2: bcftools + Snippy** | 99.0% | 96.9% | 97.9% | Much better balance; highest overall accuracy |


These results are comparing part1_results.csv to part2_validation_results.csv
- Part one only used bcftools and generated perfect precision, but it has a lower recall score compared to using both callers. Suggesting some variants were missed. 
- Using both callers generated more balanced results with high scores for both recall and precision. This suggests that using both callers can detect more variants, so would make a more reliable workflow. 

### Variant calls on real data 

```
DIAGNOSTIC: ACTUAL TRUST SCORES IN YOUR FILE
TRUST=0.5: 18,412 variants (25.3%)
TRUST=1.0: 54,361 variants (74.7%)

Total variants in file: 72,773
Variants with TRUST scores: 72,773

Summary:
Both callers (TRUST=1.0): 54,361
One caller (TRUST=0.5): 18,412
Low quality (TRUST=0.0): 0
Mean trust score: 0.87

```

From analysing Ecoli_real_pair_scored.vcf
- This data is different due to the simulated data having known mutations and it having no repeats/mapping ambiguities.
- Bcftools and snippy both agree on 74.7% of the variants detected. 
- The 25.3% score could be because only one caller agreed, this could be due to the real data has more biological diversity and may have more repetitive regions present. 

### Tview verification 

After running the following command: 

samtools tview Ecoli_real_pair_aligned_sorted.bam EcoliK12-MG1655.fasta -p NC_000913.3:58

I inspected the alignment at specific locations of the genome. This gave insight into coverage and if the reads match the reference or alternative allele. 

This confirms, for high confidence variants:
- Both variant callers agreed on single based substitution (snps)
- Consistent alignment patterns 
- The location of the variant is in non repetitive regions and good quality regions 

For low confidence variants (one caller):
- They are in low quality regions where coverage may fluctuate 
- Alignment shifts from indels can confuse the callers 
- As some reads have the alt base, and majority match the ref, there seems to inconsistent/weak support for the variant 

### Strengths/what worked and limiations of the pipeline 

Strengths:
- Bcftools and snippy both detected snps very well, with perfect precision scores and high recall scores.
- The combined vcf meant no variants are missed by using just one caller, so is better. The trust score helps interpretate how confident you could be in the results
- The results from the real E coli data show that each step of the workflow worked and high vs low confidence variants behaved as expected

Limitations:
- Minimap failed so had to be replaced with BWA-MEM.
- Results show lower recall score for indels, suggesting it is less reliable for indels
- Unlike the simulated data where known muatations were made, it is difficult to have good accuracy as there is no truth set for the real E coli data

  ### To conclude: 
- The entire variant calling pipeline worked successfully to call variants and evaluated the results on both simulated and real data 
- Combining the variant callers helped improve coverage by reducing number of missed mutations 
- This task showed how applying multiple variant callers and using samtools tview to visualise key locations can improve mutation interpretation. 




