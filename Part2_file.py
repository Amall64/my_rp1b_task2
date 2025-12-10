#!/usr/bin/env python
# coding: utf-8


#pwd


#make sure its in the correct directory 
import os
os.chdir('/shared/team/people/amal/Task_2/Part_1')
print(f"Working in: {os.getcwd()}")


# Instals and set up


import subprocess
import os
import shutil      #make sure this is working in notebook 
import pandas as pd

print(f"Working in: {os.getcwd()}")
print("All imports successful!")



#helper function
#a simpler code to help variant caller evalution part, this code was used in part1
def parse_vcf_simple(vcf_file):
    
    variants = set()  #empty for now to store variants 
    
    with open(vcf_file, 'r') as f:   #open the vcf and loop each line 
        for line in f:
            if line.startswith('#'):  #skip header during parsing 
                continue
            
            fields = line.strip().split('\t')   #splits the variant line into columns 
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            if len(ref) == 1 and len(alt) == 1:  #if alleles are single base = SNP
                var_type = 'SNP'
            elif len(ref) > len(alt):    #if the ref is longer = del
                var_type = 'DEL'
            elif len(ref) < len(alt):
                var_type = 'INS'
            else:
                var_type = 'COMPLEX'
            
            variants.add((chrom, pos, var_type))   #adds variants to the list 
    
    return variants




#to see what versions we have 
import subprocess


# Check bcftools
result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True)
print("bcftools:", result.stdout.split('\n')[0])

# Check snippy
result = subprocess.run(['snippy', '--version'], capture_output=True, text=True)
print("snippy:", result.stdout.strip())

# Check samtools
result = subprocess.run(['samtools', '--version'], capture_output=True, text=True)
print("samtools:", result.stdout.split('\n')[0])


# Part 2- reusable pipeline that takes reads and a reference, then finds variants using two different tools and combines the results

# defining the mapping function 



def run_bwa_mapping(reference, reads, output_prefix, reads2=None):  #maps reads with bwa 
    """Map reads to reference using BWA-MEM and process with samtools"""

    sam_file = f"{output_prefix}_aligned.sam"
    bam_file = f"{output_prefix}_aligned.bam"
    sorted_bam = f"{output_prefix}_aligned_sorted.bam"
    
    print("Indexing reference with BWA (if not already indexed)...")  #indexing
    subprocess.run(f"bwa index {reference}", shell=True, check=True)
    
    print("Running BWA-MEM...")
    
    #this makes sure it handle paired-end and single-end
    if reads2:
        print("Detected paired-end reads")
        cmd = f"bwa mem {reference} {reads} {reads2} > {sam_file}" #mapping
    else:
        print("Detected single-end reads")
        cmd = f"bwa mem {reference} {reads} > {sam_file}"
    
    subprocess.run(cmd, shell=True, check=True)  #makes a sam file
    
    print("Converting SAM to BAM...")
    subprocess.run(f"samtools view -bS {sam_file} > {bam_file}", 
                   shell=True, check=True)
    subprocess.run(f"samtools sort {bam_file} -o {sorted_bam}", 
                   shell=True, check=True)
    subprocess.run(f"samtools index {sorted_bam}", 
                   shell=True, check=True)
    
    print(f"Mapping complete: {sorted_bam}")
    return sorted_bam

print("BWA function defined")  


# Running bcftools 


def run_bcftools(reference, bam_file, output_vcf):   
    
    print("Calling variants with bcftools...")
    
    
    cmd = (
        f"bcftools mpileup -f {reference} {bam_file} | "  #ran the commands in notebook
        f"bcftools call -mv -Ov -o {output_vcf}"
    )
    
    subprocess.run(cmd, shell=True, check=True) #run bcftools command 
    
    print(f"bcftools complete: {output_vcf}")
    return output_vcf

print("bcftools function defined!")


#running snippy 



#calls for snippy from terminal 
def run_snippy(reference, reads, output_dir):
    import shutil
    import os
    
    print(f"Calling variants with snippy...")
    
    #removes output directory if it exists
    if os.path.exists(output_dir):
        print(f"  Removing existing directory: {output_dir}") #i had to remove old files to avoid errors as it wouldnt let me overwrite 
        shutil.rmtree(output_dir)
    
    
    cmd = f"snippy --outdir {output_dir} --ref {reference} --se {reads} --force" #calls snippy
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"\nSNIPPY ERROR:")
        print("STDERR:", result.stderr)
        print("STDOUT:", result.stdout)
        raise Exception("Snippy failed") #for errors
    
    snippy_vcf = f"{output_dir}/snps.vcf"
    print(f"snippy complete: {snippy_vcf}")
    return snippy_vcf

print("snippy function defined!")


#merge function defined 


def combine_vcfs(vcf1, vcf2, output_vcf): #combining both e coli vcfs 
    
    print(f"Combining VCF files with bcftools merge...")
    
    
    subprocess.run(f"bgzip -c {vcf1} > {vcf1}.gz", shell=True, check=True)  #compresses the files
    subprocess.run(f"bgzip -c {vcf2} > {vcf2}.gz", shell=True, check=True)
    subprocess.run(f"bcftools index {vcf1}.gz", shell=True, check=True) #indexing
    subprocess.run(f"bcftools index {vcf2}.gz", shell=True, check=True)
    
    cmd = f"bcftools merge {vcf1}.gz {vcf2}.gz -o {output_vcf}"  #merges
    subprocess.run(cmd, shell=True, check=True)
    
    print(f"Merge complete: {output_vcf}")
    return output_vcf

print("combine_vcfs function defined!")




#defining trust_score
def add_trust_scores(merged_vcf, output_vcf):  #this function will go through the combined vcf and add trust scores depending if both or just one caller is detected 
    
    print(f"Adding trust scores to variants...")  #input = test_pipeline_merged.vcf or the ecoli merged vcf
    
    variants_with_scores = []
    
    with open(merged_vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):  #handles vcf headers
                variants_with_scores.append(line)
                continue
            
            if line.startswith('#CHROM'):
                #when it reaches chrom adds trust score to header
                variants_with_scores.append('##INFO=<ID=TRUST,Number=1,Type=Float,Description="Trust score: 1.0=both callers, 0.5=one caller">\n')
                variants_with_scores.append(line)
                continue
            
            fields = line.strip().split('\t')
            info = fields[7]
            
            
            if len(fields) >= 11: #so the index isnt out of position
                bcftools_gt = fields[9].split(':')[0]  #extracts genotype from bcftools column in vcf (test_pipeline_aligned_sorted.bam)
                snippy_gt = fields[10].split(':')[0]   #does the same for snppy (test_pipeline_snippy)
                
                #checks if both have actual calls 
                bcftools_called = bcftools_gt not in ['./.', '.', '.|.', '.'] #which caller detected the variant 
                snippy_called = snippy_gt not in ['./.', '.', '.|.', '.']  #if its 0/1, 1/0, 1/1 at least one variant is called 
                
                if bcftools_called and snippy_called:
                    trust_score = 1.0  #both callers found this variant
                else:
                    trust_score = 0.5  #one caller found it
            else:
                trust_score = 0.5  #wont happen
            
            fields[7] = f"{info};TRUST={trust_score}"
            variants_with_scores.append('\t'.join(fields) + '\n')  #the trust scores are added to the new vcf 
    
    with open(output_vcf, 'w') as out:
        out.writelines(variants_with_scores)
    
    print(f"Trust scores added: {output_vcf}")
    return output_vcf

print("add_trust_scores function defined!")



#check for variants where only one caller found it
#counting function
print("Checking for single-caller variants...\n")
both_count = 0  #they all start at 0
bcftools_only = 0
snippy_only = 0
#this just counts what was detected by snippy and bcftools
with open("test_pipeline_merged.vcf", 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) >= 11:  #makes sure there are 11 columns, for the 2 callers 
            bcftools_gt = fields[9].split(':')[0]
            snippy_gt = fields[10].split(':')[0]
            
            bcftools_called = bcftools_gt not in ['./.', '.', '.|.']  #this time instead of assigning scores, its just counts how many are in each group 
            snippy_called = snippy_gt not in ['./.', '.', '.|.']
            
            if bcftools_called and snippy_called:
                both_count += 1
            elif bcftools_called:
                bcftools_only += 1
            elif snippy_called:
                snippy_only += 1

print(f"Both callers: {both_count}")
print(f"bcftools only: {bcftools_only}")
print(f"snippy only: {snippy_only}")




#this is just for reporting and making sure things are being counted/calc properly

def analyze_scored_vcf(vcf_file):
    """Extract trust scores from scored VCF with detailed diagnostics"""
    scores = []
    score_counts = {}
    total_variants = 0
    missing_trust = 0
    
    with open(vcf_file) as f: #reads vcf and extracts trust 
        for line in f:
            if line.startswith('#'):
                continue
            
            total_variants += 1  #for the percentage calc
            fields = line.strip().split('\t')
            info = fields[7]
            
            if 'TRUST=' in info:
                score = float(info.split('TRUST=')[1].split(';')[0])
                scores.append(score)
                score_counts[score] = score_counts.get(score, 0) + 1  #keeps track of how many times each trust score appears 
            else:
                missing_trust += 1  #or it adds it to missing_trust
    
    print("\n")
    print("DIAGNOSTIC: ACTUAL TRUST SCORES IN YOUR FILE")
    
    for score in sorted(score_counts.keys()):
        count = score_counts[score]
        pct = (count / total_variants * 100) if total_variants > 0 else 0
        print(f"  TRUST={score}: {count:,} variants ({pct:.1f}%)")
    
    if missing_trust > 0:
        print(f"  Missing TRUST tag: {missing_trust:,} variants")
    
    print(f"\n  Total variants in file: {total_variants:,}")
    print(f"  Variants with TRUST scores: {len(scores):,}")
    
    
    
    both_callers = sum(1 for s in scores if s == 1.0)
    one_caller = sum(1 for s in scores if s == 0.5)
    low_quality = sum(1 for s in scores if s == 0.0)
    mean_score = sum(scores) / len(scores) if scores else 0
    
    print("\nSummary:")
    print(f"  Both callers (TRUST=1.0): {both_callers:,}")
    print(f"  One caller (TRUST=0.5):   {one_caller:,}")
    print(f"  Low quality (TRUST=0.0):  {low_quality:,}")
    print(f"  Mean trust score:         {mean_score:.2f}")
    
    return scores


if __name__ == "__main__":
    scores = analyze_scored_vcf('test_pipeline_scored.vcf')


# Now everything is defined 

# Main pipeline function 



def run_pipeline(reference, reads, output_prefix, reads2=None):
    
    print(f"Starting variant calling pipeline...")
    print(f"Reference: {reference}")
    print(f"Reads: {reads}")
    if reads2:
        print(f"Reads2: {reads2}")
    print(f"Output prefix: {output_prefix}")
    print("-")
    
    print("\nStep 1: Mapping reads...")
    bam_file = run_bwa_mapping(reference, reads, output_prefix, reads2)  #pass reads2 here. alignment using BWA and produces a sorted BAM file
    
    print("\nStep 2: Calling variants with bcftools...")
    bcftools_vcf = run_bcftools(reference, bam_file, f"{output_prefix}_bcftools.vcf") #uses the bam  file to call variants and writes them to a vcf file
    
    print("\nStep 3: Calling variants with snippy...")
    snippy_vcf = run_snippy(reference, reads, f"{output_prefix}_snippy") #runs Snippy on the reference and reads
    
    print("\nStep 4: Combining VCF files...")
    merged_vcf = combine_vcfs(bcftools_vcf, snippy_vcf, f"{output_prefix}_merged.vcf")
    
    print("\nStep 5: Adding trust scores...")
    scored_vcf = add_trust_scores(merged_vcf, f"{output_prefix}_scored.vcf")
    
    print("\n")
    print("Pipeline complete!")
    print("=")
    
    results = {
        "bam": bam_file,
        "bcftools_vcf": bcftools_vcf,
        "snippy_vcf": snippy_vcf,
        "merged_vcf": merged_vcf,
        "scored_vcf": scored_vcf
    }
    
    return results

print("Main pipeline function defined")


# Now i can run this on my simulated data


#ecoli test 
#maps the reads to the reference
results = run_pipeline(
    reference="EcoliK12-MG1655.fasta", #original genome 
    reads="NC_000913.3.simulated_reads.fastq",  #simulated reads with ~320 mutations 
    output_prefix="test_pipeline"
)

print("\nOutput files:")
for key, path in results.items():
    print(f"  {key}: {path}")


# The pipeline has made: bam: test_pipeline_aligned_sorted.bam, bcftools_vcf: test_pipeline_bcftools.vcf, 
# snippy_vcf: test_pipeline_snippy/snps.vcf, merged_vcf: test_pipeline_merged.vcf, scored_vcf: test_pipeline_scored.vcf



#code for a validation analysis 

#title
print("PART 2: PIPELINE VALIDATION ON SIMULATED DATA")


#compare pipeline results against ground truth
truth_variants = parse_vcf_simple('NC_000913.3.truth_variants.vcf')
pipeline_variants = parse_vcf_simple('test_pipeline_scored.vcf')

#calculate true pos, false pos and false neg
tp = truth_variants & pipeline_variants
fp = pipeline_variants - truth_variants
fn = truth_variants - pipeline_variants

precision = len(tp) / len(pipeline_variants) if len(pipeline_variants) > 0 else 0
recall = len(tp) / len(truth_variants) if len(truth_variants) > 0 else 0
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

#create table to show results clearly
validation_data = {
    'Dataset': ['Simulated E. coli'],
    'Truth_Variants': [len(truth_variants)],
    'Called_Variants': [len(pipeline_variants)],
    'True_Positives': [len(tp)],
    'False_Positives': [len(fp)],
    'False_Negatives': [len(fn)],
    'Precision': [f"{precision:.3f}"],
    'Recall': [f"{recall:.3f}"],
    'F1_Score': [f"{f1:.3f}"]
}

validation_df = pd.DataFrame(validation_data)
print(validation_df.to_string(index=False))
print("=")

#save results as csv
validation_df.to_csv('part2_validation_results.csv', index=False)
print("\n✓ Validation results saved to: part2_validation_results.csv\n")


# Comparison by first converting to csv



import pandas as pd
#to make results easier to analyse
def vcf_to_csv(vcf_file, output_csv):
    """Convert VCF file to CSV for easier analysis"""
    
    variants = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            #skip header lines
            if line.startswith('##'):
                continue
            
            #get column names from #CHROM line
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                continue
            
            #parse variant lines
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                variant = {
                    'CHROM': parts[0],
                    'POS': parts[1],
                    'ID': parts[2],
                    'REF': parts[3],
                    'ALT': parts[4],
                    'QUAL': parts[5],
                    'FILTER': parts[6],
                    'INFO': parts[7]
                }
                
                #extract trust score if present
                if 'TRUST=' in parts[7]:
                    trust = [x for x in parts[7].split(';') if 'TRUST=' in x][0]
                    variant['TRUST'] = trust.split('=')[1]
                else:
                    variant['TRUST'] = 'N/A'
                
                variants.append(variant)
    
    #create DataFrame
    df = pd.DataFrame(variants)
    
    #convert POS to integer for sorting
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    
    #save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} variants to {output_csv}")
    
    return df

#convert vcfs to csv (simulated data)
print("Converting VCFs to CSV...")

bcftools_df = vcf_to_csv("test_pipeline_bcftools.vcf", "bcftools_variants.csv")
snippy_df = vcf_to_csv("test_pipeline_snippy/snps.vcf", "snippy_variants.csv")
merged_df = vcf_to_csv("test_pipeline_merged.vcf", "merged_variants.csv")
scored_df = vcf_to_csv("test_pipeline_scored.vcf", "scored_variants.csv")
truth_df = vcf_to_csv("NC_000913.3.truth_variants.vcf", "truth_variants.csv")

print("\nDone! CSV files created.")





#compare results
#summary 
print("=== VARIANT COUNTS ===")
print(f"Truth variants:     {len(truth_df)}")
print(f"bcftools called:    {len(bcftools_df)}")
print(f"snippy called:      {len(snippy_df)}")
print(f"Merged (total):     {len(merged_df)}")
print(f"Final scored:       {len(scored_df)}")

#trust score distribution
print("\n=== TRUST SCORE DISTRIBUTION ===")
trust_counts = scored_df['TRUST'].value_counts()
print(trust_counts)

# Show first few variants side by side
print("\n=== First 10 Scored Variants ===")
print(scored_df[['CHROM', 'POS', 'REF', 'ALT', 'TRUST']].head(10))

# Find variants only called by one tool
print("\n=== Low Confidence Variants (TRUST=0.5) ===")
low_conf = scored_df[scored_df['TRUST'] == '0.5']
print(f"Total: {len(low_conf)}")
print(low_conf[['CHROM', 'POS', 'REF', 'ALT']].head(10))

# Find variants called by both tools
print("\n=== High Confidence Variants (TRUST=1.0) ===")
high_conf = scored_df[scored_df['TRUST'] == '1.0']
print(f"Total: {len(high_conf)}")
print(high_conf[['CHROM', 'POS', 'REF', 'ALT']].head(10))



#running the pipeline on both SRR25083113 fastq , EcoliK12-MG1655.fasta will be the ref
#run_pipeline function this time on real data

results = run_pipeline(
    reference="EcoliK12-MG1655.fasta",  #current directory
    reads="../SRR25083113_1.fastq.gz",  # in parent directory
    reads2="../SRR25083113_2.fastq.gz",  # in parent directory
    output_prefix="Ecoli_real_pair"
)

print("\nOutput files:")
for key, path in results.items():
    print(f"  {key}: {path}")



#now i have Ecoli_real_pair_scored.vcf to compare variants scores 

real_scored_df = vcf_to_csv("Ecoli_real_pair_scored.vcf", "real_scored_variants.csv")


#shows the real data analysis
print("PART 2: REAL E. COLI VARIANT ANALYSIS")


#analyze the scored VCF 
scores = analyze_scored_vcf('Ecoli_real_pair_scored.vcf')



#last part for tview verification
import pandas as pd
real_scored_df = pd.read_csv('real_scored_variants.csv')


print("POSITIONS FOR MANUAL tview INSPECTION")


#5 high-confidence positions
print("\n=== Inspect these TRUST=1.0 positions (HIGH CONFIDENCE) ===")
high_conf_df = real_scored_df[real_scored_df['TRUST'] == 1.0].head(5)
print(high_conf_df[['POS', 'REF', 'ALT', 'TRUST']])

high_positions = high_conf_df['POS'].tolist()
print("\nCommands to run:")
for pos in high_positions:
    print(f"samtools tview Ecoli_real_pair_aligned_sorted.bam EcoliK12-MG1655.fasta -p NC_000913.3:{pos}")

#5 low-confidence positions  
print("\n Inspect these TRUST=0.5 positions (LOW CONFIDENCE) ")
low_conf_df = real_scored_df[real_scored_df['TRUST'] == 0.5].head(5)
print(low_conf_df[['POS', 'REF', 'ALT', 'TRUST']])

low_positions = low_conf_df['POS'].tolist()
print("\nCommands to run:")
for pos in low_positions:
    print(f"samtools tview Ecoli_real_pair_aligned_sorted.bam EcoliK12-MG1655.fasta -p NC_000913.3:{pos}")


print("Results are discussed in disscussion after inspecting visually.")

