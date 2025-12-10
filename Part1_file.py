#!/usr/bin/env python
# coding: utf-8




#pwd  just to make sure its in correct location 


# 1.1 Set up/Instals




#get_ipython().system('pip install Biopython')

from Bio import SeqIO    #Downloaded to be able to parse fasta files

genomes = {}   #Make a dictionary to store our genomes to use later 

files = ["NC_037282.1.fasta", "EcoliK12-MG1655.fasta"]

for filename in files: 
    record = next(SeqIO.parse(filename, "fasta")) 
    chrom_name = record.id
    sequence = str(record.seq)              #Turn the sequences into string, these will be the refs
    genomes[chrom_name] = {
            "chrom": chrom_name,
            "seq": str(record.seq),  #Keep track of sequence and sequence length 
            "len": len(record.seq)  
}   

import random #for the rng 
random.seed(12345)





print("Keys in genomes dictionary:")
print(list(genomes.keys()))

#checking the structure
for key in genomes.keys():
    print(f"\nKey: '{key}'")
    print(f"  Has 'seq': {'seq' in genomes[key]}")


# SNPs function code 




#Have to create copies of genomes in memory 
#The copies will be the mutated genomes 

def mutate_snps(original_genome, n_snps, chrom_name):    #no. of snps
    mutated_genome = list(original_genome)
    mutation_log = [] #empty for now 
    positions = random.sample(range(len(original_genome)), n_snps)
       
    for pos in positions:
        nucleotides = ["A", "T", "G", "C"]  #only accepts these chrs
        current_base = original_genome[pos]
        if current_base in nucleotides:
            nucleotides.remove(current_base)
            snp_base = random.choice(nucleotides)
            mutated_genome[pos] = snp_base
            mutation_log.append({
                "chrom": chrom_name, #log entry for vcf
                "position": pos + 1,  #for positioning 
                "ref": current_base,
                "alt": snp_base
            })
    mutated_seq_str = "".join(mutated_genome)  #the list comes back as string again 
    return mutated_seq_str, mutation_log


all_snp_mutations = []  #this code is just to make sure the right no. of snps are made 
snp_logs_by_genome = {}   
for chrom_name in genomes:
    original_seq = genomes[chrom_name]["seq"]
    mutated_seq, snp_log = mutate_snps(original_seq, 300, chrom_name)  #300 snps 
    
    genomes[chrom_name]["mutated_seq"] = mutated_seq
    snp_logs_by_genome[chrom_name] = snp_log
    all_snp_mutations.extend(snp_log)
    
    print(f"Applied 300 SNPs to {chrom_name}")
print(f"\nTotal SNPs recorded: {len(all_snp_mutations)}")

# INDEL function code 




def mutate_indels(seq, n_indels, chrom_name):
    mutated_genome = list(seq)  #converts to list so i can mutate 
    mutation_log = []
    positions = random.sample(range(len(seq) - 10), n_indels)  #zero based positions
    
    #works out all the mutations first before actually changing anything 
    planned_mutations = []
    for pos in positions:
        indel_type = random.choice(["insertion", "deletion"])
        indel_size = random.randint(1, 10)
        planned_mutations.append((pos, indel_type, indel_size))
    
    #goes in reverse to keep positioning from shifting 
    planned_mutations.sort(key=lambda x: x[0], reverse=True)
    
    for original_pos, indel_type, indel_size in planned_mutations:
        
        if indel_type == "insertion":
            new_bases = ''.join(random.choice(["A", "T", "G", "C"]) for x in range(indel_size))
            ref_base = mutated_genome[original_pos]  #uses the nucleotide before the indel as an anchor 
            
            
            mutated_genome[original_pos + 1:original_pos + 1] = list(new_bases)
            
            ref = ref_base
            alt = ref_base + new_bases
            
            mutation_log.append({
                "chrom": chrom_name,
                "position": original_pos + 1,  #same as snps 
                "type": "insertion",
                "size": indel_size,
                "ref": ref,
                "alt": alt
            })
        
        elif indel_type == "deletion":
            #will delete bases after 
            size_deleted = min(indel_size, len(mutated_genome) - original_pos - 1) #makes sure nothing past the genome is deleted 
            ref_base = mutated_genome[original_pos]
            deleted_bases = ''.join(mutated_genome[original_pos + 1:original_pos + 1 + size_deleted])
            
            ref = ref_base + deleted_bases
            alt = ref_base
            
            del mutated_genome[original_pos + 1:original_pos + 1 + size_deleted]  #removed base after the anchor 
            
            mutation_log.append({
                "chrom": chrom_name,
                "position": original_pos + 1,  
                "type": "deletion",
                "size": size_deleted,
                "ref": ref,
                "alt": alt
            })
    
    mutated_genome_str = ''.join(mutated_genome)
    return mutated_genome_str, mutation_log


# Now I have applied the indels to the genomes


indel_logs_by_genome = {}
all_indel_mutations = []

for chrom_name in genomes:
    # start from the SNP-mutated sequence
    seq_with_snps = genomes[chrom_name]["mutated_seq"]
    
    # apply 20 indels
    seq_with_indels, indel_log = mutate_indels(seq_with_snps, 20, chrom_name)
    
    # update the stored mutated sequence to include indels too
    genomes[chrom_name]["mutated_seq"] = seq_with_indels
    
    # store logs
    indel_logs_by_genome[chrom_name] = indel_log
    all_indel_mutations.extend(indel_log)
    
    print(f"Applied {len(indel_log)} indels to {chrom_name}")


#to make sure 2 files were made 

print(list(genomes.keys()))
print(len(snp_logs_by_genome["NC_037282.1"]))
print(len(snp_logs_by_genome["NC_000913.3"]))





#To check if the snps works, test on the random bases  
#skipping the verification to keep it simple 
Bases_test = "ATCGATCGATCG"  #Had to use on the valid letters this code would work on 
print("Original_bases_test:", Bases_test)
n_mutations = 3
chrom_name = "test_chromosome"

mutated_bases_test, bases_test_log = mutate_snps(Bases_test, 3, "test_chromosome")

print("Mutated_bases_test: ", mutated_bases_test) #returns the new mutated bases back 

#show what happened
print("\nMutations:")
for mutation in bases_test_log:
    print(f"  {mutation}")  


# To see if the snps function works 






#defining the genomes 
NC_snp_log = snp_logs_by_genome['NC_037282.1']
Ecoli_snp_log = snp_logs_by_genome['NC_000913.3']





print("NC genome:")
print(f"  SNPs created: {len(NC_snp_log)}")
print(f"  All unique: {len(NC_snp_log) == len(set(snp['position'] for snp in NC_snp_log))}") #another check to  make sure snps are one base and that its unique 

print("E. coli genome:")
print(f"  SNPs created: {len(Ecoli_snp_log)}")
print(f"  All unique: {len(Ecoli_snp_log) == len(set(snp['position'] for snp in Ecoli_snp_log))}")





# To test if the snps and indels code work 



#similar test to before on indels
Bases_test_2 = "ATCGATCGATCGATCGATCG"  #string has to be longer for the test to work
print("Original bases test:", Bases_test_2)
mutated_bases_test, bases_test_log = mutate_indels(Bases_test_2, 3, "test_chromosome")
print("Mutated bases test:", mutated_bases_test)

print("\nMutations:")
for mutation in bases_test_log:
    print(f"  {mutation}") #prints the variable we put through at the time i.e. the alphabet 


#according to the  test, the indels code works 



#Now i need to get 2 fasta for the mutated genomes

def write_fasta(chrom_name, sequence, filename):
    with open(filename, "w") as fasta:
        fasta.write(f">{chrom_name}\n")
        fasta.write(sequence + "\n")
        
for chrom_name in genomes:
    seq = genomes[chrom_name]["mutated_seq"]
    filename = f"{chrom_name}.mutated.fasta"
    
    write_fasta(chrom_name, seq, filename)
    print("Wrote", filename)


# Now i want to create vcf from mutation_logs of snps and indels 




#using snps + indels 
#and make it into one big list  and then makes it into a vcf 
#but first define the function 

def write_vcf(all_mutations, output_file="truth_variants.vcf"):
    all_mutations.sort(key=lambda x: (x["chrom"], x["position"]))  #sorts so the vcf is in order 
    with open(output_file, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=SimulatedMutations\n")
        vcf.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of variant (SNP/INS/DEL)\">\n")
        vcf.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n")
        vcf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")  #vcf headers
        
        for i, mut in enumerate(all_mutations, 1):  #1 variant per line
            chrom = mut["chrom"]  #from the mutation logs
            pos = mut["position"]    
            ref = mut["ref"]
            alt = mut["alt"]
            
            vcf.write(f"{chrom}\t{pos}\tvar{i}\t{ref}\t{alt}\t.\tPASS\t.\n") #writes the vcf 
    print(f" Wrote {len(all_mutations)} variants to {output_file}")  #checks if the code works 


for chrom_name in genomes:
    genome_snps = snp_logs_by_genome[chrom_name]
    genome_indels = indel_logs_by_genome[chrom_name]
    genome_mutations = genome_snps + genome_indels    #combines all the snps and indels of that specific genome into a list 
    
    genome_mutations.sort(key=lambda x: (x["chrom"], x["position"]))
    
    output_file = f"{chrom_name}.truth_variants.vcf"
    write_vcf(genome_mutations, output_file)

#then combines both genome truths 
all_mutations = all_snp_mutations + all_indel_mutations
all_mutations.sort(key=lambda x: (x["chrom"], x["position"]))

print(f"Total mutations: {len(all_mutations)}")
print(f"First mutation: {all_mutations[0]}")
print(f"Last mutation:  {all_mutations[-1]}")

write_vcf(all_mutations, "truth_variants.vcf")  #muatations from both genomes 



# Simulating sequence reads (precision/recall) 



#making fake reads 
#part 1 calculate the number of reads. this should have the coverage = 30, and the read length = 100
#define generate_reads

def generate_read(genome_seq, read_length):   #returns a random read of that length 
    
    genome_length = len(genome_seq)  #needs to know length to be able to pick random position 
    max_start = genome_length - read_length   #works out highest starting position 
    start_pos = random.randint(0, max(0, max_start))  #makes sure the random no. generated isnt neg
    seq = genome_seq[start_pos:start_pos + read_length]  #gives the actual genome_read 
    return seq, start_pos



def simulate_reads(genome_seq, read_length=100, coverage=30):  #define the parameters
#this will generate reads at specific coverage and read length
    genome_length = len(genome_seq)
    number_reads = int((coverage * genome_length) / read_length) #int to get whole no.
    #number_reads = int((30*genome_length)/100)
    reads = [] #this is for 1 chromosome

    for i in range(number_reads):
        seq, start_pos = generate_read(genome_seq, read_length) #gives back the bases starting from whatever position 
        qual = 'I' * read_length   #confidence values for the fastq
        reads.append((seq, qual))
    
    return reads

#use it for all chromosomes
def write_fastq(reads, filename):  
    with open(filename, 'w') as f:
        for i, (seq, qual) in enumerate(reads):
            f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")  #actually writes the fastq

for chrom_name in genomes:
    mutated_seq = genomes[chrom_name]["mutated_seq"]   #gets the bases from the mutated genome from a specific chromosome 
    
    print(f"Generating reads for {chrom_name}...")
    
    reads = simulate_reads(mutated_seq, read_length=100, coverage=30)
    
    print(f"  Generated {len(reads)} reads")
    
    fastq_filename = f"{chrom_name}.simulated_reads.fastq"
    write_fastq(reads, fastq_filename)
    print(f"  Wrote {fastq_filename}")
           
                     


# To test the simulate code


test_genome = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

print("Testing with alphabet:")
print(f"Genome: {test_genome}")
print(f"Genome length: {len(test_genome)}\n")

#Test generate_read
print("Testing generate_read:")
for i in range(5):
    seq, start_pos = generate_read(test_genome, read_length=5)
    print(f"  Read {i+1}: '{seq}' (starts at position {start_pos})")

print("\nTesting simulate_reads:")
reads = simulate_reads(test_genome, read_length=5, coverage=3)
print(f"Generated {len(reads)} reads")
print(f"Expected: {int((3 * 26) / 5)} reads")  #3 coverage, 26 bases, 5 read_length this is random and should end up with 15

print("\nFirst 10 reads:")
for i, (seq, qual) in enumerate(reads[:10]):
    print(f"  {i+1}. Seq: '{seq}', Qual: '{qual}'")

#Test write_fastq
print("\nTesting write_fastq:")
write_fastq(reads[:5], "test_alphabet.fastq")
print("Wrote test_alphabet.fastq - check the file!")

#Read it back to verify
print("\nContents of test_alphabet.fastq:")
with open("test_alphabet.fastq", 'r') as f:
    print(f.read())


# mapping 

# Mapping and comparing vcf, starting with Ecoli genomes 
#commands were run in terminal 




#need to now compare simulated variants with variants detected 
def parse_vcf_simple(vcf_file):
    variants = set()  #empty for now to store variants 
    
    with open(vcf_file, 'r') as f:   #open the vcf and loop each line 
        for line in f:
            if line.startswith('#'):  #skips header during parsing 
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

#loads both VCFs
truth_variants = parse_vcf_simple('NC_000913.3.truth_variants.vcf')
called_variants = parse_vcf_simple('Ecoli_variants.vcf')
#first calculates accuracy for nc_000913.3 e coli
true_positives = truth_variants & called_variants  #real mutation detected
false_positives = called_variants - truth_variants  #whatever bcftools called but isnt in truth 
false_negatives = truth_variants - called_variants  #real mutations that were missed 
#variant caller evalution 
print(f"Truth variants: {len(truth_variants)}")
print(f"Called variants: {len(called_variants)}")
print(f"True Positives: {len(true_positives)}")
print(f"False Positives: {len(false_positives)}")
print(f"False Negatives: {len(false_negatives)}")

recall = len(true_positives) / len(truth_variants) if len(truth_variants) > 0 else 0  #out of all the real mutations, how many were found, also avoided /0
precision = len(true_positives) / len(called_variants) if len(called_variants) > 0 else 0 #how many were real 
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0  #combines precision and recall into f1, the higher the better 

print(f"\nRecall: {recall:.3f}")
print(f"Precision: {precision:.3f}")
print(f"F1 Score: {f1:.3f}")


#use the same code for the other genome so this part just repeats 



def parse_vcf_simple(vcf_file):
    variants = set()  #empty for now to store variants 
    
    with open(vcf_file, 'r') as f:   #open the vcf and loop each line 
        for line in f:
            if line.startswith('#'):  #skips header during parsing 
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

#loads both VCFs
truth_variants = parse_vcf_simple('NC_037282.1.truth_variants.vcf')
called_variants = parse_vcf_simple('NC_037282_variants.vcf')

true_positives = truth_variants & called_variants  #real mutation detected
false_positives = called_variants - truth_variants  #whatever bcftools called but isnt in truth 
false_negatives = truth_variants - called_variants  #real mutations that were missed 

print(f"Truth variants: {len(truth_variants)}")
print(f"Called variants: {len(called_variants)}")
print(f"True Positives: {len(true_positives)}")
print(f"False Positives: {len(false_positives)}")
print(f"False Negatives: {len(false_negatives)}")

recall = len(true_positives) / len(truth_variants) if len(truth_variants) > 0 else 0  #out of all the real mutations, how many were found, also avoided /0
precision = len(true_positives) / len(called_variants) if len(called_variants) > 0 else 0 #how many were real 
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0  #combines precision and recall into f1, the higher the better 

print(f"\nRecall: {recall:.3f}")
print(f"Precision: {precision:.3f}")
print(f"F1 Score: {f1:.3f}")




#to compare results 
#this metric calc is repeated from previous 
import pandas as pd

#calculate accuracy metrics for E. coli
truth_variants_ecoli = parse_vcf_simple('NC_000913.3.truth_variants.vcf')
called_variants_ecoli = parse_vcf_simple('Ecoli_variants.vcf')
true_positives_ecoli = truth_variants_ecoli & called_variants_ecoli
false_positives_ecoli = called_variants_ecoli - truth_variants_ecoli
false_negatives_ecoli = truth_variants_ecoli - called_variants_ecoli

recall_ecoli = len(true_positives_ecoli) / len(truth_variants_ecoli) if len(truth_variants_ecoli) > 0 else 0
precision_ecoli = len(true_positives_ecoli) / len(called_variants_ecoli) if len(called_variants_ecoli) > 0 else 0
f1_ecoli = 2 * (precision_ecoli * recall_ecoli) / (precision_ecoli + recall_ecoli) if (precision_ecoli + recall_ecoli) > 0 else 0

#calculate metrics for Plasmid
truth_variants_plasmid = parse_vcf_simple('NC_037282.1.truth_variants.vcf')
called_variants_plasmid = parse_vcf_simple('plasmid_variants.vcf')
true_positives_plasmid = truth_variants_plasmid & called_variants_plasmid
false_positives_plasmid = called_variants_plasmid - truth_variants_plasmid
false_negatives_plasmid = truth_variants_plasmid - called_variants_plasmid

recall_plasmid = len(true_positives_plasmid) / len(truth_variants_plasmid) if len(truth_variants_plasmid) > 0 else 0
precision_plasmid = len(true_positives_plasmid) / len(called_variants_plasmid) if len(called_variants_plasmid) > 0 else 0
f1_plasmid = 2 * (precision_plasmid * recall_plasmid) / (precision_plasmid + recall_plasmid) if (precision_plasmid + recall_plasmid) > 0 else 0

#create a summary table
results_data = {
    'Genome': ['E. coli (NC_000913.3)', 'Plasmid (NC_037282.1)'],
    'Truth_Variants': [len(truth_variants_ecoli), len(truth_variants_plasmid)],
    'Called_Variants': [len(called_variants_ecoli), len(called_variants_plasmid)],
    'True_Positives': [len(true_positives_ecoli), len(true_positives_plasmid)],
    'False_Positives': [len(false_positives_ecoli), len(false_positives_plasmid)],
    'False_Negatives': [len(false_negatives_ecoli), len(false_negatives_plasmid)],
    'Precision': [precision_ecoli, precision_plasmid],
    'Recall': [recall_ecoli, recall_plasmid],
    'F1_Score': [f1_ecoli, f1_plasmid]
}

results_df = pd.DataFrame(results_data)


#display table on screen
print("PART 1: variant calling for both mutations")

print(results_df.to_string(index=False))


#save to CSV file
results_df.to_csv('part1_results.csv', index=False)
print("\n✓ Results saved to: part1_results.csv")


def evaluate_snps(truth_vcf, called_vcf, genome_name):
    """Calculate for snps only"""
    
    truth_variants = parse_vcf_simple(truth_vcf)
    called_variants = parse_vcf_simple(called_vcf)
    
    #filters for snps only
    truth_snps = {v for v in truth_variants if v[2] == 'SNP'}
    called_snps = {v for v in called_variants if v[2] == 'SNP'}
    
    #same as previous calc
    tp = truth_snps & called_snps
    fp = called_snps - truth_snps
    fn = truth_snps - called_snps
    
    precision = len(tp) / len(called_snps) if len(called_snps) > 0 else 0
    recall = len(tp) / len(truth_snps) if len(truth_snps) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {
        'Genome': genome_name,
        'Truth_SNPs': len(truth_snps),
        'Called_SNPs': len(called_snps),
        'True_Positives': len(tp),
        'False_Positives': len(fp),
        'False_Negatives': len(fn),
        'Precision': precision,
        'Recall': recall,
        'F1_Score': f1
    }

# snps for both genomes

print("How good was the SNPS caller?")


ecoli_snp_results = evaluate_snps(
    'NC_000913.3.truth_variants.vcf',
    'Ecoli_variants.vcf',
    'E. coli (NC_000913.3)'
)

plasmid_snp_results = evaluate_snps(
    'NC_037282.1.truth_variants.vcf',
    'plasmid_variants.vcf',
    'Plasmid (NC_037282.1)'
)

#create and display table
snp_results_df = pd.DataFrame([ecoli_snp_results, plasmid_snp_results])
print(snp_results_df.to_string(index=False))
print("="*80)

#save to csv
snp_results_df.to_csv('part1_snp_results.csv', index=False)
print("\n✓ SNP results saved to: part1_snp_results.csv")








