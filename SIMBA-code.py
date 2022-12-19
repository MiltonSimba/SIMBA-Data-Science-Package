from Bio import SeqIO

ref_seq = SeqIO.read('GISAID Ref Seq.fasta', 'fasta')    # reading the GISAID Reference Sequence into memory

def spike_gene_cutter_cleaner(input_file, file_format):

    # from Bio import SeqIO    # We use Biopython's SeqIO parser to load our sequences
    import sys

    # The  next section loads the SARS-CoV-2 sequences from file.

    input_file = input("Enter full 'input' filename and filetype eg input_file.fasta")
    file_format = input("Enter 'file format' eg fasta")

    # Extraction of the approximate S gene
    spike_genes = []        # All extracted spike genes will be stored in this list

    for seq_record in SeqIO.parse(input_file, file_format):
        spike_cut = seq_record[21500:25500]
        spike_genes.append(spike_cut)
    

    spike_file = SeqIO.write(spike_genes, 'Spike_gene_file.fasta','fasta')

    # Next part removes spike genes with a number of Unidentified nucleotides more than 2

    clean_spike_genes = [] # This list will store cleaned spike genes

    for spike_record in SeqIO.parse('Spike_gene_file.fasta','fasta'):
        # first its wiser to use a string object of the sequence for counting additionally all letters should be ensured to be in upper case

        string_spike = (str(spike_record.seq)).upper()

        if string_spike.count('N') > 2:              #Unidentified nucleotides more than 2
            if spike_record.seq not in clean_spike_genes:  # Removing repeat sequences
                clean_spike_genes.append(spike_record)      #spike record is appended instead of string_spike so it does not lose Sequence object attributes


print('Please check current folder for the Spike gene sequences which were extracted')


# In case the user prefers to see extracted sequences
print_sample = input('Do you want to view snips of extracted sequences: Y/N ')

if print_sample.upper() == 'Y':
    for spike in SeqIO.parse('Spike_gene_file.fasta','fasta'):
            print(f'Length of approximate spike_gene = {len(spike)}')
            print(f'Representative spike sequence = {repr(spike.seq)}')
            print(f'Sequence id = {spike.id}')

elif print_sample.upper() == 'N':
        pass

else:
    print("Please enter 'Y' for 'Yes' and 'N' for 'No'")

# User can use alignment tool of choice and it is recommended to align against the GISAID Ref Seq for more accurate
# downstream processes


####################################################################################################
####################################################################################################


aligned_sequences = []
def triplet_code_table(triplet_code):

    '''Triplet codes that code for specific amino acids based on the DNA triplet code table'''

    for triplet_code in aligned_sequences:

        # Gene sequence is a fasta file from a codon-aligned sequences, loaded using AlignIO of Biopython.
        # Reminder that all loaded DNA sequences must be in upper case and in string format

        if triplet_code == 'TTT' or 'TTC':
            amino_acid = 'F'

        elif triplet_code == 'TTA' or 'TTG' or 'CTT' or 'CTC'or 'CTA' or 'CTG':
            amino_acid = 'L'

        elif triplet_code == 'ATT' or 'ATC' or 'ATA':
            amino_acid = 'I'

        elif triplet_code == 'ATG':
            amino_acid = 'M'

        elif triplet_code == 'GTT' or 'GTC' or 'GTA' or 'GTG':
            amino_acid = 'V'

        elif triplet_code == 'TCT' or 'TCC' or 'TCA' or 'TCG' or 'AGT' or 'AGC':
            amino_acid = 'S'

        elif triplet_code == 'CCT' or 'CCC' or 'CCA' or 'CCG':
            amino_acid = 'P'

        elif triplet_code == 'ACT' or 'ACC' or 'ACA' or 'ACG':
            amino_acid = 'T'

        elif triplet_code == 'GCT' or 'GCC' or 'GCA' or 'GCG':
            amino_acid = 'A'

        elif triplet_code == 'TAT' or 'TAC':
            amino_acid = 'Y'
        
        elif triplet_code == 'TAA' or 'TAG' or 'TGA':
            amino_acid = '*'

        elif triplet_code == 'CAT' or 'CAC':
            amino_acid = 'H'

        elif triplet_code == 'CAA' or 'CAG':
            amino_acid = 'Q'

        elif triplet_code == 'AAA' or 'AAG':
            amino_acid = 'K'

        elif triplet_code == 'GAA' or 'GAG':
            amino_acid = 'E'

        elif triplet_code == 'TGT' or 'TGC':
            amino_acid = 'C'

        elif triplet_code == 'TGG':
            amino_acid = 'W'

        elif triplet_code == 'CGT' or 'CGC' or 'CGA' or 'CGG' or 'AGA' or 'AGG':
            amino_acid = 'R'

        elif triplet_code == 'GGT' or 'GGC' or 'GGA' or 'GGG':
            amino_acid = 'G'
        
        elif 'N' in triplet_code:
            amino_acid = '?'

        else:
            print('Please check if all bases are A, C, T or G')

############################################################################################################################################
base = ref_seq

def sum_of_codon_synonymous_substitutions(triplet_code, lower=0, upper=3):
    for i in range(lower=0, upper=3):      #Lower represents the lower bound in our sum and upper is the upper bound. For python index, the first value would be 0
        
        if i = 0:                      
            if base != ref_seq:            
                fraction_of_change1 = 0.05            #Fraction of change denotes the chances of a codon for coding another amino acid when there is a mutation eg 
            else:
                fraction_of_change1= 0                                    #changing the first codon nucleotide results in 5% chance of change in amino acid coded.                     
        
        elif i = 1:
            if base != ref_seq:
                fraction_of_change2 = 1
            else:
                fraction_of_change2 = 0

        elif i = 2:
            if base != ref_seq:
                fraction_of_change3 = 0.72
            else:
                fraction_of_change3 = 0

        sum_of_synonymous_fractions = sum(fraction_of_change1 + fraction_of_change2 + fraction_of_change3)              #gives use the sum of fraction of synonymous changes at ith position.
        
        return sum_of_synonymous_fractions


def sum_of_non_synonymous_substitutions(sum_of_synonymous_fractions):
    sum_of_non_synonymous_fractions = 3 - sum_of_synonymous_fractions   # to find sum of non-synonymous fraction we subtract the sum of 
    return sum_of_non_synonymous_fractions                               #synonymous fraction from 3



