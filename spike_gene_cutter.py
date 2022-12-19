def spike_gene_cutter_cleaner(input_file, file_format):

    from Bio import SeqIO     # We use Biopython's SeqIO parser to load our sequences

    # The  next section loads the SARS-CoV-2 sequences from file.

    input_file = input("Enter full 'input' filename and filetype eg input_file.fasta")
    file_format = input("Enter 'file format' eg 'fasta")

    # Extraction of the approximate S gene
    spike_genes = []        # All extracted spike genes will be stored in this list

    for seq_record in SeqIO.parse(input_file, file_format):
        spike_cut = seq_record[21500:25500]
        spike_genes.append(spike_cut)
    

    spike_file = SeqIO.write(spike_genes, 'Spike_gene_file.fasta', 'fasta')
    # print('Please check current folder for the Spike gene sequences which were extracted')


    # print_sample = input('Do you want to view snips of extracted sequences: Y/N ')

    # if print_sample.upper() == 'Y':
    #     for spike in SeqIO.parse('Spike_gene_file.fasta','fasta'):
    #         print(f'Length of approximate spike_gene = {len(spike)}')
    #         print(f'Representative spike sequence = {repr(spike.seq)}')
    #         print(f'Sequence id = {spike.id}')

    # elif print_sample.upper() == 'N':
    #     pass

    # else:
    #     print("Please enter 'Y' for 'Yes' and 'N' for 'No'")


####################################################################################################
####################################################################################################













