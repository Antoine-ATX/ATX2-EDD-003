from itertools import combinations, product
from Bio import SeqIO
import numpy as np
import random

input_file = 'ATX2/ATX2-EDD-003/Initial_peptides/9_initial_points_from_iteration_1.fasta'

amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
# Excluding Selenocystein - very rare and not really drug-like
nb_AA = len(amino_acids)

nb_modifications_exhaustive_search = 3
nb_modifications_non_exhaustive_search = 5 
nb_peptides_non_exhaustive_search = 1000



def modify_peptides(initial_peptide: str, indexes_to_modify: tuple, new_AA: list):
    generated_peptides=[]
    for AAs in new_AA:
        modified_peptide = list(initial_peptide)
        for i in range(len(indexes_to_modify)):
            modified_peptide[int(indexes_to_modify[i])]=AAs[i]
        generated_peptides.append(''.join(modified_peptide))
    return(generated_peptides)

def exhaustive_exploration(initial_peptide: str, nb_modifications: int):
    generated_list=[]
    positions = np.linspace(0,len(initial_peptide)-1,num=len(initial_peptide))
    list_of_indexes = list(combinations(positions, nb_modifications))
    new_AA = list(product(amino_acids, repeat=nb_modifications))
    for indexes_to_modify in list_of_indexes:
        generated_peptides = modify_peptides(initial_peptide, indexes_to_modify, new_AA)
        generated_list.extend(generated_peptides)
    
    return(generated_list)

def random_exploration(initial_peptide: str, nb_modifications: int, nb_peptides: int):
    generated_list=[]
    positions = np.linspace(0,len(initial_peptide)-1,num=len(initial_peptide))
    list_of_indexes = list(combinations(positions, nb_modifications))
    new_AA = list(product(amino_acids, repeat=nb_modifications))

    while len(generated_list)<nb_peptides:
        random_indexes = random.choice(list_of_indexes)
        random_AAs = random.choice(new_AA)
        modified_peptide = list(initial_peptide)   
        for i in range(len(random_indexes)):
            modified_peptide[int(random_indexes[i])]=random_AAs[i] 
        generated_list.append(''.join(modified_peptide))
        
    generated_list = list(set(generated_list))
    
    return(generated_list)

fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

for peptide in fasta_sequences:
    name, sequence = peptide.id, str(peptide.seq)
    with open('ATX2/ATX2-EDD-003/Sequences/derived_from_'+name+'.fasta','w') as out_file:
        exhaustive_peptides = exhaustive_exploration(sequence, nb_modifications_exhaustive_search)
        non_exhaustive_peptides = random_exploration(sequence, nb_modifications_non_exhaustive_search, nb_peptides_non_exhaustive_search)
        exhaustive_peptides.extend(non_exhaustive_peptides)
        all_derived_peptides = list(set(exhaustive_peptides))
        for i in range(len(all_derived_peptides)):
            out_file.write('>parent_'+name+'_number_'+str(i)+'\n')
            out_file.write(all_derived_peptides[i]+'\n')


