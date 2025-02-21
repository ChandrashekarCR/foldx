#! /usr/bin/python3
''' This script creates all 19 different mutations for each amino acid at every position 
    with respect to the PDB structure.  It generates a file listing all possible single amino acid mutations.
'''

import os
import pandas as pd
import argparse
import subprocess


# Dictionaries for converting between one-letter and three-letter amino acid codes
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C', 'UNK':'X'}

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS', 'X':'UNK'}

aa_list = list(three_letter.keys())[:-1] # List of one-letter amino acid codes, excluding 'X' (UNK)

class mutation_list_builder:
    def __init__(self,input_file,output_file):
        self.input_file = input_file  # Path to the input PDB file
        self.output_file = output_file # Path to the output directory

    def mutation_list(self):
        amino_acid_list=[] # List to store information about each amino acid residue
        with open(self.input_file,'r') as f:
            for line in f.readlines():
                # Extract amino acid information from ATOM records where alpha-carbon (CA) is present
                if line.startswith('ATOM') and 'CA' in line: 
                    # Create a string representation of the amino acid: One-letter code + Chain ID + Residue Number
                    amino_acid_list.append(one_letter[line[17:20]].strip()+line[21].strip()+line[22:26].strip())
        
        # Create a Pandas DataFrame for easier manipulation of the amino acid data
        df = pd.DataFrame(amino_acid_list,columns=['aa'])
        df['wild_aa'] = df['aa'].apply(lambda x:x[0]) # Extract the wild-type amino acid
        df['chain'] = df['aa'].apply(lambda x:x[1]) # Extract the chain ID
        df['pos'] = df['aa'].apply(lambda x:x[2:]).astype(int) # Extract and convert residue number to integer
        
        # Group by residue number to handle cases where multiple chains have the same residue number
        df = df.groupby(by='pos').agg({'aa':','.join}).reset_index()
        df['no_mutation'] = df['aa'].apply(lambda x:len(x.split(','))) # Count how many chains have the same residue

        temp = [] # List to store all possible mutations
        for index,row in df.iterrows():
            for aa in aa_list: # Iterate through all possible amino acids
                if row['aa'][0]!=aa: # Ensure we are not mutating to the same amino acid
                    # Create a list of mutated residues for each original residue at a given position
                    temp.append(list(str(x+aa) for x in row['aa'].split(',')))
        
        # Join the list of mutated residues into a comma-separated string
        joined_data = [",".join(sublist) for sublist in temp]
    
        # Write the list of mutations to the output file.  Each line represents a position and all possible mutations at that position.
        with open(f'{self.output_file}/individual_list_{self.input_file.split("/")[-1].split(".")[0]}.txt','w') as f:
            for line in joined_data:
                f.write(line+';'+'\n')
    
   

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='foldx_mutations.py',usage='python3 <path to script> -i <pdb_file> -o <path to output directory>',description='A script to generate mutation files for the given PDB structure')
    parser.add_argument('-i',nargs = '+',dest='pdb_file',help='Path to the input PDB file')
    parser.add_argument('-o',nargs = '+',dest='output_directory',help='Path to the output directory')
    args_ = parser.parse_args()
    input_files = args_.pdb_file
    output_files = args_.output_directory
    print('Making mutation files...')
    mutation_list_builder(input_file=input_files[0],output_file=output_files[0]).mutation_list() # Create and use instance of the class
    print('Done!')