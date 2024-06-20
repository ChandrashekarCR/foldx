#! /usr/bin/python3
''' This script creates all 19 different mutations for each amino acid at every position 
    with respect to the PDB structure.
'''

import os
import pandas as pd
import argparse
import subprocess


# python Amino Acids Thee2One & One2Three Dictionary
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C', 'UNK':'X'}

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS', 'X':'UNK'}

aa_list = list(three_letter.keys())[:-1]

class mutation_list_builder:
    def __init__(self,input_file,output_file):
        self.input_file = input_file
        self.output_file = output_file

    def mutation_list(self):
        amino_acid_list=[]
        with open(self.input_file,'r') as f:
            for line in f.readlines():
                if line.startswith('ATOM') and 'CA' in line: 
                    amino_acid_list.append(one_letter[line[17:20]].strip()+line[21].strip()+line[22:26].strip())
        
        df = pd.DataFrame(amino_acid_list,columns=['aa'])
        df['wild_aa'] = df['aa'].apply(lambda x:x[0])
        df['chain'] = df['aa'].apply(lambda x:x[1])
        df['pos'] = df['aa'].apply(lambda x:x[2:]).astype(int)
        df = df.groupby(by='pos').agg({'aa':','.join}).reset_index()
        df['no_mutation'] = df['aa'].apply(lambda x:len(x.split(',')))
    
        temp = []
        for index,row in df.iterrows():
            for aa in aa_list:
                if row['aa'][0]!=aa:
                    temp.append(list(str(x+aa) for x in row['aa'].split(',')))
        
        joined_data = [",".join(sublist) for sublist in temp]
    
        with open(f'{self.output_file}/individual_list_{self.input_file.split("/")[-1].split(".")[0]}.txt','w') as f:
            for line in joined_data:
                f.write(line+';'+'\n')
    
   

if __name__=='__main__':
    parser = argparse.ArgumentParser('A script to generate mutation files')
    parser.add_argument('-i',nargs = '+',dest='input_arguement',default='something',help='help file here')
    parser.add_argument('-o',nargs = '+',dest='output_arguement',default='something',help='help file here')
    args_ = parser.parse_args()
    input_files = args_.input_arguement
    output_files = args_.output_arguement
    print('Making mutation files...')
    mutation_list_builder(input_file=input_files[0],output_file=output_files[0]).mutation_list()
    print('Done!')
    
