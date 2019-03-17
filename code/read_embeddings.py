from pathlib import Path
import numpy as np
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 12:08:08 2019

@author: jjnun
"""
data_folder = Path("../data")
devine_embeds = data_folder / "devine_etal_200.txt"



def read_embedding_matrix(filename):
    ##concept_to_cui, cui_to_concept = get_cui_concept_mappings() # comment out this after fix input
    with open(filename, 'r') as infile:
        embedding_num, dimension = map(int, infile.readline().strip().split(' '))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines() 
        ##idx_to_cui = {}
        ##cui_to_idx = {}
        for idx in range(embedding_num-1):
            datum = data[idx+1].strip().split(' ')
            ##cui = datum[0]
            ##if cui[0] != 'C':
            ##    if cui in concept_to_cui:
            ##        cui = concept_to_cui[cui]
            ##embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
            embedding_matrix[idx,:] = np.array(datum[1:])
            ##idx_to_cui[idx] = cui
            ##cui_to_idx[cui] = idx
        return embedding_matrix##, idx_to_cui, cui_to_idx
    
print(read_embedding_matrix(devine_embeds))

## TO DO: Continue incorporating the Choi code. You'll want to to repurpose and translate it to Python 3.
## Probably just want to use the Pederson values because they're available lol, and compare some of those. 
## And then maybe next go on to the michigian stuff. And then by the hopefully you'll get that silly missing file and can 
## Deploy some of their stuff