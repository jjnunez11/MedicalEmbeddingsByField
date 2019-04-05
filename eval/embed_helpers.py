# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:08:59 2019

@author: jjnun
"""
from __future__ import division
##import argparse
import numpy as np
##import scipy as sp
##from scipy.spatial.distance import cdist, cosine
##from scipy.stats import spearmanr
from icd9 import ICD9
from pathlib import Path
##import re
from cui_icd9_helpers import get_cui_concept_mappings, get_icd9_cui_mappings


tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")

# JJN changed to support csv embedding files
def read_embedding_cui(filename):
    concept_to_cui, cui_to_concept = get_cui_concept_mappings() # comment out this after fix input
    cui_to_icd9, icd9_to_cui = get_icd9_cui_mappings()
    
    # Assume embedding files with .txt are deliminated with ' ' and ',' if .csv
    if filename.endswith('.txt'):
        delim = ' '
    elif filename.endswith('.csv'):
        delim = ','
    else:
        raise Exception('embedding file must be .txt or .csv depending on deliminator')
    
    with open(str(data_folder / filename), 'r') as infile:
        first_line = infile.readline().strip().split(delim)
        embedding_num = int(first_line[0])
        dimension = int(first_line[1])
        #embedding_num, dimension = map(int, infile.readline().strip().split(' '))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines() 
        idx_to_name = {}
        name_to_idx = {}
        embedding_type_to_indices = {}
        embedding_type_to_indices['IDX'] = []
        embedding_type_to_indices['O'] = []
        for idx in xrange(embedding_num-1):
            datum = data[idx+1].strip().split(delim)
            cui = datum[0]
            if cui[0] != 'C':
                if cui in concept_to_cui:
                    cui = concept_to_cui[cui].strip() # JJN Added to remove \n
            ##print len(datum[1:]) TODO REMOVE
            embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
            # potential bug here
            if cui in cui_to_icd9:
                idx_to_name[idx] = cui_to_icd9[cui]
                name_to_idx[cui_to_icd9[cui]] = idx
                embedding_type_to_indices['IDX'].append(idx)
            else:
                idx_to_name[idx] = cui
                name_to_idx[cui] = idx
                embedding_type_to_indices['O'].append(idx)
        return embedding_matrix, embedding_type_to_indices, name_to_idx, idx_to_name
    

def read_embedding_codes(filename):
    with open(str(data_folder / filename), 'r') as infile:
        embedding_num, dimension = map(int, infile.readline().strip().split(' '))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines()
        embedding_type_to_indices = {}
        name_to_idx = {}
        idx_to_name = {}
        for idx in xrange(embedding_num-1):
            datum = data[idx+1].strip().split(' ')
            embedding_name = datum[0]
            embedding_type, embedding_value = embedding_name.split('_')
            name_to_idx[embedding_value] = idx
            idx_to_name[idx] = embedding_value 
            if embedding_type not in embedding_type_to_indices:
                embedding_type_to_indices[embedding_type] = []
            embedding_type_to_indices[embedding_type].append(idx)
            embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
        return embedding_matrix, embedding_type_to_indices, name_to_idx, idx_to_name

def generate_overlapping_sets_icd9(filenames_type):
    embedding_idx_icd9 = {} # a dictionary of (embedding_matrix, idx_to_icd9, icd9_to_idx)
    overlapping_icd9s = set([])
    start = 1

    for filename, embedding_type, _ in filenames_type:
        #print filename
        #print embedding_type
        if embedding_type == 'codes':
            embedding_matrix, embedding_type_to_indices, icd9_to_idx, idx_to_icd9 = read_embedding_codes(filename)
            embedding_idx_icd9[filename] = (embedding_matrix, idx_to_icd9, icd9_to_idx)
            if start == 1:
                start = 0
                overlapping_icd9s.update(set(icd9_to_idx.keys()))
            else:
                overlapping_icd9s.intersection_update(set(icd9_to_idx.keys()))
                #print len(overlapping_icd9s)
        elif embedding_type == 'cui':
            embedding_matrix, embedding_type_to_indices, icd9_to_idx, idx_to_icd9 = read_embedding_cui(filename)
            embedding_idx_icd9[filename] = (embedding_matrix, idx_to_icd9, icd9_to_idx)
            if start == 1:
                start = 0
                overlapping_icd9s.update(set(icd9_to_idx.keys()))
            else:
                overlapping_icd9s.intersection_update(set(icd9_to_idx.keys()))
                #print len(overlapping_icd9s)
    overlapping_icd9s = list(overlapping_icd9s)

    idx_of_overlapping_icd9s = {}
    for filename, embedding_type, _ in filenames_type:
        idx_of_overlapping_icd9s[filename] = []

    idx_to_icd9 = {}
    icd9_to_idx = {}
    for idx, icd9 in enumerate(overlapping_icd9s):
        idx_to_icd9[idx] = icd9
        icd9_to_idx[icd9] = idx
        for filename, embedding_type, _ in filenames_type:
            idx_of_overlapping_icd9s[filename].append(embedding_idx_icd9[filename][2][icd9])

    filename_to_embedding_matrix = {}
    for filename, embedding_type, _ in filenames_type:
        idx_of_overlapping_icd9s[filename] = np.array(idx_of_overlapping_icd9s[filename]).astype(int) #JJN added due to change in numpy
        filename_to_embedding_matrix[filename] = embedding_idx_icd9[filename][0][idx_of_overlapping_icd9s[filename]]
    return filename_to_embedding_matrix, idx_to_icd9, icd9_to_idx

def generate_overlapping_sets_cui(filenames_type):    
    filenames = [filename for filename, _, _2 in filenames_type]
    
    embedding_idx_cui = {} # a dictionary of (embedding_matrix, idx_to_cui, cui_to_idx)
    overlapping_cuis = set([])

    if len(filenames) == 1:
        embedding_matrix, idx_to_cui, cui_to_idx = read_embedding_matrix_cui(filenames[0])
        filename_to_embedding_matrix = {}
        filename_to_embedding_matrix[filenames[0]] = embedding_matrix
        return filename_to_embedding_matrix, idx_to_cui, cui_to_idx

    for fileidx, filename in enumerate(filenames):
        embedding_matrix, idx_to_cui, cui_to_idx = read_embedding_matrix_cui(filename)
        embedding_idx_cui[filename] = (embedding_matrix, idx_to_cui, cui_to_idx)
        if fileidx == 0:
            overlapping_cuis.update(set(cui_to_idx.keys()))
        else:
            overlapping_cuis.intersection_update(set(cui_to_idx.keys()))
    overlapping_cuis = list(overlapping_cuis)
    
    idx_of_overlapping_cuis = {}
    for filename in filenames:
        idx_of_overlapping_cuis[filename] = []

    idx_to_cui = {}
    cui_to_idx = {}
    for idx, cui in enumerate(overlapping_cuis):
        idx_to_cui[idx] = cui
        cui_to_idx[cui] = idx
        for filename in filenames:
            idx_of_overlapping_cuis[filename].append(embedding_idx_cui[filename][2][cui])
    filename_to_embedding_matrix = {}
    for filename in filenames:
        idx_of_overlapping_cuis[filename] = np.array(idx_of_overlapping_cuis[filename]).astype(int) #JJN added due to change in numpy
        filename_to_embedding_matrix[filename] = embedding_idx_cui[filename][0][idx_of_overlapping_cuis[filename]]
    return filename_to_embedding_matrix, idx_to_cui, cui_to_idx

# JJN changed to support csv embedding files
def read_embedding_matrix_cui(filename):
    concept_to_cui, cui_to_concept = get_cui_concept_mappings() # comment out this after fix input
    
    # Assume embedding files with .txt are deliminated with ' ' and ',' if .csv
    if filename.endswith('.txt'):
        delim = ' '
    elif filename.endswith('.csv'):
        delim = ','
    else:
        raise Exception('embedding file must be .txt or .csv depending on deliminator')
    

    with open(str(data_folder / filename), 'r') as infile:            
        first_line = infile.readline().strip().split(delim)
        embedding_num = int(first_line[0])
        dimension = int(first_line[1])
        #embedding_num, dimension = map(int, infile.readline().strip().split(delim))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines() 
        idx_to_cui = {}
        cui_to_idx = {}
        for idx in xrange(embedding_num-1):
            datum = data[idx+1].strip().split(delim)
            cui = datum[0]
            if cui[0] != 'C':
                if cui in concept_to_cui:
                    cui = concept_to_cui[cui].strip() # JJN added to remove \n
            embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
            idx_to_cui[idx] = cui
            cui_to_idx[cui] = idx
        return embedding_matrix, idx_to_cui, cui_to_idx