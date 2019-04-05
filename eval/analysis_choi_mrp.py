# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:04:31 2019

@author: jjnun
"""
from __future__ import division
##import argparse
import numpy as np
##import scipy as sp
from scipy.spatial.distance import cdist
##from scipy.stats import spearmanr
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_icd9
from cui_icd9_helpers import get_coarse_icd9_pairs, get_icd9_pairs, get_icd9_to_description

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")


def get_choi_mrp_by_system(filenames_type, num_of_neighbor, start, end, type='f'):
    filename_to_embedding_matrix, idx_to_icd9, icd9_to_idx = generate_overlapping_sets_icd9(filenames_type)
    #print len(icd9_to_idx.keys())
    if type == 'c':
        icd9_pairs = get_coarse_icd9_pairs(set(icd9_to_idx.keys()))
    else:
        icd9_pairs = get_icd9_pairs(set(icd9_to_idx.keys()))

    #print len(icd9_pairs)
    icd9_to_check = set(icd9_pairs.keys())
    icd9_to_check.intersection_update(set(icd9_to_idx.keys()))

    print 'Overlapping icd9 number: ' + str(len(icd9_to_check))
    
    icd9_to_description = get_icd9_to_description()
    for icd9 in icd9_to_idx.keys():
        if icd9 not in icd9_to_description:
            if tree.find(icd9):
                icd9_to_description[icd9] = tree.find(icd9).description.encode('utf-8')
            else:
                icd9_to_description[icd9] = ''
    
    # Remove ICD9 codes that contain an alphanumeric charector, these are 'supplemental/misc', not applicable
    icd9_to_check_noV = [x for x in icd9_to_check if not(any(char.isalpha() for char in x))]
    
    # JJN: Select ICD9 codes pertaining to a specific system according to proivded start, end integer values 
    icd9_in_system = [ x for x in icd9_to_check_noV if start <= float(x) < end + 1]

    filename_all = [] 
    value_all = []
    for filename, embedding_type, _ in filenames_type:
        icd9_embeddings = filename_to_embedding_matrix[filename]
        Y = cdist(icd9_embeddings, icd9_embeddings, 'cosine')
        ranks = np.argsort(Y)
    
        cumulative_ndcgs = []

        for icd9 in icd9_in_system:
            target = ranks[icd9_to_idx[icd9], 1:num_of_neighbor+1]
            num_of_possible_hits = 0
            
            icd9_to_remove = set()

            for val in icd9_pairs[icd9]:
                if val not in icd9_to_idx:
                    icd9_to_remove.add(val)
            icd9_pairs[icd9].difference(icd9_to_remove)

            num_of_possible_hits = min(len(icd9_pairs[icd9]), num_of_neighbor)
            #print icd9 + '(' + str(num_of_possible_hits) + ')',
            #if icd9 in icd9_to_description:
            #    print '(' + icd9_to_description[icd9] + ')',
            #print ''
            #print '-------------------------------------------'
            dcg = 0
            best_dcg = np.sum(np.reciprocal(np.log2(range(2, num_of_possible_hits+2))))
            for i in xrange(num_of_neighbor):
                if idx_to_icd9[target[i]] in icd9_pairs[icd9]:
                    dcg += np.reciprocal(np.log2(i+2))
                    #print 'hit: ',
                #else:
                    #print '     ',
                #print idx_to_icd9[target[i]],
                #if idx_to_icd9[target[i]] in icd9_to_description:
                    #print icd9_to_description[idx_to_icd9[target[i]]],
                #print ''
            #print dcg/best_dcg
            #print ''
            cumulative_ndcgs.append(dcg/best_dcg)
        filename_all.append((filename))
        value_all.append(np.mean(np.array(cumulative_ndcgs)))
    return filename_all, value_all, len(icd9_in_system)




# JJN: Prints the Medical Relatedness Property by ICD9 system
def print_choi_mrp(filenames, num_of_nn=40):
    # csv file to write results to
    choi_mrp_by_system = 'choi_mrp_by_system.csv'
    o = open(str(results_folder / choi_mrp_by_system ), 'w')
    o.write('ICD9 System,')
    # Write headers from 3rd entry in orig_files_all.txt
    o.write(",".join(list(map(lambda x: x[2], filenames))))
    # Write Examples per System
    o.write(", Examples in System")
    
    # Text file containing the system, start, end. Note that 'end' is an integer, so will end up to next integer
    icd9_systems_file = 'icd9_systems.txt'
    # Parse above file to get the system names, starts, ends
    icd9_systems = []
    with open(icd9_systems_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            icd9_systems.append(row.strip().split('|'))
    
    print 'Choi Medical Relatedness Property by ICD9 system'
    for system in icd9_systems:
        system_name = system[0]
        start = float(system[1])
        end = float(system[2])
        
        filename_to_print, ndcgs_to_print, comparisons_in_cuis = get_choi_mrp_by_system(filenames, num_of_nn, start, end, 'c')
        # Write ncdgs to file
        ndcgs_rounded = [round(x*100,2) for x in ndcgs_to_print]
        ncdgs_str = ','.join(map(str, ndcgs_rounded))
        o.write('\n' + re.sub(",", " ", system_name) + ',') # Replace commas with space to use as csv
        o.write(ncdgs_str)
        o.write(", " + str(comparisons_in_cuis))
        # Print ncdfs 
        print '\n' + system_name
        for file_name, ndcg in zip(filename_to_print, ndcgs_to_print):
            print '%s & %.2f \\\\' %(file_name.split('/')[-1], ndcg*100)
        print "Number of examples: " + str(comparisons_in_cuis)