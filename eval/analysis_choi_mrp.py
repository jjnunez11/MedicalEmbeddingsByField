# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:04:31 2019

@author: jjnun
"""
from __future__ import division
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import cosine_similarity
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_icd9
from cui_icd9_helpers import get_coarse_icd9_pairs, get_icd9_pairs, get_icd9_to_description
from alt_cosdist import cosine_vectorized_v2


tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")

#def mem_sav_cdist_rank(array, num_of_neighbor):
#    """JJN: When embedding set is large, original code will not work as call to cdist
#    creates too large a matrix (eg 100,000^2). This function is slower but will save memory.
#    Calculates cosine distance one row at a time, and only saves the top num_of_neighbor + 1
#    entries per row, so the rank matrix will only contain embeddings*(num_of_neighbor + 1) entries
#    """
#    
#    r, c = array.shape
#    
#    ranks = np.zeros((r,num_of_neighbor + 1))
#
#    for i in range(r):
#        row_cdist = cdist(array[i].reshape((1,c)),array,'cosine')[0]
#        row_ranks = np.argsort(row_cdist)[0:num_of_neighbor + 1]
#        ranks[i] = row_ranks
#        if i % 1000 == 0: print i
#    return ranks


def get_choi_mrp_by_system(filenames_type, num_of_neighbor, start, end, type='f'):
    """ JJN: Calculates Choi et al's Medical Relatedness Property by ICD9 system.
    Heavily based upon similiar code from this author's Github 
    """
    
    filename_to_embedding_matrix, idx_to_icd9, icd9_to_idx = generate_overlapping_sets_icd9(filenames_type)

    if type == 'c':
        icd9_pairs = get_coarse_icd9_pairs(set(icd9_to_idx.keys()))
    else:
        icd9_pairs = get_icd9_pairs(set(icd9_to_idx.keys()))

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
    raw_all = [] # Store the raw scores for each filename
    for filename, embedding_type, _ in filenames_type:
        icd9_embeddings = filename_to_embedding_matrix[filename]
        r,c = icd9_embeddings.shape
        
        # Following replaced as won't work with large embeddings. Is faster though
        Y = cdist(icd9_embeddings, icd9_embeddings, 'cosine')        
        ranks = np.argsort(Y) 
        # End of what was replace
        cumulative_ndcgs = []

        for icd9 in icd9_in_system:
            # Changed for neighbours to be found "as needed" to memory cost, though is slower
            ##a_idx = icd9_to_idx[icd9]
            ##a_cdist = cdist(icd9_embeddings[a_idx].reshape((1,c)),icd9_embeddings,'cosine')[0]
            ##target = np.argsort(a_cdist)[1:num_of_neighbor + 1]
            # end of change to allow larger embeddings to be used

            # Uncomment following and block before for loop to return to faster but memory-needing version
            target = ranks[icd9_to_idx[icd9], 1:num_of_neighbor+1]
            num_of_possible_hits = 0
            
            icd9_to_remove = set()

            for val in icd9_pairs[icd9]:
                if val not in icd9_to_idx:
                    icd9_to_remove.add(val)
            icd9_pairs[icd9].difference(icd9_to_remove)

            num_of_possible_hits = min(len(icd9_pairs[icd9]), num_of_neighbor)

            dcg = 0
            best_dcg = np.sum(np.reciprocal(np.log2(range(2, num_of_possible_hits+2))))
            for i in xrange(num_of_neighbor):
                if idx_to_icd9[target[i]] in icd9_pairs[icd9]:
                    dcg += np.reciprocal(np.log2(i+2))

            cumulative_ndcgs.append(dcg/best_dcg)
        filename_all.append((filename))
        raw_all.append(cumulative_ndcgs)
        value_all.append(np.mean(np.array(cumulative_ndcgs)))
    return filename_all, value_all, len(icd9_in_system), raw_all


def print_choi_mrp(filenames, num_of_nn=40):
    """ JJN: Prints and writes Choi's Medical Relatednes Property by ICD9 system
    """
    
    # csv file to write results to
    choi_mrp_by_system = 'choi_mrp_by_system_withoutfix.csv'
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
    with open(str(data_folder / icd9_systems_file), 'r') as infile:
        data = infile.readlines()
        for row in data:
            icd9_systems.append(row.strip().split('|'))
    
    
    ## TODO: DELETE THIS, ONLY FOR "ALL" SYSTEM
    ##icd9_systems_just_all = []
    ##icd9_systems_just_all.append(icd9_systems[0])
    ##icd9_systems = icd9_systems_just_all
    
    print 'Choi Medical Relatedness Property by ICD9 system'
    for system in icd9_systems:
        system_name = system[0]
        start = float(system[1])
        end = float(system[2])
        
        filename_to_print, ndcgs_to_print, comparisons_in_cuis, raw_all = get_choi_mrp_by_system(filenames, num_of_nn, start, end, 'c')
        # Write ncdgs to file
        ndcgs_rounded = [round(x*100,2) for x in ndcgs_to_print]
        ncdgs_str = ','.join(map(str, ndcgs_rounded))
        o.write('\n' + re.sub(",", " ", system_name) + ',') # Replace commas with space to use as csv
        o.write(ncdgs_str)
        o.write(", " + str(comparisons_in_cuis))
        # Print ncdfs 
        print '\n' + system_name
        for filename, ndcg in zip(filename_to_print, ndcgs_to_print):
            print '%s & %.2f \\\\' %(filename.split('/')[-1], ndcg*100)
        print "Number of examples: " + str(comparisons_in_cuis)
        
        
        # New addition: Print raw scores
        # csv file to write raw results to
        system_name_compact = re.sub(",", "", system_name)
        choi_mrp_raw_system = 'choi_mrp_raw_' + system_name_compact + '.csv'
        o_raw = open(str(results_folder / choi_mrp_raw_system ), 'w')
        
        # Print raw scores
        for filename, raw_scores in zip(filename_to_print, raw_all):
            o_raw.write(filename + ",")
            o_raw.write( ','.join(map(str, raw_scores)))
            o_raw.write('\n')
            
        o_raw.close()
            
    o.close()
    