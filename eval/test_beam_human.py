# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 09:43:18 2019

@author: jjnun
"""
from __future__ import division
## from scipy.spatial.distance import cosine
from scipy.stats import spearmanr
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_cui
from cui_icd9_helpers import cui_in_system, get_icd9_cui_mappings_rangeok, get_cui_may_treat_prevent_icd9
from sklearn.metrics.pairwise import cosine_similarity


tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")

filenames_type = [['DeVine_etal_200.txt','cui','DeVine200']]

filename_to_embedding_matrix, idx_to_cui, cui_to_idx = generate_overlapping_sets_cui(filenames_type)
    
print 'Number of overlapping cuis between embeddings: ' + str(len(idx_to_cui))
    
#umnsrs_filename = 'UMNSRS_relatedness_mod458_word2vec.csv'
umnsrs_filename = 'UMNSRS_similarity_mod449_word2vec.csv'
with open(str(data_folder / umnsrs_filename), 'rU') as f:
    umnsrs_rows = f.readlines()[1:]

# Find all cuis in the overlapping matrix
cuis = cui_to_idx.keys()

# Start arrays to store filenames, and to store the comparison values
filename_all = []
value_all = []

## TODO REMOVE THIS JUST FOR DEBUGGING
o = open(str(results_folder / 'overlappingcuis.txt' ), 'w')
for cui in cuis: o.write(str(cui)+'\n')
o.close()
## -----------------------------------

# Use a dict to keep track of how many of the UMNSRS comparisons are used in this system
# To avoid extra code, this will be reset with each file, but will all be the same
compares = {}
compares['possible'] = 0 # Total # of UMNSRS that could be used, as have both cuis found in the overlapping matrix for the embeddings
##compares['total'] = 0 # Total # of UMNSRS comparisons used in this system
##compares['diags'] = 0 # Total # with a diagnosis from this system
##compares['drugs'] = 0 # Total # with a drug from this system
files_processed = 0 # Only set above for first file, again, same for all

missing_cuis = []
  
for filename, embedding_type, _ in filenames_type:
    #Contains the unmsrs scores from comparisons that have to do with this system
    unmsrs_scores = []
    #Contains the vector cosine similarity between the embedding vectors
    veccos_scores = []
    # Matrix to convert cui to positions in the relevant filename
    embedding_matrix = filename_to_embedding_matrix[filename]

    for row_str in umnsrs_rows:
        row = row_str.strip().split(',')
        umnsrs_rating = row[0]
        cui_1 = row[4]
        cui_2 = row[5]
        
        if (cui_1 in cuis) and (cui_2 in cuis):
            if files_processed == 0: compares['possible'] += 1
            ## in_system_diag_1, in_system_drug_1 = cui_in_system(cui_1, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
            ## in_system_diag_2, in_system_drug_2 = cui_in_system(cui_2, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
            if 1==1:
                # Calculate Cosine similiarity between vectors
                vec_1 = embedding_matrix[cui_to_idx[cui_1],:]
                vec_2 = embedding_matrix[cui_to_idx[cui_2],:]
                cos_sim = cosine_similarity([vec_1], [vec_2])[0,0]
                
                # Recrod scores
                veccos_scores.append(cos_sim)
                unmsrs_scores.append(float(umnsrs_rating))
                
# =============================================================================
                ## Debugging/
                if cui_1 == 'C0000970' and cui_2 == 'C0020740':
                    print 'This is cui_1: ' + cui_1
                    print 'This is cui_2: ' + cui_2
                    print 'This is the cos_sim: ' + str(cos_sim)
                    print 'This is umnsrs rating: ' + str(umnsrs_rating)
                    # print 'Part of 1st vec' + str(vec_1[0:5])
                    # print 'Part of 2nd vec' + str(vec_2[0:5])
                    # print 'This is the cui_to_idx ' + str(cui_to_idx[cui_1])
                    # print 'This is start of vector above: ' + str(embedding_matrix[cui_to_idx[cui_1],0:10])
                     
                     ## /Debugging----
# =============================================================================
                
                
                # Keep track of # of UMNSRS comparisons for this system
# =============================================================================
#                 if files_processed == 0:
#                     in_system_diag = in_system_diag_1 or in_system_diag_2
#                     in_system_drug = in_system_drug_1 or in_system_drug_2
#                     
#                     if in_system_diag: 
#                         compares['diags'] += 1
#                         compares['total'] += 1
#                     if in_system_drug:
#                         compares['drugs'] += 1
#                         compares['total'] += 1
#                     if in_system_drug and in_system_diag:
#                         compares['total'] -= 1 # Don't double count
# =============================================================================
                
        elif cui_1 in cuis:
            missing_cuis.append(cui_2)
        else:
            missing_cuis.append(cui_1)
    ## print(set(missing_cuis))
    print unmsrs_scores[0:5]
    print veccos_scores[0:5]
        
    rho, pval = spearmanr(unmsrs_scores,veccos_scores)
    print 'Rho for all is: ' + str(rho)
    print 'Based on this many comparisons: ' + str(compares['possible'])
    filename_all.append((filename))
    value_all.append(rho)
    files_processed += 1
    
    