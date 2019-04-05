# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:06:34 2019

@author: jjnun
"""
from __future__ import division
from scipy.spatial.distance import cosine
from scipy.stats import spearmanr
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_cui
from cui_icd9_helpers import cui_in_system, get_icd9_cui_mappings_rangeok, get_cui_may_treat_prevent_icd9

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")



# JJN Calculates the Spearman Correlation Coefficient between UMNSRS ratings and vector cosines when a pair contains
    # Either a diagnosis in the ICD 9 category, or treats or prevents a condition in that ICD 9 category
def get_yu_umnsrs_cor_by_system(filenames_type, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9):
    #filename_to_embedding_matrix, idx_to_icd9, icd9_to_idx = generate_overlapping_sets_cui(filenames_type)
    filename_to_embedding_matrix, idx_to_cui, cui_to_idx = generate_overlapping_sets_cui(filenames_type)
    
    print 'Number of overlapping cuis between embeddings: ' + str(len(idx_to_cui))
        
    umnsrs_filename = 'UMNSRS_relatedness_mod458_word2vec.csv'
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
    compares['total'] = 0 # Total # of UMNSRS comparisons used in this system
    compares['diags'] = 0 # Total # with a diagnosis from this system
    compares['drugs'] = 0 # Total # with a drug from this system
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
                in_system_diag_1, in_system_drug_1 = cui_in_system(cui_1, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
                in_system_diag_2, in_system_drug_2 = cui_in_system(cui_2, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
                if in_system_diag_1 or in_system_drug_1 or in_system_diag_2 or in_system_drug_2:
                    # Calculate Cosine similiarity between vectors
                    vec_1 = embedding_matrix[cui_to_idx[cui_1],:]
                    vec_2 = embedding_matrix[cui_to_idx[cui_2],:]
                    cos_sim = cosine(vec_1, vec_2)
                    
                    # Recrod scores
                    veccos_scores.append(cos_sim)
                    unmsrs_scores.append(umnsrs_rating)
                    
# =============================================================================
#                     ## Debugging/
#                     if cui_1 == 'C0019270':
#                         print 'This is cui_1: ' + cui_1
#                         print 'This is cui_2: ' + cui_2
#                         print 'This is the cos_sim: ' + str(cos_sim)
#                         print 'This is umnsrs rating: ' + str(umnsrs_rating)
#                         print 'Part of 1st vec' + str(vec_1[0:5])
#                         print 'Part of 2nd vec' + str(vec_2[0:5])
#                     # print 'This is the cui_to_idx ' + str(cui_to_idx[cui_1])
#                     # print 'This is start of vector above: ' + str(embedding_matrix[cui_to_idx[cui_1],0:10])
#                     
#                     ## /Debugging----
# =============================================================================
                    
                    
                    # Keep track of # of UMNSRS comparisons for this system
                    if files_processed == 0:
                        in_system_diag = in_system_diag_1 or in_system_diag_2
                        in_system_drug = in_system_drug_1 or in_system_drug_2
                        
                        if in_system_diag: 
                            compares['diags'] += 1
                            compares['total'] += 1
                        if in_system_drug:
                            compares['drugs'] += 1
                            compares['total'] += 1
                        if in_system_drug and in_system_diag:
                            compares['total'] -= 1 # Don't double count
                    
            elif cui_1 in cuis:
                missing_cuis.append(cui_2)
            else:
                missing_cuis.append(cui_1)
        ## print(set(missing_cuis))
        
            
        rho, pval = spearmanr(unmsrs_scores,veccos_scores)
        filename_all.append((filename))
        value_all.append(rho)
        files_processed += 1
    
    ## print 'Here is value_all: ' + str(value_all)        
    return filename_all, value_all, compares
    



# JJN: Prints the Spearman Correlation with Relevant Comparisons in the UMNSRS database by ICD9 system
def print_yu_umnsrs_cor(filenames):
    # Cui_to_icd9 mappings will be used
    cui_to_icd9 = get_icd9_cui_mappings_rangeok()
    # Create dictionaries linking drug cuis to the icd9 conditions they prevent or treat
    cui_icd9_tr, cui_icd9_pr = get_cui_may_treat_prevent_icd9(cui_to_icd9)
    
    # csv file to write results to
    yu_umnsrs_cor_by_system = 'yu_umnsrs_cor_by_system.csv'
    o = open(str(results_folder / yu_umnsrs_cor_by_system ), 'w')
    o.write('ICD9 System,')
    # Write headers from 3rd entry in orig_files_all.txt
    o.write(",".join(list(map(lambda x: x[2], filenames))))
    
    # Write headings for the # of UMNSRS comparisons used
    o.write(", UMNSRS Comparisons with cuis found")
    o.write(", Total UMNSRS comparisons for this system")
    o.write(", Diag comparisons")
    o.write(", Drug comparisons")
    
    # Text file containing the system, start, end. Note that 'end' is an integer, so will end up to next integer
    icd9_systems_file = 'icd9_systems.txt'
    # Parse above file to get the system names, starts, ends
    icd9_systems = []
    with open(icd9_systems_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            icd9_systems.append(row.strip().split('|'))
    
    print 'Yu Spearman Correlation with UMNSRS ratings by ICD9 system'
    for system in icd9_systems:
        system_name = system[0]
        start = float(system[1])
        end = float(system[2])
        
        filename_to_print, ndcgs_to_print, compares = get_yu_umnsrs_cor_by_system(filenames, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
        # Write ncdgs to file
        ndcgs_rounded = [round(x*100,2) for x in ndcgs_to_print]
        ncdgs_str = ','.join(map(str, ndcgs_rounded))
        o.write('\n' + re.sub(",", " ", system_name) + ',') # Replace commas with space to use as csv
        o.write(ncdgs_str)
        o.write(", " + str(compares['possible']))
        o.write(", " + str(compares['total']))
        o.write(", " + str(compares['diags']))
        o.write(", " + str(compares['drugs']))
        # Print ncdfs 
        print '\n' + system_name
        for file_name, ndcg in zip(filename_to_print, ndcgs_to_print):
            print '%s & %.2f \\\\' %(file_name.split('/')[-1], ndcg*100)
        print "Number of comparisons with both cuis present: " + str(compares['possible'])
        print "Number of comparisons involving this system: " + str(compares['total'])
        print "Number of comparisons with a drug from this system: " + str(compares['drugs'])
        print "Number of comparisons with a diag from this system: " + str(compares['diags'])
