# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 20:30:55 2019

@author: jjnun
"""
from __future__ import division
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import cdist
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_cui
from cui_icd9_helpers import get_coarse_icd9_pairs, get_icd9_pairs, get_icd9_to_description, get_cui_to_systems
from cui_icd9_helpers import cui_in_system, get_icd9_cui_mappings_rangeok, get_cui_may_treat_prevent_icd9, get_cui_to_icd9_drug_or_diag
import random

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")


def write_results_to_file(results,f):
    for row in results:
        f.write(','.join(row) + '\n')


def get_new_centroid_by_system(filenames_type, icd9_systems, cui_to_icd9_dicts, results):
    filename_to_embedding_matrix, idx_to_cui, cui_to_idx, cui_to_icd9_types = generate_overlapping_sets_cui(filenames_type, True, cui_to_icd9_dicts)
    
    # Obtain dictionary between cuis and the ICD9 disease systems they're part of/related to
    cui_to_systems = get_cui_to_systems(cui_to_icd9_types, icd9_systems)
    
    # Obtain dictionary between drug cui and the cuis they treat or prevent
    #cui_to_tr_pr_cuis = {}
    #for cui in cui_to_icd9_types.keys():
    #    if cui_to_icd9_types[cui]['icd9_type'] == 'drug':
    #        tr_pr_cuis = cui_to_icd9_types[cui]['cuis']
    #        cui_to_tr_pr_cuis[cui] = tr_pr_cuis
    
    

    
#    # And also build a similiar idx_to_system, and a list of idx that are drugs, or diags
#    idx_to_tr_pr_systems = {}
#    idx_to_tr_pr_cuis = {}
#    drug_idxs = [] # idx's of cuis repersenting a drug
#    diag_idxs = [] # '                        ' a diagnosis
#    for idx in idx_to_cui:
#        cui = idx_to_cui[idx]
#        systems = cui_to_systems[cui]
#        idx_to_tr_pr_systems[idx] = systems
#        if cui_to_icd9_types[cui]['icd9_type'] == 'diag':
#            diag_idxs.append(idx)
#        elif cui_to_icd9_types[cui]['icd9_type'] == 'drug':
#            drug_idxs.append(idx)
#            tr_pr_cuis = cui_to_icd9_types[cui]['cuis']
#            idx_to_tr_pr_cuis[idx] = tr_pr_cuis
#            assert len(tr_pr_cuis) == len(systems), 'While building dictionary, cuis and systems had different length for cui: ' + cui
#        else:
#            raise Exception('Each cui used must repersent a ICD9 diagnosis, or drug that treats one')
#    
    # Get list of ICD9 system names in this analysis
    icd9_systems_names = []
    
    for row in icd9_systems: icd9_systems_names.append(row[0])    
    n_of_systems = len(icd9_systems_names)
    
    
    
        
    
    
    #test_cuis = cui_to_systems.keys()[0:5]
    #for cui in test_cuis:
    #    print 'Here is the type of some cuis: ' + cui_to_icd9_types[cui]['icd9_type']
    #print 'Number of idx that are diags: ' + str(len(diag_idxs))
    #print 'Number of idx that are drugs: ' + str(len(drug_idxs))
    
    ##print "Here are some keys" + str(cui_to_systems.keys()[0:5])
    ##print "Here are some values" + str(cui_to_systems.values()[0:5])
    ##print "Here's one example of a cui to system: " + str(cui_to_systems['C0206180'])
    
    filename_index = 0
    for filename, embedding_type, _ in filenames_type:
#        # Matrix to convert cui to positions in the relevant filename
         embedding_matrix = filename_to_embedding_matrix[filename]
         
         # Make centroids
         systems_centroid = {}
         systems_n = {}
         systems_correct = {}
         for system in icd9_systems_names:
             systems_centroid[system] = []
             systems_n[system] = 0
             systems_correct[system] = 0

         for cui in cui_to_idx.keys():
             cui_vec = embedding_matrix[cui_to_idx[cui],:]
             for system in cui_to_systems[cui]:
                 systems_centroid[system].append(cui_vec)
                 ##systems_n[system] += 1
        
         for system in icd9_systems_names:
             systems_centroid[system] = np.array(systems_centroid[system])
             systems_centroid[system] = np.mean(systems_centroid[system], axis=0)
         
         # Calculate accuracies using centroids. 
         for cui in cui_to_idx.keys():
             cui_vec = embedding_matrix[cui_to_idx[cui],:]
             true_systems = cui_to_systems[cui]
             
             cos_sims = np.zeros(n_of_systems)
             
             # Generate list of cos similarities with the system centroids
             for i in range(n_of_systems):
                 system = systems_centroid.keys()[i]
                 system_vec = systems_centroid[system]
                 cos_sim = cosine_similarity([cui_vec], [system_vec])[0,0]
                 cos_sims[i] = cos_sim
             
             n = len(true_systems) # Number of systems this cui treats or prevents or 1 if diagnosis
             pred_systems = [icd9_systems_names[i] for i in np.argsort(cos_sims)[-n:]]
             
             for system in true_systems:
                 systems_n[system] += 1
                 if system in pred_systems: systems_correct[system] += 1

             #TODO: Iterate through systems, make a new dict for cos_sims
             #Then, find a way to find top n systems
             # Lastly, go through the systems, +1 for n, if correct, +1 for correct
             
             
         
         
         
         
# =============================================================================
#          
#          
#          # Create a null distribution by a bootstrap sample involving a set number 
#          # of cosine similiarities between random drug and disease pairs    
#          n_bootstrap = 1000
#          null_cos_sims = []
#          p_value = 0.05
#          for j in range(n_bootstrap):
#              rand_drug = random.choice(drug_idxs)
#              rand_diag = random.choice(diag_idxs)
#              vec_1 = embedding_matrix[rand_drug,:]
#              vec_2 = embedding_matrix[rand_diag,:]
#              cos_sim = cosine_similarity([vec_1], [vec_2])[0,0]
#              null_cos_sims.append(cos_sim)
#         
#          # Calculate the threshold for p < 0.05 significance in cosine similarity
#          pos_threshold = int(n_bootstrap*p_value)
#          sig_threshold = sorted(null_cos_sims)[-1*pos_threshold]        
#          print 'Done with bootstrap. Have this many examples: ' + str(len(null_cos_sims))
#          print 'Significance threshold is: ' + str(sig_threshold)
#          #print null_cos_sims[77:107]
#             
# #        Y = cdist(embedding_matrix, embedding_matrix, 'cosine')
# #        ranks = np.argsort(Y)
# #        systems_idx_dcg_err = {}
#          systems_n = {}
#          systems_sig = {}
# #        systems_n = {} # Count of examples found for this system
# #        print 'done calcualting distance'
# #        
#          
#         # Test all drug-relations that have cuis in this system
#          for idx in drug_idxs:
#              cui_drug = idx_to_cui[idx]
#              tr_pr_systems = idx_to_tr_pr_systems[idx]
#              tr_pr_cuis    = idx_to_tr_pr_cuis[idx]
#              assert len(tr_pr_systems) == len(tr_pr_cuis), 'Length must be same as they correspond'
#              
#              #assert len(tr_pr_systems) == len(tr_pr_cuis), 'The systems and cuis a drug treats or prevents must be the same length (so refer to same list)'    
#              for i in range(len(tr_pr_systems)):
#                  tr_pr_system = tr_pr_systems[i]
#                  tr_pr_cui = tr_pr_cuis[i]
#                  
#                  # Ignore this treated or prevented disease if it doesn't have a vector repersentation. 
#                  # Not checked until this point as for other analysis we just want to know the system a drug treats
#                  # And our cui_to_icd9 dictionary has more entries than embedings
#                  if tr_pr_cui in cui_to_idx.keys(): 
#                      systems_n[tr_pr_system] += 1    
#                      vec_1 = embedding_matrix[cui_to_idx[cui_drug],:]
#                      vec_2 = embedding_matrix[cui_to_idx[tr_pr_cui],:]
#                      cos_sim = cosine_similarity([vec_1], [vec_2])[0,0]
#                      
#                      if cos_sim > sig_threshold: systems_sig[tr_pr_system] += 1
#              
#              
#              
#              
# #            target = ranks[idx, 1:num_of_neighbor+1]
#                  # 
# #                for i in xrange(num_of_neighbor):
# #                    neighbor_idx = target[i]
# #                    neighbor_systems = idx_to_tr_pr_systems[neighbor_idx]
# #                    if system in neighbor_systems:
# #                        dcg += np.reciprocal(np.log2(i+2))
# #                        if err == 0:
# #                            err = 1/(1+i)
# #                systems_idx_dcg_err[system].append((idx, dcg, err))
# #                systems_dcg[system].append(dcg)
# #                systems_err[system].append(err)
# =============================================================================
#        
#        ##Temporary
         system_index = 0
         for system in icd9_systems_names:
            results[system_index + 1][0] = re.sub(",", " ", system)
            results[system_index + 1][filename_index + 1] = '%2.2f'  %(100*systems_correct[system]/systems_n[system]) ##, np.std(np.array(systems_dcg[system])))
            results[system_index + 1][-1] = str(systems_n[system]) # Number of examples used for this calculation. Will be re-written by each file but that's okay as always same
#            ##print '%50s (DCG) %2.5f %2.5f' %(system, np.mean(np.array(systems_dcg[system])), np.std(np.array(systems_dcg[system])))
#            ##print '%50s (ERR) %2.5f %2.5f' %(system, np.mean(np.array(systems_err[system])), np.std(np.array(systems_err[system])))
#            ##f.write('%50s (DCG) %2.5f %2.5f\n' %(type, np.mean(np.array(type_dcg[type])), np.std(np.array(type_dcg[type]))))
#            ##f.write('%50s (ERR) %2.5f %2.5f\n' %(type, np.mean(np.array(type_err[type])), np.std(np.array(type_err[type]))))
            system_index += 1
         filename_index += 1
        
    return results



# JJN: Prints the Medical Relatedness Property by ICD9 system
def print_new_centroid(filenames):
    # Cui_to_icd9 mappings will be used
    cui_to_icd9 = get_icd9_cui_mappings_rangeok()
    # Create dictionaries linking drug cuis to the icd9 conditions they prevent or treat
    cui_icd9_tr, cui_icd9_pr = get_cui_may_treat_prevent_icd9(cui_to_icd9)
    cui_to_icd9_drug_or_diag = get_cui_to_icd9_drug_or_diag(cui_to_icd9, cui_icd9_tr, cui_icd9_pr)
    # Store in a dict to pass
    cui_to_icd9_dicts = {}
    cui_to_icd9_dicts['cui_to_icd9'] = cui_to_icd9
    ##cui_to_icd9_dicts['cui_to_icd9_may_treat'] = cui_icd9_tr
    ##cui_to_icd9_dicts['cui_to_icd9_may_prevent'] = cui_icd9_pr
    cui_to_icd9_dicts['cui_to_icd9_drug_or_diag'] = cui_to_icd9_drug_or_diag
    
    
    # Text file containing the system, start, end. Note that 'end' is an integer, so will end up to next integer
    icd9_systems_file = 'icd9_systems.txt'
    # Parse above file to get the system names, starts, ends
    icd9_systems = []
    with open(icd9_systems_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            row_str = row.strip().split('|')
            if row_str[0] != 'all':
                icd9_systems.append([row_str[0], float(row_str[1]), float(row_str[2])])        
    
    # Make a numpy matrix to store results for easy printing
    cols = len(filenames) + 2 # System, Each Filename's MCSP, Examples per System
    rows = len(icd9_systems) + 1 # One for each system plus the header
    empty_results = [['-' for x in range(cols)] for y in range(rows)] 
    empty_results[0][0] = 'ICD9 Systems'
    empty_results[0][1:len(filenames) + 1] = list(map(lambda x: x[2], filenames))
    empty_results[0][-1] = 'Examples in System'
    
    
    print 'New Centroid Method Quantifying Accuracy of a System Centroid'
    results = get_new_centroid_by_system(filenames, icd9_systems, cui_to_icd9_dicts, empty_results)
    for line in results: print line
    
    out_filename = 'new_centroid_by_system.csv'
    o = open(str(results_folder / out_filename ), 'w')
    write_results_to_file(results,o)
    o.close()    
