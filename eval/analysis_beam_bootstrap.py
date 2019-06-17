# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 20:30:55 2019

@author: jjnun
"""
from __future__ import division
from sklearn.metrics.pairwise import cosine_similarity
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_cui
from cui_icd9_helpers import get_cui_to_systems, get_icd9_cui_mappings_rangeok, get_cui_may_treat_prevent_icd9, get_cui_to_icd9_drug_or_diag
import random

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")


def write_results_to_file(results,f):
    for row in results:
        f.write(','.join(row) + '\n')

def get_beam_bootstrap_by_systems(filenames_type, icd9_systems, cui_to_icd9_dicts, results):
    """JJN: Using Beam et al's boostrap method, calculates the percentage of known drug-disease relations for a given system
    are within the top 5% of relations in a null distribution
    """
    filename_to_embedding_matrix, idx_to_cui, cui_to_idx, cui_to_icd9_types = generate_overlapping_sets_cui(filenames_type, True, cui_to_icd9_dicts)
    
    # Obtain dictionary between cuis and the ICD9 disease systems they're part of/related to
    cui_to_systems = get_cui_to_systems(cui_to_icd9_types, icd9_systems)

    
    # And also build a similiar idx_to_system, and a list of idx that are drugs, or diags
    idx_to_tr_pr_systems = {}
    idx_to_tr_pr_cuis = {}
    drug_idxs = [] # idx's of cuis repersenting a drug
    diag_idxs = [] # '                        ' a diagnosis
    for idx in idx_to_cui:
        cui = idx_to_cui[idx]
        systems = cui_to_systems[cui]
        idx_to_tr_pr_systems[idx] = systems
        if cui_to_icd9_types[cui]['icd9_type'] == 'diag':
            diag_idxs.append(idx)
        elif cui_to_icd9_types[cui]['icd9_type'] == 'drug':
            drug_idxs.append(idx)
            tr_pr_cuis = cui_to_icd9_types[cui]['cuis']
            idx_to_tr_pr_cuis[idx] = tr_pr_cuis
            assert len(tr_pr_cuis) == len(systems), 'While building dictionary, cuis and systems had different length for cui: ' + cui
        else:
            raise Exception('Each cui used must repersent a ICD9 diagnosis, or drug that treats one')
    
    # Get list of ICD9 system names in this analysis
    icd9_systems_names = []
    for row in icd9_systems: icd9_systems_names.append(row[0])    
    
    filename_index = 0
    for filename, embedding_type, _ in filenames_type:
#        # Matrix to convert cui to positions in the relevant filename
         embedding_matrix = filename_to_embedding_matrix[filename]
         
         # Create a null distribution by a bootstrap sample involving a set number 
         # of cosine similiarities between random drug and disease pairs 
         n_bootstrap = 10000
         null_cos_sims = []
         p_value = 0.05
         for j in range(n_bootstrap):
             rand_drug = random.choice(drug_idxs)
             rand_diag = random.choice(diag_idxs)
             vec_1 = embedding_matrix[rand_drug,:]
             vec_2 = embedding_matrix[rand_diag,:]
             cos_sim = cosine_similarity([vec_1], [vec_2])[0,0]
             null_cos_sims.append(cos_sim)
        
         # Calculate the threshold for p < 0.05 significance in cosine similarity
         pos_threshold = int(n_bootstrap*p_value)
         sig_threshold = sorted(null_cos_sims)[-1*pos_threshold]        
         print 'Done with bootstrap. Have this many examples: ' + str(len(null_cos_sims))
         print 'Significance threshold is: ' + str(sig_threshold)
         
         systems_n = {}
         systems_sig = {}

         # Set up
         system_index = 0
         for system in icd9_systems_names:
             systems_n[system] = 0.0000001
             systems_sig[system] = 0
        # Test all drug-relations that have cuis in this system
         for idx in drug_idxs:
             cui_drug = idx_to_cui[idx]
             tr_pr_systems = idx_to_tr_pr_systems[idx]
             tr_pr_cuis    = idx_to_tr_pr_cuis[idx]
             assert len(tr_pr_systems) == len(tr_pr_cuis), 'Length must be same as they correspond'
             
             for i in range(len(tr_pr_systems)):
                 tr_pr_system = tr_pr_systems[i]
                 tr_pr_cui = tr_pr_cuis[i]
                 
                 # Ignore this treated or prevented disease if it doesn't have a vector repersentation. 
                 # Not checked until this point as for other analysis we just want to know the system a drug treats
                 # And our cui_to_icd9 dictionary has more entries than embedings
                 if tr_pr_cui in cui_to_idx.keys(): 
                     systems_n[tr_pr_system] += 1    
                     vec_1 = embedding_matrix[cui_to_idx[cui_drug],:]
                     vec_2 = embedding_matrix[cui_to_idx[tr_pr_cui],:]
                     cos_sim = cosine_similarity([vec_1], [vec_2])[0,0]
                     
                     if cos_sim > sig_threshold: systems_sig[tr_pr_system] += 1
             
         # Store results
         system_index = 0
         for system in icd9_systems_names:
            results[system_index + 1][0] = re.sub(",", " ", system)
            results[system_index + 1][filename_index + 1] = '%2.5f' %(systems_sig[system]/systems_n[system]) ##, np.std(np.array(systems_dcg[system])))
            results[system_index + 1][-1] = str(systems_n[system]) # Number of examples used for this calculation. Will be re-written by each file but that's okay as always same
            system_index += 1
         filename_index += 1
        
    return results




def print_beam_bootstrap(filenames):
    """JJN: Prints and writes the result from using Beam et al's boostrap method to evaluate
    whether known drug-disease relations are in the top 5% against a null distribution
    """
    
    # Cui_to_icd9 mappings will be used
    cui_to_icd9 = get_icd9_cui_mappings_rangeok()
    # Create dictionaries linking drug cuis to the icd9 conditions they prevent or treat
    cui_icd9_tr, cui_icd9_pr = get_cui_may_treat_prevent_icd9(cui_to_icd9)
    cui_to_icd9_drug_or_diag = get_cui_to_icd9_drug_or_diag(cui_to_icd9, cui_icd9_tr, cui_icd9_pr)
    # Store in a dict to pass
    cui_to_icd9_dicts = {}
    cui_to_icd9_dicts['cui_to_icd9'] = cui_to_icd9
    cui_to_icd9_dicts['cui_to_icd9_drug_or_diag'] = cui_to_icd9_drug_or_diag
    
    
    # Text file containing the system, start, end. Note that 'end' is an integer, so will end up to next integer
    icd9_systems_file = 'icd9_systems.txt'
    # Parse above file to get the system names, starts, ends
    icd9_systems = []
    with open(str(data_folder / icd9_systems_file), 'r') as infile:
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
    
    
    print 'Beam Boostrap Method Examining Drug-Disease Relations'
    results = get_beam_bootstrap_by_systems(filenames, icd9_systems, cui_to_icd9_dicts, empty_results)
    for line in results: print line
    
    beam_bootstrap_by_system = 'beam_bootstrap_by_system_beamonly.csv'
    o = open(str(results_folder / beam_bootstrap_by_system ), 'w')
    write_results_to_file(results,o)
    o.close()    
