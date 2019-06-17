# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 20:30:55 2019

@author: jjnun
"""
from __future__ import division
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize
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


def get_new_sysvec_by_system(filenames_type, icd9_systems, cui_to_icd9_dicts, results):
    """ For each medical system, calculates a repersentative vector based on the mean of normalized vecotrs repersenting medical conditions 
    in the system. Then caluclates the percentage of drugs known to treat or prevent a disease in this system who are closest to the 
    correct system's vector. For drugs that treat/prevent conditions in n>1 systems, scores as a correct prediction if
    correct system vector is one of closest n system vectors. 
    """
    
    filename_to_embedding_matrix, idx_to_cui, cui_to_idx, cui_to_icd9_types = generate_overlapping_sets_cui(filenames_type, True, cui_to_icd9_dicts)
    
    # Obtain dictionary between cuis and the ICD9 disease systems they're part of/related to
    cui_to_systems = get_cui_to_systems(cui_to_icd9_types, icd9_systems)

    # Get list of ICD9 system names in this analysis
    icd9_systems_names = []
    
    for row in icd9_systems: icd9_systems_names.append(row[0])    
    n_of_systems = len(icd9_systems_names)
    
    filename_index = 0
    for filename, embedding_type, _ in filenames_type:
#        # Matrix to convert cui to positions in the relevant filename
         embedding_matrix = filename_to_embedding_matrix[filename]
         
         # Make System Vectors
         systems_sysvec = {}
         systems_n = {}
         systems_correct = {}
         for system in icd9_systems_names:
             systems_sysvec[system] = []
             systems_n[system] = 0.0001
             systems_correct[system] = 0

         for cui in cui_to_idx.keys():
             if cui_to_icd9_types[cui]['icd9_type'] == 'diag':
                 cui_vec = embedding_matrix[cui_to_idx[cui],:]
                 cui_vec = normalize(cui_vec.reshape(1, -1))[0]
                 for system in cui_to_systems[cui]:
                     systems_sysvec[system].append(cui_vec)
                     # Below for generating random set for negative control
                     ##rand_system = random.choice(icd9_systems_names)
                     ##systems_sysvec[rand_system].append(cui_vec)
        
         for system in icd9_systems_names:
             systems_sysvec[system] = np.array(systems_sysvec[system])
             systems_sysvec[system] = np.mean(systems_sysvec[system], axis=0)
         
         # Calculate accuracies using System Vectors. 
         for cui in cui_to_idx.keys():
             if  cui_to_icd9_types[cui]['icd9_type'] == 'drug':
                 cui_vec = embedding_matrix[cui_to_idx[cui],:]
                 #cui_vec = cui_vec/np.linalg.norm(cui_vec) #Normalize
                 cui_vec = normalize(cui_vec.reshape(1, -1))[0]
                 true_systems = cui_to_systems[cui]
                 #print true_systems
                 cos_sims = np.zeros(n_of_systems)
             
                # Generate list of cos similarities with the system vectors
                 for i in range(n_of_systems):
                    #system = systems_sysvec.keys()[i]
                    system = icd9_systems_names[i] 
                    system_vec = systems_sysvec[system]
                    cos_sim = cosine_similarity([cui_vec], [system_vec])[0,0]
                    cos_sims[i] = cos_sim
             
                 n = len(true_systems) # Number of systems this cui treats or prevents or 1 if diagnosis
                 pred_systems = [icd9_systems_names[i] for i in np.argsort(cos_sims)[-n:]]
                 assert len(pred_systems) == n, 'wrong number of systems predicted: ' + str(n) + 'vs ' + str(len(pred_systems))
                
                
                 for system in true_systems:
                    systems_n[system] += 1
                    #print 'Here is the system: ' + str(system)
                    #print 'Here is the pred_systems: ' + str(pred_systems)
                    #print 'Here is the true_systems: ' + str(true_systems)
                    #if system in pred_systems: print 'pred correct!'
                    if system in pred_systems: systems_correct[system] += 1
                    #rand_system = random.choice(icd9_systems_names)
                    #if rand_system in pred_systems: systems_correct[system] += 1
                 
         system_index = 0
         for system in icd9_systems_names:
            results[system_index + 1][0] = re.sub(",", " ", system)
            results[system_index + 1][filename_index + 1] = '%2.2f'  %(100*systems_correct[system]/systems_n[system]) ##, np.std(np.array(systems_dcg[system])))
            results[system_index + 1][-1] = str(int(systems_n[system])) # Number of examples used for this calculation. Will be re-written by each file but that's okay as always same
            system_index += 1
         filename_index += 1
        
    return results

def print_new_sysvec(filenames):
    """ Prints and writes the results for each medical system from our new method
    based on caluclating a medical systen's mean vector. 
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
    
    print 'New System Vector Method Quantifying Accuracy of a Systems Mean Vector to Predict a Drugs Related System'
    results = get_new_sysvec_by_system(filenames, icd9_systems, cui_to_icd9_dicts, empty_results)
    for line in results: print line
    
    out_filename = 'new_sysvec_by_system_diag_to_drug_beamonly.csv'
    o = open(str(results_folder / out_filename ), 'w')
    write_results_to_file(results,o)
    o.close()    
