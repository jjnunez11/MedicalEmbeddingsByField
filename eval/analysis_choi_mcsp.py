# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 11:57:27 2019

@author: jjnun
"""

from __future__ import division
import numpy as np
from scipy.spatial.distance import cdist
from icd9 import ICD9
from pathlib import Path
import re
from embed_helpers import generate_overlapping_sets_cui
from cui_icd9_helpers import get_cui_to_systems, get_icd9_cui_mappings_rangeok, get_cui_may_treat_prevent_icd9, get_cui_to_icd9_drug_or_diag

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")

def write_results_to_file(results,f):
    for row in results:
        f.write(','.join(row) + '\n')

def get_choi_mcsp_by_systems(filenames_type, num_of_neighbor, icd9_systems, cui_to_icd9_dicts, results):
    """ JJN: Calculates Choi et al's Medical Conceptual Relatedness Property, using includsion in a ICD9 system as the
    'concept type', that is, we examine whether embeddings that are related to an ICD9 medical system do indeed cluster
    """
    
    filename_to_embedding_matrix, idx_to_cui, cui_to_idx, cui_to_icd9_types = generate_overlapping_sets_cui(filenames_type, True, cui_to_icd9_dicts)
    
    # Obtain dictionary between cuis and the ICD9 disease systems they're part of/related to
    cui_to_systems = get_cui_to_systems(cui_to_icd9_types, icd9_systems)
    # And also build a similiar idx_to_system
    idx_to_systems = {}
    for idx in idx_to_cui:
        cui = idx_to_cui[idx]
        systems = cui_to_systems[cui]
        idx_to_systems[idx] = systems
    
    # Get list of ICD9 system names in this analysis
    icd9_systems_names = []
    for row in icd9_systems: icd9_systems_names.append(row[0])
    
    o = open(str(results_folder / 'debug_cui_to_systems' ), 'w')
    for cui in cui_to_systems.keys():
        o.write(cui + ': ')
        o.write(str(cui_to_systems[cui])+ '\n')
    
    filename_index = 0
    for filename, embedding_type, _ in filenames_type:
        # Matrix to convert cui to positions in the relevant filename
        embedding_matrix = filename_to_embedding_matrix[filename]
        Y = cdist(embedding_matrix, embedding_matrix, 'cosine')
        ranks = np.argsort(Y)
        systems_idx_dcg_err = {}
        systems_dcg = {}
        systems_err = {}
        systems_n = {} # Count of examples found for this system
        print 'done calcualting distance'
        
        system_index = 0
        for system in icd9_systems_names:
            systems_dcg[system] = []
            systems_err[system] = []
            systems_idx_dcg_err[system] = []
            systems_n[system] = 0
        for idx in idx_to_systems.keys():
            target = ranks[idx, 1:num_of_neighbor+1]
            for system in idx_to_systems[idx]:
                systems_n[system] += 1
                dcg = 0
                err = 0 
                for i in xrange(num_of_neighbor):
                    neighbor_idx = target[i]
                    neighbor_systems = idx_to_systems[neighbor_idx]
                    if system in neighbor_systems:
                        dcg += np.reciprocal(np.log2(i+2))
                        if err == 0:
                            err = 1/(1+i)
                systems_idx_dcg_err[system].append((idx, dcg, err))
                systems_dcg[system].append(dcg)
                systems_err[system].append(err)
        
        # Folowing to output the actual dcgs for further statistical analysis. Similiar could be used for stats on raw schools in other methods
        ## dcg_filename = 'choi_mcsp_dcg_' + filename + '.csv' 
        ##o2 = open(str(results_folder / dcg_filename), 'w')
        ## array_dcg = []
        
        system_index = 0
        for system in icd9_systems_names:
            results[system_index + 1][0] = re.sub(",", " ", system)
            results[system_index + 1][filename_index + 1] = '%2.5f +/-  %2.5f' %(np.mean(np.array(systems_dcg[system])), np.std(np.array(systems_dcg[system])))
            results[system_index + 1][-1] = str(systems_n[system]) # Number of examples used for this calculation. Will be re-written by each file but that's okay as always same
            ## array_dcg.append([str(x) for x in systems_dcg[system]])
            system_index += 1
        filename_index += 1
        
        ## write_results_to_file(array_dcg,o2)
        ##o2.close()  
        
    return results

def print_choi_mcsp(filenames, num_of_nn=40):
    """ JJN: Prints Choi et al's Medical Conceptual Relatedness Property by ICD9 system"""
    
    # Cui_to_icd9 mappings will be used
    cui_to_icd9 = get_icd9_cui_mappings_rangeok()
    # Create dictionaries linking drug cuis to the icd9 conditions they prevent or treat
    cui_icd9_tr, cui_icd9_pr = get_cui_may_treat_prevent_icd9(cui_to_icd9)
    # Store in a dict to pass
    cui_to_icd9_dicts = {}
    cui_to_icd9_dicts['cui_to_icd9'] = cui_to_icd9
    cui_to_icd9_dicts['cui_to_icd9_may_treat'] = cui_icd9_tr
    cui_to_icd9_dicts['cui_to_icd9_may_prevent'] = cui_icd9_pr
    cui_to_icd9_drug_or_diag = get_cui_to_icd9_drug_or_diag(cui_to_icd9, cui_icd9_tr, cui_icd9_pr)
    cui_to_icd9_dicts['cui_to_icd9_drug_or_diag'] = cui_to_icd9_drug_or_diag

    # Text file containing the system, start, end. Note that 'end' is an integer, so will end up to next integer
    icd9_systems_file = 'icd9_systems.txt'
    # Parse above file to get the system names, starts, ends
    icd9_systems = []
    with open(data_folder / icd9_systems_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            row_str = row.strip().split('|')
            if row_str[0] != 'all':
                icd9_systems.append([row_str[0], float(row_str[1]), float(row_str[2])])        
    
    # Make a numpy matrix to store results for easy printing
    cols = len(filenames) + 2 # System, Each Filename's MCSP, Examples per System
    rows = len(icd9_systems) + 1 # One for each system plus the header
    empty_results = [[0 for x in range(cols)] for y in range(rows)] 
    empty_results[0][0] = 'ICD9 Systems'
    empty_results[0][1:len(filenames) + 1] = list(map(lambda x: x[2], filenames))
    empty_results[0][-1] = 'Examples in System'
        
    print 'Choi Medical Conceptual Similarity Property by ICD9 system'
    results = get_choi_mcsp_by_systems(filenames, num_of_nn, icd9_systems, cui_to_icd9_dicts, empty_results)
    for line in results: print line
    
    choi_mrp_by_system = 'choi_mcsp_by_system.csv'
    o = open(str(results_folder / choi_mrp_by_system ), 'w')
    write_results_to_file(results,o)
    o.close()    
