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
from cui_icd9_helpers import get_coarse_icd9_pairs, get_icd9_pairs, get_icd9_to_description, get_cui_to_systems
from cui_icd9_helpers import cui_in_system, get_icd9_cui_mappings_rangeok, get_cui_may_treat_prevent_icd9

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")


def write_results_to_file(results,f):
    for row in results:
        f.write(','.join(row) + '\n')


def get_choi_mcsp_by_systems(filenames_type, num_of_neighbor, icd9_systems, cui_to_icd9_dicts, results):
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
    
    
    ##print "Here are some keys" + str(cui_to_systems.keys()[0:5])
    ##print "Here are some values" + str(cui_to_systems.values()[0:5])
    ##print "Here's one example of a cui to system: " + str(cui_to_systems['C0206180'])
    
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
        systems_examples = {} # Count of examples found for this system
        print 'done calcualting distance'
        
        system_index = 0
        for system in icd9_systems_names:
            systems_dcg[system] = []
            systems_err[system] = []
            systems_idx_dcg_err[system] = []
            systems_examples[system] = 0
        for idx in idx_to_systems.keys():
            target = ranks[idx, 1:num_of_neighbor+1]
            for system in idx_to_systems[idx]:
                systems_examples[system] += 1
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
        
        ##Temporary
        system_index = 0
        for system in icd9_systems_names:
            results[system_index + 1][0] = re.sub(",", " ", system)
            results[system_index + 1][filename_index + 1] = '%2.5f +/-  %2.5f' %(np.mean(np.array(systems_dcg[system])), np.std(np.array(systems_dcg[system])))
            results[system_index + 1][-1] = str(systems_examples[system]) # Number of examples used for this calculation. Will be re-written by each file but that's okay as always same
            ##print '%50s (DCG) %2.5f %2.5f' %(system, np.mean(np.array(systems_dcg[system])), np.std(np.array(systems_dcg[system])))
            ##print '%50s (ERR) %2.5f %2.5f' %(system, np.mean(np.array(systems_err[system])), np.std(np.array(systems_err[system])))
            ##f.write('%50s (DCG) %2.5f %2.5f\n' %(type, np.mean(np.array(type_dcg[type])), np.std(np.array(type_dcg[type]))))
            ##f.write('%50s (ERR) %2.5f %2.5f\n' %(type, np.mean(np.array(type_err[type])), np.std(np.array(type_err[type]))))
            system_index += 1
        filename_index += 1
        
    return results



# JJN: Prints the Medical Relatedness Property by ICD9 system
def print_choi_mcsp(filenames, num_of_nn=40):
    # Cui_to_icd9 mappings will be used
    cui_to_icd9 = get_icd9_cui_mappings_rangeok()
    # Create dictionaries linking drug cuis to the icd9 conditions they prevent or treat
    cui_icd9_tr, cui_icd9_pr = get_cui_may_treat_prevent_icd9(cui_to_icd9)
    # Store in a dict to pass
    cui_to_icd9_dicts = {}
    cui_to_icd9_dicts['cui_to_icd9'] = cui_to_icd9
    cui_to_icd9_dicts['cui_to_icd9_may_treat'] = cui_icd9_tr
    cui_to_icd9_dicts['cui_to_icd9_may_prevent'] = cui_icd9_pr
    
    
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
    empty_results = [[0 for x in range(cols)] for y in range(rows)] 
    empty_results[0][0] = 'ICD9 Systems'
    empty_results[0][1:len(filenames) + 1] = list(map(lambda x: x[2], filenames))
    empty_results[0][-1] = 'Examples in System'
    ## print empty_results
    # csv file to write results to
    
    ##o.write('ICD9 System,')
    # Write headers from 3rd entry in orig_files_all.txt
    ##o.write(",".join(list(map(lambda x: x[2], filenames))))
    # Write Examples per System
    ##o.write(", Examples in System")
    # Write DCG, Accuracy
    ##o.write(", DCG")
    ##o.write(", Accuracy")
    
    
    
    print 'Choi Medical Conceptual Similarity Property by ICD9 system'
    results = get_choi_mcsp_by_systems(filenames, num_of_nn, icd9_systems, cui_to_icd9_dicts, empty_results)
    for line in results: print line
    
    choi_mrp_by_system = 'choi_mcsp_by_system.csv'
    o = open(str(results_folder / choi_mrp_by_system ), 'w')
    write_results_to_file(results,o)
    o.close()    
    
    
    #for system in icd9_systems:
    #    system_name = system[0]
    #    start = float(system[1])
    #    end = float(system[2])
        
        #filename_to_print, ndcgs_to_print, comparisons_in_cuis = get_choi_mcsp_by_system(filenames, num_of_nn, start, end, cui_to_icd9_dicts)
        # Write ncdgs to file
#        ndcgs_rounded = [round(x*100,2) for x in ndcgs_to_print]
#        ncdgs_str = ','.join(map(str, ndcgs_rounded))
#        o.write('\n' + re.sub(",", " ", system_name) + ',') # Replace commas with space to use as csv
#        o.write(ncdgs_str)
#        o.write(", " + str(comparisons_in_cuis))
#        # Print ncdfs 
#        print '\n' + system_name
#        for file_name, ndcg in zip(filename_to_print, ndcgs_to_print):
#            print '%s & %.2f \\\\' %(file_name.split('/')[-1], ndcg*100)
#        print "Number of examples: " + str(comparisons_in_cuis)