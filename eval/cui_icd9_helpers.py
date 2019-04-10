# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:11:42 2019

@author: jjnun
"""
from __future__ import division
##import argparse
##import numpy as np
##import scipy as sp
##from scipy.spatial.distance import cdist, cosine
##from scipy.stats import spearmanr
from icd9 import ICD9
from pathlib import Path
##import re

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")

def get_icd9_to_description():
    icd9_to_description = {}
    with open(str(data_folder / 'CMS32_DESC_LONG_DX.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            icd9 = row.strip()[:6].strip()
            if len(icd9) > 3:
                icd9 = icd9[:3] + '.' + icd9[3:]
            description = row.strip()[6:].strip()
            icd9_to_description[icd9] = description
    return icd9_to_description

def get_icd9_cui_mappings_rangeok():
    """ Modified version of original. Some cui correspond to a range of icd9 codes
    eg C0041296 is Tuberculosis, which is ICD9 range 010-018.99. This original method
    just ignored these codes. For our purpose of figuring out if a cui is in a broad system
    cateogry, we can simply choose a value within this range (now, midpoint). This is useful for cui --> icd9
    but not the reverse, so we'll now only export the useful direction. 
    """
    cui_to_icd9 = {}
    with open(str(data_folder / 'cui_icd9.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            ele = row.strip().split('|')
            if ele[11] == 'ICD9CM':
                cui = ele[0]
                icd9 = ele[10]
                if cui not in cui_to_icd9 and icd9 != '' and not any(c.isalpha() for c in icd9): #Ignore icd9's with alphabet in them not part of our study
                    # Check if this cui is a range of ICD9 
                    if '-' in icd9:
                        icd9range = icd9.split('-')
                        beg = float(icd9range[0])
                        end = float(icd9range[1])
                        cui_to_icd9[cui] = '%.2f' % ((beg+end)/2)
                        # else: # A select few ICD9 start with a letter, if so just choose first part of range
                        #    icd9range = icd9.split('-')
                        #    cui_to_icd9[cui] = icd9range[0]
                    else:
                        cui_to_icd9[cui] = icd9
    return cui_to_icd9

def get_cui_may_treat_prevent_icd9(cui_to_icd9):
    """JJN: Builds two dictionaries. Keys are the cuis repersenting drugs. 
    Values are arrays containing dicts which have the cui and ICD9 code of a disease treated
    or prevented (for the two dictionaries respectively. Takes in a dict of cui_to_icd9 repersenting
    the list of cuis of diagnoses that have known icd9 relations). 
    """
    cui_to_icd9_may_treat = {}
    cui_to_icd9_may_prevent = {}
    
    # Find diagonses that may be treated by drugs
    with open(str(data_folder / 'may_treat_cui.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            drug, diseases = row.strip().split(':')
            diseases_cui = diseases.split(',')[:-1]
            treated_diseases = []
            for disease_cui in diseases_cui:
                if disease_cui in cui_to_icd9.keys():
                    treated_disease = {}
                    treated_disease['cui'] = disease_cui
                    treated_disease['icd9'] = cui_to_icd9[disease_cui]
                    treated_diseases.append(treated_disease)
                else:
                    print "Warning, this treated cui not found in cui_to_icd9: " + disease_cui
            if len(treated_diseases) != 0: # Only add entries that treat an ICD9 code
                cui_to_icd9_may_treat[drug] = treated_diseases
            
    
    # Find diagnoses that may be prevent by drugs
    with open(str(data_folder / 'may_prevent_cui.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            drug, diseases = row.strip().split(':')
            diseases_cui = diseases.split(',')[:-1]
            prevented_diseases = []
            ##if drug == 'C0066685':
                ##print 'inisde get_cui_may_treat_prevent'
                ##print 'here is the drug: ' + drug
                ##print 'here is the diseases: ' + diseases
                ##print 'here is the diseases split: ' + str(diseases_cui)
            for disease_cui in diseases_cui:
                if disease_cui in cui_to_icd9.keys():
                    prevented_disease = {}
                    prevented_disease['cui'] = disease_cui
                    prevented_disease['icd9'] = cui_to_icd9[disease_cui]
                    prevented_diseases.append(prevented_disease)
                    ##if drug == 'C0066685':
                    ##    print 'here is the prevented_disease: ' + str(prevented_disease)
                else:
                    print "Warning, this prevented cui not found in cui_to_icd9: " + disease_cui
            if len(prevented_diseases) != 0:
                cui_to_icd9_may_prevent[drug] = prevented_diseases

    return cui_to_icd9_may_treat, cui_to_icd9_may_prevent

def cui_in_system(cui, start, end, cui_icd9_treat, cui_icd9_prevent, cui_to_icd9):
    """JJN: finds whether a given CUI is in an ICD9 system. First, assumes its a diagnossis, looking up if CUI is a ICD9 diagnoses, if so returns whether in system. 
    If not, assumes it is a drug, looks up if CUI is found in a list of may-prevent or may-treat relations.
    If found, looks into whether these may-treat or may-prevent lists include a diagnosis with an ICD9 code within given range"""
    
    in_system_diag = False
    in_system_drug = False
    
    # Try to see if cui is a diagnosis  
    if cui in cui_to_icd9.keys():    
        icd9 = cui_to_icd9[cui]
        in_system_diag = start <= float(icd9) < end + 1
        
    # If not, determines if it may treat or prevent a disease in ICD9 system
    if cui in cui_icd9_treat.keys():
        for icd9 in cui_icd9_treat[cui]:
            if start <= float(icd9) < end + 1: in_system_drug = True
    
    if cui in cui_icd9_prevent.keys():
        for icd9 in cui_icd9_prevent[cui]:
            if start <= float(icd9) < end + 1: in_system_drug = True
    
    return in_system_diag, in_system_drug

def get_cui_to_icd9_drug_or_diag(cui_to_icd9, cui_icd9_tr, cui_icd9_pr):
    """ JJN: Makes a dictionary with keys all possible cuis, and the values
    a dictionary for each cui, with the icd9_type (diag if ICD9 diagnosis, drug
    if it treats or prevents an ICD9 diagnosis), and an array of the ICD9(s)
    and an array of the treated/prevent cuis (empty if it's a diag)
    """
    cui_to_icd9_drug_or_diag = {}
    
    # Unpack cui to icd9 relations from dict
    ##cui_to_icd9 = cui_to_icd9_dicts['cui_to_icd9']
    ##cui_to_icd9_may_treat = cui_to_icd9_dicts['cui_to_icd9_may_treat']
    ##cui_to_icd9_may_prevent = cui_to_icd9_dicts['cui_to_icd9_may_prevent']
    
    ##print 'Len of cui icd9 may tr: ' + str(len(cui_to_icd9_may_treat))
    ##print 'Len of cui icd9 may pr: ' + str(len(cui_to_icd9_may_prevent)) 
    #if cui == 'C0066685':
    #    print 'Heres the treat dict at the start'
    #    print cui_to_icd9_may_treat[cui]
    #    print 'Heres the prevent dict at the start'
    #    print cui_to_icd9_may_prevent[cui]
    #    print 'end of start'

    # First, make an entry for every cui that's a diagnosis (in cui_to_icd9)
    for diag_cui in cui_to_icd9.keys():
        cui_dict = {}
        cui_dict['icd9_type'] = 'diag'
        cui_dict['icd9s'] = [cui_to_icd9[diag_cui]]
        cui_dict['cuis'] = []
        cui_to_icd9_drug_or_diag[diag_cui] = cui_dict
    
    tr_cuis = cui_icd9_tr.keys() # 
    pr_cuis = cui_icd9_pr.keys()
    all_cuis = list(set(tr_cuis + pr_cuis))
    
    for drug_cui in all_cuis:
        if drug_cui == 'C0066685XXX':
            print 'Heres the treat dict at the start'
            print cui_icd9_tr[drug_cui]
            print 'Heres the prevent dict at the start'
            print cui_icd9_pr[drug_cui]
            print 'end of start'
        
        if drug_cui in tr_cuis and drug_cui in pr_cuis:
            tr_or_pr_diseases = combine_dicts(cui_icd9_tr[drug_cui], cui_icd9_pr[drug_cui])
            icd9s = []
            cuis = []
            for tr_or_pr_disease in tr_or_pr_diseases:
                cuis.append(tr_or_pr_disease['cui'])
                icd9s.append(tr_or_pr_disease['icd9'])
            cui_dict = {}
            cui_dict['icd9_type'] = 'drug'
            cui_dict['icd9s'] = icd9s
            cui_dict['cuis'] = cuis
            cui_to_icd9_drug_or_diag[drug_cui] = cui_dict
            
            if drug_cui == 'C0066685XXX':
                print 'Heres the icd9s after combination: '
                print icd9s
                print 'Heres the cuis after combination: '
                print cuis
                print 'Heres the whole dict: '
                print cui_dict
                print "end------"
            
        elif drug_cui in tr_cuis:
            icd9s = []
            cuis = []
            for tr_disease in cui_icd9_tr[drug_cui]:
                cuis.append(tr_disease['cui'])
                icd9s.append(tr_disease['icd9'])            
            cui_dict = {}
            cui_dict['icd9_type'] = 'drug'
            cui_dict['icd9s'] = icd9s
            cui_dict['cuis'] = cuis
            cui_to_icd9_drug_or_diag[drug_cui] = cui_dict
        elif drug_cui in pr_cuis:
            cui_dict = {}
            
            icd9s = []
            cuis = []
            for pr_disease in cui_icd9_pr[drug_cui]:
                cuis.append(pr_disease['cui'])
                icd9s.append(pr_disease['icd9'])
            cui_dict['icd9_type'] = 'drug'
            cui_dict['icd9s'] = icd9s
            cui_dict['cuis'] = cuis
            cui_to_icd9_drug_or_diag[drug_cui] = cui_dict
        
    return cui_to_icd9_drug_or_diag

# =============================================================================
#         
#     
#     #Merge the lists of which cuis treat or prevent which icd9s
#     cui_to_icd9_may_treat_or_prevent = cui_to_icd9_may_treat
#     for a_cui in cui_to_icd9_may_prevent.keys():
#         if a_cui in cui_to_icd9_may_treat.keys():
#             if a_cui == 'C0066685XXX':
#                 print 'We found this cui in both: - ' + a_cui
#                 print 'Len of cui icd9 may tr: ' + str(len(cui_to_icd9_may_treat))
#                 print 'Len of cui icd9 may pr: ' + str(len(cui_to_icd9_may_prevent))
#                 print cui_to_icd9_may_treat[a_cui]
#                 print cui_to_icd9_may_prevent[a_cui]
#                 print 'end of this stuff'
#                 combined_tr_pr_diseases = combine_dicts(cui_to_icd9_may_treat[a_cui], cui_to_icd9_may_prevent[a_cui])
#             combined_tr_pr_diseases = []
#             cui_to_icd9_may_treat_or_prevent[a_cui] = combined_tr_pr_diseases
#         else:
#             cui_to_icd9_may_treat_or_prevent[a_cui] = cui_to_icd9_may_prevent[a_cui]
#         
#     elif cui in cui_to_icd9_may_treat_or_prevent:
#         # This step done is unison to ensure position of cuis and icd9s are same
#         # so they repersent the same disease being treated/prevented
#         diseases_tr_or_pr = cui_to_icd9_may_treat_or_prevent[cui]
#         ##print 'Length of diseases_tr_or_pr: ' + str(len(diseases_tr_or_pr))
#         icd9s = []
#         cuis = []
#         for disease_tr_or_pr in diseases_tr_or_pr:
#             cuis.append(disease_tr_or_pr['cui'])
#             icd9s.append(disease_tr_or_pr['icd9'])
#         assert len(icd9s) != 0, 'No icd9s found for: ' + cui
#         return 'drug', icd9s, cuis  # Keep only uniques
#     else: 
#         return 'none', [], []
# =============================================================================

def combine_dicts(list_dicts_1, list_dicts_2):
    """ JJN: Quick helper function that takes in two lists of dicts and returns
    a list of dicts. If a dict is only in one list, is brought to new list. If 
    a dict is in both, for each key, values are a unique set of values from 
    both dicts. In this implementation, the two lists are dicts of a disease that is
    treated or prevented by a given drug, and so I will be explcit for keys. Could
    later be changed to be general if needed. 
    """
    combined_dicts = list_dicts_2
    
    for dict_1 in list_dicts_1:
        cui_1 = dict_1['cui']
        icd9_1 = dict_1['icd9']
        
        in_other_dict = False
        # Check if in other list of dict
        for dict_2 in list_dicts_2:
            cui_2 = dict_2['cui']
            icd9_2 = dict_2['icd9']
            
            if cui_1 == cui_2:
                assert icd9_1 == icd9_2, 'A disease dict should have icd9 and cui the same'
                in_other_dict = True
        
        if not in_other_dict:
            combined_dicts.append(dict_1)
    
    # Check
    if len(combined_dicts) > 99:
        print 'Check to see this was combined correctly: '
        print '---'
        print list_dicts_1
        print '---'
        print list_dicts_2
        print '---'
        print combined_dicts
        print '---'
    
    #--
    return combined_dicts    

def get_icd9_pairs(icd9_set):
    icd9_pairs = {}
    with open(str(data_folder / 'icd9_grp_file.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            codes, name = row.strip().split('#')
            name = name.strip()
            codes = codes.strip().split(' ')
            new_codes = set([])
            for code in codes:
                if code in icd9_set:
                    new_codes.add(code)
                elif len(code) > 5 and code[:5] in icd9_set:
                    new_codes.add(code[:5])
                elif len(code) > 4 and code[:3] in icd9_set:
                    new_codes.add(code[:3])
            codes = list(new_codes)

            if len(codes) > 1:
                for idx, code in enumerate(codes):
                    if code not in icd9_pairs:
                        icd9_pairs[code] = set([])
                    icd9_pairs[code].update(set(codes[:idx]))
                    icd9_pairs[code].update(set(codes[idx+1:]))
    return icd9_pairs


def get_coarse_icd9_pairs(icd9_set): 
    icd9_pairs = {}
    ccs_to_icd9 = {}
    with open(str(data_folder / 'ccs_coarsest.txt'), 'r') as infile:
        data = infile.readlines()
        currect_ccs = ''
        for row in data:
            if row[:10].strip() != '':
                current_ccs = row[:10].strip()
                ccs_to_icd9[current_ccs] = set([])
            elif row.strip() != '':
                ccs_to_icd9[current_ccs].update(set(row.strip().split(' ')))

    ccs_coarse = {}
    for ccs in ccs_to_icd9.keys():
        ccs_eles = ccs.split('.')
        if len(ccs_eles) >= 2:
            code = ccs_eles[0] + '.' + ccs_eles[1]
            if code not in ccs_coarse:
                ccs_coarse[code] = set([])
            ccs_coarse[code].update(ccs_to_icd9[ccs])
    
    for ccs in ccs_coarse.keys():
        new_codes = set([])
        for code in ccs_coarse[ccs]:
            if len(code) > 3:
                new_code = code[:3] + '.' + code[3:]
            code = new_code
            if code in icd9_set:
                new_codes.add(code)
            elif len(code) > 5 and code[:5] in icd9_set:
                new_codes.add(code[:5])
            elif len(code) > 4 and code[:3] in icd9_set:
                new_codes.add(code[:3])
        codes = list(new_codes)
        if len(codes) > 1:
            for idx, code in enumerate(codes):
                if code not in icd9_pairs:
                    icd9_pairs[code] = set([])
                icd9_pairs[code].update(set(codes[:idx]))
                icd9_pairs[code].update(set(codes[idx+1:]))
    return icd9_pairs

def get_cui_concept_mappings():
    concept_to_cui_hdr = str(data_folder / '2b_concept_ID_to_CUI.txt')
    concept_to_cui = {}
    cui_to_concept = {}
    with open(concept_to_cui_hdr, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            concept = line.split('\t')[0]
            cui = line.split('\t')[1].split('\r')[0]
            concept_to_cui[concept] = cui 
            cui_to_concept[cui] = concept
    return concept_to_cui, cui_to_concept


def get_icd9_cui_mappings():
    cui_to_icd9 = {}
    icd9_to_cui = {}
    with open(str(data_folder / 'cui_icd9.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            ele = row.strip().split('|')
            if ele[11] == 'ICD9CM':
                cui = ele[0]
                icd9 = ele[10]
                if cui not in cui_to_icd9 and icd9 != '' and '-' not in icd9:
                    cui_to_icd9[cui] = icd9
                    icd9_to_cui[icd9] = cui
    return cui_to_icd9, icd9_to_cui

def get_cui_to_systems(cui_to_icd9_types, icd9_systems):
    """
    Takes in a dictionary of cuis containing their ICD9 relationship (diag is they're
    directly a ICD9 diagnosis, drug is that they may prevent or may treat an ICD9 diagnosis)
    and outputs a new cui dictionary containing the systems these ICD9 relationships are within
    """
    cui_to_systems = {}

    for cui in cui_to_icd9_types.keys():
        ## print "This should be cui: " + cui
        systems = []
        cui_dict = cui_to_icd9_types[cui] 
        ## print cui_dict
        icd9s = cui_dict['icd9s']
        ## print 'Here are some icd9s: ' + str(icd9s)
        ##if len(icd9s) == 0: print '0 len icd9s for: ' + cui
        for icd9 in icd9s:
            systems_len = len(systems)
            for system in icd9_systems:
                system_name  = system[0]
                system_start = system[1]
                system_end   = system[2]
                if system_start <= float(icd9) < system_end + 1:
                    systems.append(system_name)
            if systems_len == len(systems): # Catch if no systems added
                print 'icd9 not found in a system: ' + icd9
        cui_to_systems[cui] = list(set(systems)) # Keep only unique
    
    return cui_to_systems
        
