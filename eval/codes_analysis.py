from __future__ import division
import argparse
import numpy as np
import scipy as sp
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from icd9 import ICD9
from pathlib import Path
import re

tree = ICD9('codes.json')
data_folder = Path("../data")
results_folder = Path("../results")

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


# JJN changed to support csv embedding files
def read_embedding_cui(filename):
    concept_to_cui, cui_to_concept = get_cui_concept_mappings() # comment out this after fix input
    cui_to_icd9, icd9_to_cui = get_icd9_cui_mappings()
    
    # Assume embedding files with .txt are deliminated with ' ' and ',' if .csv
    if filename.endswith('.txt'):
        delim = ' '
    elif filename.endswith('.csv'):
        delim = ','
    else:
        raise Exception('embedding file must be .txt or .csv depending on deliminator')
    
    with open(str(data_folder / filename), 'r') as infile:
        first_line = infile.readline().strip().split(delim)
        embedding_num = int(first_line[0])
        dimension = int(first_line[1])
        #embedding_num, dimension = map(int, infile.readline().strip().split(' '))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines() 
        idx_to_name = {}
        name_to_idx = {}
        embedding_type_to_indices = {}
        embedding_type_to_indices['IDX'] = []
        embedding_type_to_indices['O'] = []
        for idx in xrange(embedding_num-1):
            datum = data[idx+1].strip().split(delim)
            cui = datum[0]
            if cui[0] != 'C':
                if cui in concept_to_cui:
                    cui = concept_to_cui[cui].strip() # JJN Added to remove \n
            ##print len(datum[1:]) TODO REMOVE
            embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
            # potential bug here
            if cui in cui_to_icd9:
                idx_to_name[idx] = cui_to_icd9[cui]
                name_to_idx[cui_to_icd9[cui]] = idx
                embedding_type_to_indices['IDX'].append(idx)
            else:
                idx_to_name[idx] = cui
                name_to_idx[cui] = idx
                embedding_type_to_indices['O'].append(idx)
        return embedding_matrix, embedding_type_to_indices, name_to_idx, idx_to_name
    

def read_embedding_codes(filename):
    with open(str(data_folder / filename), 'r') as infile:
        embedding_num, dimension = map(int, infile.readline().strip().split(' '))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines()
        embedding_type_to_indices = {}
        name_to_idx = {}
        idx_to_name = {}
        for idx in xrange(embedding_num-1):
            datum = data[idx+1].strip().split(' ')
            embedding_name = datum[0]
            embedding_type, embedding_value = embedding_name.split('_')
            name_to_idx[embedding_value] = idx
            idx_to_name[idx] = embedding_value 
            if embedding_type not in embedding_type_to_indices:
                embedding_type_to_indices[embedding_type] = []
            embedding_type_to_indices[embedding_type].append(idx)
            embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
        return embedding_matrix, embedding_type_to_indices, name_to_idx, idx_to_name

def generate_overlapping_sets_icd9(filenames_type):
    embedding_idx_icd9 = {} # a dictionary of (embedding_matrix, idx_to_icd9, icd9_to_idx)
    overlapping_icd9s = set([])
    start = 1

    for filename, embedding_type, _ in filenames_type:
        #print filename
        #print embedding_type
        if embedding_type == 'codes':
            embedding_matrix, embedding_type_to_indices, icd9_to_idx, idx_to_icd9 = read_embedding_codes(filename)
            embedding_idx_icd9[filename] = (embedding_matrix, idx_to_icd9, icd9_to_idx)
            if start == 1:
                start = 0
                overlapping_icd9s.update(set(icd9_to_idx.keys()))
            else:
                overlapping_icd9s.intersection_update(set(icd9_to_idx.keys()))
                #print len(overlapping_icd9s)
        elif embedding_type == 'cui':
            embedding_matrix, embedding_type_to_indices, icd9_to_idx, idx_to_icd9 = read_embedding_cui(filename)
            embedding_idx_icd9[filename] = (embedding_matrix, idx_to_icd9, icd9_to_idx)
            if start == 1:
                start = 0
                overlapping_icd9s.update(set(icd9_to_idx.keys()))
            else:
                overlapping_icd9s.intersection_update(set(icd9_to_idx.keys()))
                #print len(overlapping_icd9s)
    overlapping_icd9s = list(overlapping_icd9s)

    idx_of_overlapping_icd9s = {}
    for filename, embedding_type, _ in filenames_type:
        idx_of_overlapping_icd9s[filename] = []

    idx_to_icd9 = {}
    icd9_to_idx = {}
    for idx, icd9 in enumerate(overlapping_icd9s):
        idx_to_icd9[idx] = icd9
        icd9_to_idx[icd9] = idx
        for filename, embedding_type, _ in filenames_type:
            idx_of_overlapping_icd9s[filename].append(embedding_idx_icd9[filename][2][icd9])

    filename_to_embedding_matrix = {}
    for filename, embedding_type, _ in filenames_type:
        idx_of_overlapping_icd9s[filename] = np.array(idx_of_overlapping_icd9s[filename]).astype(int) #JJN added due to change in numpy
        filename_to_embedding_matrix[filename] = embedding_idx_icd9[filename][0][idx_of_overlapping_icd9s[filename]]
    return filename_to_embedding_matrix, idx_to_icd9, icd9_to_idx

def generate_overlapping_sets_cui(filenames_type):    
    filenames = [filename for filename, _, _2 in filenames_type]
    
    embedding_idx_cui = {} # a dictionary of (embedding_matrix, idx_to_cui, cui_to_idx)
    overlapping_cuis = set([])

    if len(filenames) == 1:
        embedding_matrix, idx_to_cui, cui_to_idx = read_embedding_matrix_cui(filenames[0])
        filename_to_embedding_matrix = {}
        filename_to_embedding_matrix[filenames[0]] = embedding_matrix
        return filename_to_embedding_matrix, idx_to_cui, cui_to_idx

    for fileidx, filename in enumerate(filenames):
        embedding_matrix, idx_to_cui, cui_to_idx = read_embedding_matrix_cui(filename)
        embedding_idx_cui[filename] = (embedding_matrix, idx_to_cui, cui_to_idx)
        if fileidx == 0:
            overlapping_cuis.update(set(cui_to_idx.keys()))
        else:
            overlapping_cuis.intersection_update(set(cui_to_idx.keys()))
    overlapping_cuis = list(overlapping_cuis)
    
    idx_of_overlapping_cuis = {}
    for filename in filenames:
        idx_of_overlapping_cuis[filename] = []

    idx_to_cui = {}
    cui_to_idx = {}
    for idx, cui in enumerate(overlapping_cuis):
        idx_to_cui[idx] = cui
        cui_to_idx[cui] = idx
        for filename in filenames:
            idx_of_overlapping_cuis[filename].append(embedding_idx_cui[filename][2][cui])
    filename_to_embedding_matrix = {}
    for filename in filenames:
        idx_of_overlapping_cuis[filename] = np.array(idx_of_overlapping_cuis[filename]).astype(int) #JJN added due to change in numpy
        filename_to_embedding_matrix[filename] = embedding_idx_cui[filename][0][idx_of_overlapping_cuis[filename]]
    return filename_to_embedding_matrix, idx_to_cui, cui_to_idx

# JJN changed to support csv embedding files
def read_embedding_matrix_cui(filename):
    concept_to_cui, cui_to_concept = get_cui_concept_mappings() # comment out this after fix input
    
    # Assume embedding files with .txt are deliminated with ' ' and ',' if .csv
    if filename.endswith('.txt'):
        delim = ' '
    elif filename.endswith('.csv'):
        delim = ','
    else:
        raise Exception('embedding file must be .txt or .csv depending on deliminator')
    

    with open(str(data_folder / filename), 'r') as infile:            
        first_line = infile.readline().strip().split(delim)
        embedding_num = int(first_line[0])
        dimension = int(first_line[1])
        #embedding_num, dimension = map(int, infile.readline().strip().split(delim))
        # -1 for remove </s>
        embedding_matrix = np.zeros((embedding_num-1, dimension))
        data = infile.readlines() 
        idx_to_cui = {}
        cui_to_idx = {}
        for idx in xrange(embedding_num-1):
            datum = data[idx+1].strip().split(delim)
            cui = datum[0]
            if cui[0] != 'C':
                if cui in concept_to_cui:
                    cui = concept_to_cui[cui].strip() # JJN added to remove \n
            embedding_matrix[idx,:] = np.array(map(float, datum[1:]))
            idx_to_cui[idx] = cui
            cui_to_idx[cui] = idx
        return embedding_matrix, idx_to_cui, cui_to_idx


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


# type == f: fine grain, c: coarse grain
def get_css_analysis(filenames_type, num_of_neighbor, type='f'):
    filename_to_embedding_matrix, idx_to_icd9, icd9_to_idx = generate_overlapping_sets_icd9(filenames_type)
    print len(icd9_to_idx.keys())
    if type == 'c':
        icd9_pairs = get_coarse_icd9_pairs(set(icd9_to_idx.keys()))
    else:
        icd9_pairs = get_icd9_pairs(set(icd9_to_idx.keys()))

    print len(icd9_pairs)
    icd9_to_check = set(icd9_pairs.keys())
    icd9_to_check.intersection_update(set(icd9_to_idx.keys()))

    print len(icd9_to_check)

    icd9_to_description = get_icd9_to_description()
    for icd9 in icd9_to_idx.keys():
        if icd9 not in icd9_to_description:
            if tree.find(icd9):
                icd9_to_description[icd9] = tree.find(icd9).description.encode('utf-8')
            else:
                icd9_to_description[icd9] = ''


    filename_all = []
    value_all = []
    for filename, embedding_type, _ in filenames_type:
        #print filename
        icd9_embeddings = filename_to_embedding_matrix[filename]
        Y = cdist(icd9_embeddings, icd9_embeddings, 'cosine')
        ranks = np.argsort(Y)
    
        cumulative_ndcgs = []

        for icd9 in icd9_to_check:
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
    return filename_all, value_all


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

def get_cui_may_treat_prevent_icd9(cui_to_icd9):
    """JJN: Builds two dictionaries. Keys are the cuis repersenting drugs. 
    Values are the ICD9 codes of the conditions they may treat, or may prevent 
    (for the two dictionaries respectively)
    """
    cui_to_icd9_may_treat = {}
    cui_to_icd9_may_prevent = {}
    
    # Find diagonses that may be treated by drugs
    with open(str(data_folder / 'may_treat_cui.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            drug, diseases = row.strip().split(':')
            diseases_cui = diseases.split(',')[:-1]
            diseases_icd9 = []
            for disease_cui in diseases_cui:
                if disease_cui in cui_to_icd9.keys():
                    diseases_icd9.append(cui_to_icd9[disease_cui])
                else:
                    print "Warning, this cui not found in cui_to_icd9: " + disease_cui
            cui_to_icd9_may_treat[drug] = diseases_icd9
    
    # Find diagnoses that may be prevent by drugs
    with open(str(data_folder / 'may_prevent_cui.txt'), 'r') as infile:
        data = infile.readlines()
        for row in data:
            diseases_cui = diseases.split(',')[:-1]
            drug, diseases = row.strip().split(':')
            diseases_icd9 = []
            for disease_cui in diseases_cui:
                if disease_cui in cui_to_icd9.keys():
                    diseases_icd9.append(cui_to_icd9[disease_cui])
                else:
                    print "Warning, this cui not found in cui_to_icd9: " + disease_cui
            cui_to_icd9_may_prevent[drug] = diseases_icd9

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

    filename_all = []
    value_all = []
    
    ## TODO REMOVE THIS JUST FOR DEBUGGING
    o = open(str(results_folder / 'overlappingcuis.txt' ), 'w')
    for cui in cuis: o.write(str(cui)+'\n')
    o.close()
    ## -----------------------------------
    n_drugs_in_system = 0 # Number of comparisons with a drug found in this system
    n_diags_in_system = 0 # Number of comparisons with a diag found in this system
    missing_cuis = []
    
    for filename, embedding_type, _ in filenames_type:
        #Contains the unmsrs scores from comparisons that have to do with this system
        unmsrs_scores = []
        #Contains the vector cosine similarity between the embedding vectors
        veccos_scores = []
        #Number of relevant judgements
        comparisons_in_cuis = 0

        for row_str in umnsrs_rows:
            row = row_str.strip().split(',')
            umnsrs_rating = row[0]
            cui_1 = row[4]
            cui_2 = row[5]
            
            if (cui_1 in cuis) and (cui_2 in cuis):
                comparisons_in_cuis += 1
                in_system_diag_1, in_system_drug_1 = cui_in_system(cui_1, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
                in_system_diag_2, in_system_drug_2 = cui_in_system(cui_2, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
                if in_system_diag_1 or in_system_drug_1 or in_system_diag_2 or in_system_drug_2:
                    cos_sim = 0 # Calculate Cosine similiarity between vectors
                    
                    
                    
                    
                    
                    veccos_scores.append(cos_sim)
                    unmsrs_scores.append(umnsrs_rating)
                    
                    # Keep track of hits from drugs or diag
                    if in_system_diag_1 or in_system_diag_2:
                        n_diags_in_system += 1
                    if in_system_drug_1 or in_system_drug_2:
                        n_drugs_in_system += 1
                    
            elif cui_1 in cuis:
                missing_cuis.append(cui_2)
            else:
                missing_cuis.append(cui_1)
        #print(set(missing_cuis))
        
            
        rho, pval = spearmanr(unmsrs_scores,veccos_scores)
        filename_all.append((filename))
        value_all.append(rho)        
            
    return filename_all, value_all, comparisons_in_cuis, n_diags_in_system, n_drugs_in_system
    


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
    
    print 'Yu Spearman Correlation with UMNSRS ratings by ICD9 system'
    for system in icd9_systems:
        system_name = system[0]
        start = float(system[1])
        end = float(system[2])
        
        filename_to_print, ndcgs_to_print, comparisons_in_cuis, n_diags_in_system, n_drugs_in_system = get_yu_umnsrs_cor_by_system(filenames, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9)
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
        print "Number of comparisons with both cuis present: " + str(comparisons_in_cuis)
        print "Number of lines with a drug found in this system: " + str(n_drugs_in_system)
        print "Number of lines with a diagnoses found in this system: " + str(n_diags_in_system)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filenames", default='orig_files_all.txt')
    args = parser.parse_args()
    filenames_file = args.filenames

    filenames = []
    with open(filenames_file, 'r') as infile:
        data = infile.readlines()
        for row in data:
            filenames.append(row.strip().split(','))
    
    
    
    ##filename_to_print, finegrain_ndcgs = get_css_analysis(filenames, num_of_nn, 'f')
    ##filename_to_print, coarsegrain_ndcgs = get_css_analysis(filenames, num_of_nn, 'c')
    ##for name, finegrain, coarsegrain in zip(filename_to_print, finegrain_ndcgs, coarsegrain_ndcgs):
    ##    print '%s & %.2f & %.2f \\\\' %(name.split('/')[-1], finegrain*100, coarsegrain*100)
    
    #print_choi_mrp(filenames)
    
    print_yu_umnsrs_cor(filenames)
    
    
    ##filename_to_print, psych_ndcgs = get_choi_mrp_by_system(filenames, num_of_nn, 'c')
    ##filename_to_print, coarsegrain_ndcgs = get_css_analysis(filenames, num_of_nn, 'c')
    ##for name, psychs, coarsegrain in zip(filename_to_print, psych_ndcgs, coarsegrain_ndcgs):
    ##    print '%s & %.2f & %.2f \\\\' %(name.split('/')[-1], psychs*100, coarsegrain*100)
