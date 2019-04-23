# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:16:47 2019

@author: jjnun
"""

def get_new_sysvec_by_system(filenames_type, icd9_systems, cui_to_icd9_dicts, results):
    """ For each medical system, calculates a repersentative vector based on the
    mean of normalized vecotrs repersenting medical conditions in the system. 
    Then caluclates the percentage of drugs known to treat or prevent a disease in
    this system who are closest to the correct system's vector. For drugs that 
    treat/prevent conditions in n>1 systems, scores as a correct prediction if
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
             
                 cos_sims = np.zeros(n_of_systems)
             
                # Generate list of cos similarities with the system vectors
                 for i in range(n_of_systems):
                    system = systems_sysvec.keys()[i]
                    system_vec = systems_sysvec[system]
                    cos_sim = cosine_similarity([cui_vec], [system_vec])[0,0]
                    cos_sims[i] = cos_sim
             
                 n = len(true_systems) # Number of systems this cui treats or prevents or 1 if diagnosis
                 pred_systems = [icd9_systems_names[i] for i in np.argsort(cos_sims)[n:]]
                
                 for system in true_systems:
                    systems_n[system] += 1
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

#-----------------------------------------------------------------------------

def get_yu_umnsrs_cor_by_system(filenames_type, start, end, cui_icd9_tr, cui_icd9_pr, cui_to_icd9):
    """JJN: Calculates the Spearman Correlation Coefficient between UMNSRS ratings a
    nd vector cosines when a pair contains Either a diagnosis in the ICD 9 category,
    or treats or prevents a condition in that ICD 9 category
    """
    filename_to_embedding_matrix, idx_to_cui, cui_to_idx = generate_overlapping_sets_cui(filenames_type)
    
    
    print 'Number of overlapping cuis between embeddings: ' + str(len(idx_to_cui))
        
    #umnsrs_filename = 'UMNSRS_relatedness_mod458_word2vec.csv' #Can use judged relatedness instead of similiarity if desired
    umnsrs_filename = 'UMNSRS_similarity_mod449_word2vec.csv'
    with open(str(data_folder / umnsrs_filename), 'rU') as f:
        umnsrs_rows = f.readlines()[1:]
    
    # Find all cuis in the overlapping matrix
    cuis = cui_to_idx.keys()

    # Start arrays to store filenames, and to store the comparison values
    filename_all = []
    value_all = []
    
    # Following code to output which cuis are overlapping between all sets 
    ## o = open(str(results_folder / 'overlappingcuis.txt' ), 'w')
    ## for cui in cuis: o.write(str(cui)+'\n')
    ## o.close()

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
                    cos_sim = cosine_similarity([vec_1], [vec_2])[0,0]
                    
                    # Recrod scores
                    veccos_scores.append(cos_sim)
                    unmsrs_scores.append(float(umnsrs_rating))
                    
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
        
        # Carryout and append Spearman Rank Correlation     
        rho, pval = spearmanr(unmsrs_scores,veccos_scores)
        filename_all.append((filename))
        value_all.append(rho)
        files_processed += 1
   
    return filename_all, value_all, compares

#-----------------------------------------------------------------------------

def get_beam_bootstrap_by_systems(filenames_type, icd9_systems, cui_to_icd9_dicts, results):
    """JJN: Using Beam et al's boostrap method, calculates the percentage of known 
    drug-disease relations for a given system are within the top 5% of relations
    in a null distribution
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