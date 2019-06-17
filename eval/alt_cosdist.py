# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 19:29:34 2019

@author: jjnun
"""
import numpy as np

def cosine_vectorized_v2(array1, array2):
    sumyy = np.einsum('ij,ij->i',array2,array2)
    sumxx = np.einsum('ij,ij->i',array1,array1)[:,None]
    sumxy = array1.dot(array2.T)
    return (sumxy/np.sqrt(sumxx))/np.sqrt(sumyy)



#Y = cdist(icd9_embeddings, icd9_embeddings, 'cosine')

#n = len(icd9_embeddings)
#Y = np.zeros((n,n))

#for i in range(n):
#    Y[i,:] = cdist(icd9_embeddings)
    
    
#    (14988L, 500L)
#    (109052L, 500L)
    # (2110L, 200L) Whem other files