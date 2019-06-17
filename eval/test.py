# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 19:59:20 2019

@author: jjnun
"""

import numpy as np
from scipy.spatial.distance import cdist

num_of_neighbor = 2
randy = np.random.rand(4,5)
onesy =  np.random.rand(4,5)

print 'onesy'
print onesy
print 'randy'
print randy
print 'cdist off shelf'
print cdist(randy,onesy,'cosine')

loopy = np.zeros((4,4))

r, c = randy.shape
ranks = np.zeros((r,num_of_neighbor + 1))


print 'cdist per line, 1'
print cdist(randy[0].reshape((1,c)),onesy, 'cosine')[0]
print 'just one line'
print loopy[1]

for i in range(r):
    #row = cdist(randy[i].reshape((1,c)),onesy,'cosine')
    #print row
    #print onesy.shape
    #print randy[i].reshape((1,c)).shape
    loopy[i] = cdist(randy[i].reshape((1,c)),onesy,'cosine')[0]
    row_ranks = np.argsort(loopy[i])[0:num_of_neighbor + 1]
    ranks[i] = row_ranks
    
print loopy

print loopy == cdist(randy,onesy,'cosine')

print ranks