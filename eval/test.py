# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 22:57:08 2019

@author: jjnun
"""

import numpy as np
x = np.arange(10)
print("Original array:")
print(x)
np.random.shuffle(x)
n = 5
cat = ['hi','h','a','b','q','q','w','t','y','u']
cat2 = [1,2,3,4,5,6,5,5,5,5,5,5,5]
cat3 = np.zeros(15)
print cat3[np.argsort(x)[-n:]]

print [cat2[x] for x in np.argsort(x)[-n:]]

print 'lets test normazlie'

a = np.zeros((3, 3))
a[0,2] = 1
b = a[0,:]
print normalize(b.reshape(1,-1))[0]