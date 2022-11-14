# -*- coding: utf-8 -*-
"""
Elegant Scipy

@author: Sarvenaz
"""
import os
import numpy as np
import pandas as pd

def rmkm(counts, length):
    '''
    Calculate reads per kilobase transcript per million reads
    '''
    N = np.sum(counts, axis = 0)
    
    normed = 1e9 * counts / (N[np.newaxis, :] * length[:, np.newaxis])
    
    return normed

if __name__ == '__main__':
    fldr_curr = os.path.dirname(os.path.realpath(__file__))
    filename = fldr_curr+'\data\counts.txt'
    
    with open(filename, 'rt') as f:
        data_table = pd.read_csv(f, index_col = 0)
        
    samples = list(data_table.columns)