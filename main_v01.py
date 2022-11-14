# -*- coding: utf-8 -*-
"""
Elegant Scipy

@author: Sarvenaz
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('style/elegant.mplstyle')
from scipy import stats

# %%

def rmkm(counts, length):
    '''
    Calculate reads per kilobase transcript per million reads
    '''
    N = np.sum(counts, axis = 0)
    
    normed = 1e9 * counts / (N[np.newaxis, :] * length[:, np.newaxis])
    
    return normed

if __name__ == '__main__':
    fldr_curr = os.path.dirname(os.path.realpath(__file__))
    file_counts = fldr_curr+'\data\counts.txt'
    
    with open(file_counts, 'rt') as f:
        data_table = pd.read_csv(f, index_col = 0)
        
    samples = list(data_table.columns)
    
    file_genes = fldr_curr+'\data\genes.csv'
    gene_info = pd.read_csv(file_genes)
    
    
    print(f'loaded the right data? {data_table.shape[0]==gene_info.shape[0]}')
    
    # %%
    # Lets match the two dataset we have
    gene_info.index = gene_info.GeneSymbol
    matched_index = pd.Index.intersection(data_table.index, gene_info.index)
    
    counts = np.asarray(data_table.loc[matched_index], dtype = int)
    gene_names = np.array(matched_index)
    
    gene_length = np.asarray(gene_info.loc[matched_index]['GeneLength'], dtype = int)
    
    total_counts = np.sum(counts, axis = 0)
    
    # Use Gaussian smoothing to estimate the density
    density = stats.kde.gaussian_kde(total_counts)
    
    # Make values for ehivh to estimate the density for plotting
    x = np.arange(min(total_counts), max(total_counts), 10000)
    
    # display
    fig, ax = plt.subplot()
    ax.plt(x, density(x))
    ax.set_xlabel('Total counts per individual')
    ax.set_ylabel('Density')    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    