# -*- coding: utf-8 -*-
"""
Elegant Scipy

@author: Sarvenaz
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
from collections import defaultdict
plt.style.use('style/elegant.mplstyle')
from scipy import stats

# %% Functions

def rmkm(counts, length):
    '''
    Calculate reads per kilobase transcript per million reads
    '''
    N = np.sum(counts, axis = 0)
    
    normed = 1e9 * counts / (N[np.newaxis, :] * length[:, np.newaxis])
    
    return normed

def reduce_xaxis_labels(ax, factor):
    ''' Show every ith label to prevent crowding on x-axis'''
    plt.setp(ax.xaxis.get_ticklabels(), visible = False)
    for label in ax.xaxis.get_ticklabels()[factor-1::factor]:
        label.set_visible(True)

def class_boxplot(data, classes, color=None, **kwargs):
    ''' Make boxplots with boxes colored according to their class'''
    all_classes = sorted(set(classes))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    class2color = dict(zip(all_classes, it.cycle(colors)))
    
    # map classes to data vector
    class2data = defaultdict(list)
    for distrib, cls in zip(data, classes):
        for c in all_classes:
            class2data[c].append([])
        class2data[cls][-1] = distrib
    
    # Plotting each class with appropriate color
    fig, ax = plt.subplots()
    lines = []
    for cls in all_classes:
        for key in ['boxprops', 'whiskerprops', 'flierprops']:
            kwargs.setdefault(key, {}).update(color = class2color[cls])
        box = ax.boxplot(class2data[cls], **kwargs)
        lines.append(box['whiskers'][0])
    ax.legend(lines, all_classes)
    
    return ax

# %%

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
    fig, ax = plt.subplots()
    ax.plot(x, density(x))
    ax.set_xlabel('Total counts per individual')
    ax.set_ylabel('Density')    
    
    # %%
    sample_index = np.random.choice(range(counts.shape[1]), size = 70, replace = False)
    counts_subset = counts[:, sample_index]
    
    fig, ax = plt.subplots(figsize = (4.8, 2.4))
    
    with plt.style.context('style/thinner.mplstyle'):
        ax.boxplot(counts_subset)
        ax.set_xlabel('Individuals')
        ax.set_ylabel('Gene Expression Counts')
        reduce_xaxis_labels(ax, 5)
        
    # %% Since above figure was too clustered around 0, lets look at he log n of the data
    fig, ax = plt.subplots(figsize = (4.8, 2.4))
    
    with plt.style.context('style/thinner.mplstyle'):
        ax.boxplot(np.log(counts_subset + 1))
        ax.set_xlabel('Individuals')
        ax.set_ylabel('log gene expression counts')
        reduce_xaxis_labels(ax, 5)
    # %% Let's normalize the data
    counts_lib_norm = (counts / total_counts) * 1e6
    counts_subset_lib_norm = counts_lib_norm[:, sample_index]
    
    fig, ax = plt.subplots(figsize = (4.8, 2.4))
    
    with plt.style.context('style/thinner.mplstyle'):
        ax.boxplot(np.log(counts_subset_lib_norm + 1))
        ax.set_xlabel('Individuals')
        ax.set_ylabel('log gene expression counts')
        reduce_xaxis_labels(ax, 5)
    
    # %% 
    log_counts_3 = list(np.log(counts.T[:3]+1))
    log_ncounts_3 = list(np.log(counts_lib_norm.T[:3]+1))
    ax = class_boxplot(log_counts_3 + log_ncounts_3, 
                       ['raw counts']*3+['normalized by library size']*3,
                       labels = [1,2,3,1,2,3])
    ax.set_xlabel('Sample Num')
    ax.set_ylabel('log gene expression counts')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    