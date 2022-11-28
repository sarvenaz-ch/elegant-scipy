# -*- coding: utf-8 -*-

"""
Created on Mon Nov 14 14:41:03 2022

@author: sarve
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage # hierchical clustering
plt.style.use('style/elegant.mplstyle')


def quantile_norm(x):
    ''' Normalize the columns of x, so each have the same distribution'''
    quantiles = np.mean(np.sort(x, axis = 0), axis = 1) # sort along column, take mean along each row
    ranks = np.apply_along_axis(stats.rankdata, 0, x)
    
    rank_indecis = ranks.astype(int)-1
    
    Xn = quantiles[rank_indecis]
    
    return Xn


def quantile_norm_log(X):
    ''' log transformation of input data before using the quantile_norm'''
    logX = np.log(X+1)
    logXn = quantile_norm(logX)
    return logXn


def plot_col_density(data):
    '''density plot for each column'''
    density_per_col = [stats.gaussian_kde(col) for col in data.T]
    x = np.linspace(np.min(data), np.max(data), 100) # x-axis
    
    _, ax = plt.subplots()
    for density in density_per_col:
        ax.plot(x, density(x))
    ax.set_xlabel('Data values (per column)')
    ax.set_ylabel('Density')
    

def most_variable_rows(data, *, n= 1000):
    ''' Subset data to the n most variable rows'''
    # compute variance along the columns
    rowvar = np.var(data, axis = 1)
    sort_indecis = np.argsort(rowvar)[-n:]
    
    return data[sort_indecis, :]


def bicluster(data, linkage_method = 'average', distance_metric = 'correlation'):
    ''' Clusters the rows and the columns of a matrix'''
    y_rows = linkage(data, method = linkage_method, metric = distance_metric) #clustering of the rows of the input data
    y_cols = linkage(data.T, method = linkage_method, metric = distance_metric) # clustering of the columns of the input data
    
    return y_rows, y_cols
# %%

if __name__ == "__main__":
    filename = 'data/counts.txt'
    data_table = pd.read_csv(filename, index_col = 0)
    
    # print(data_table.head())
    
    counts = data_table.values # dataframe to 2D array
    
    ''' PLOTTING THE DATA '''
    log_counts = np.log(counts+1)
    plot_col_density(log_counts)
    
    # %%
    ''' Normalizing the input using quantile'''
    log_counts_normalized = quantile_norm_log(counts)
    plot_col_density(log_counts_normalized)
    
    
    #%%
    '''________________   CLUSTERING ___________________________________'''
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    