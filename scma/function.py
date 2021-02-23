#!/usr/bin/env python
"""
# Author: Kang Tian
# File Name: function.py
# Description:
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def _peaks(meta):
    """
    Find peaks of signals
    """
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(meta.iloc[:,0].values)
    return meta.iloc[peaks,:]

def _ppm(x1, x2):
    """
    Calculate ppm
    """
    return 1000000 * (x2-x1)/x1

def _merge_peaks(l):
    """
    Merge signals if the difference of Nuclear to cytoplasmic ratio is 1
    """
    idx = []
    while len(l)>0:
        first, *rest = l
        first = set(first)

        lf = -1
        while len(first)>lf:
            lf = len(first)
            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)     
            rest = rest2
        first = list(first)
        first.sort()
        idx.append(first)
        l = rest    
    return idx

def _merge(l,
           ppm_threshold=10):
    """
    Merge index if ppm<ppm_threshold
    """
    index = [0]
    i=0
    for j in np.arange(i+1, len(l)):
        ppm = _ppm(l[i], l[j])
        if ppm <= ppm_threshold:
            index.append(i)
        else:
            i=j
            index.append(i)
    return index

def _filter1(meta,
             ppm_threshold=5):
    """
    Filter signals by ppm
    """
    index = []
    for i in range(meta.shape[0]-1):
        for j in np.arange(i+1, meta.shape[0]):
            ppm = _ppm(meta.index[i]+1, meta.index[j])
            if abs(ppm) <= ppm_threshold:
                index.append([i,j])
            if ppm > ppm_threshold:
                index.append([i])
                break
    return index

def _filter2(l,
             decrease=True):
    """
    Filter signals if decrease
    """
    idx = []
    if decrease:
        for i, item in enumerate(l):
            if(i == 0):
                _max = item
                idx.append(i)
            else:
                if item > _max:
                    idx.append(i)
                _max = item
    else:
        idx = [0]
    return idx

def _filter3(meta, filter_peak=0.5):
    """
    Filter peaks based on number of cells.
    """
    meta = meta[(meta>0).sum(axis=1)>meta.shape[1]*filter_peak]
    return meta.T

def _distance(list1, list2):
    """
    """
    array = np.asarray(list2)
    pair = [list2[(np.abs(array - value)).argmin()] for value in list1]
    return pair


def _index(ppm_threshold=20, begin=300, end=1010):
    """
    """
    from math import ceil

    c = np.int(begin)
    end =  ceil(end)
    index = []
    while c <= end:
        index.append(c)
        c = c*(1 + (ppm_threshold/1e6))
    return index

def _process_cell(meta, 
                  benchmark,
                  peak=False,
                  ppm_threshold_peak=10,
                  decrease=True):
    """
    Pre-process each cell
    """
    from functools import reduce
    from math import ceil
    
    if not peak:
        meta = _peaks(meta)
    index = _filter1(meta, ppm_threshold=ppm_threshold_peak)
    index = _merge_peaks(index)
    
    index_ = []
    for item in index:
        idx = _filter2(meta.iloc[item,0].values, decrease=decrease)
        index_.append([item[i] for i in idx])
        
    index = list(set(reduce(lambda x,y:x+y, index_)))
    index.sort()
    meta = meta.iloc[index,:]
    
    merged_index = _merge(meta.index, ppm_threshold=ppm_threshold_peak)
    benchmark = _distance(meta.index[merged_index], benchmark)
    meta.loc[:,'benchmark'] = benchmark
    meta = meta.groupby(by='benchmark').max()
    return meta

def _process(df,
             groups=None,
             peak=False,
             ppm_threshold_peak=10,
             ppm_threshold_cell=20,
             decrease=True):
    """
    """
    df = df.dropna(axis=0,how='all')  
    df = df.dropna(axis=1,how='all')
    if groups is not None:
        print(groups)
        for group in groups:
            print(group)
            df.columns = df.columns.str.replace(group, group+'-')
    n_cells = np.int(df.shape[1]/2)
    cell_names = [name for i,name in enumerate(df.columns) if i%2==0]
    group_names = [item.split('-')[0] for item in cell_names]
    
    m_z = df.iloc[:,[i for i in range(n_cells) if i%2==0]]
    begin = m_z.min().min()
    end = m_z.max().max()
    benchmark = _index(ppm_threshold=ppm_threshold_cell,
                       begin=begin,
                       end=end)
    
    Meta = []
    for i in range(n_cells):
        meta = df.iloc[:,2*i:2*(i+1)]
        meta = meta.dropna(axis=0,how='all') 
        meta.index = meta.iloc[:,0]
        meta = pd.DataFrame(meta.iloc[:,1])
        Meta.append(_process_cell(meta, 
                                  benchmark=benchmark,
                                  peak=peak,
                                  ppm_threshold_peak=ppm_threshold_peak, 
                                  decrease=decrease))
    meta = pd.concat(Meta, axis=1)
    meta.columns = cell_names
    return meta

def _knn(df, n_neighbors=5):
    """
    """
    from sklearn.impute import KNNImputer
    imputer = KNNImputer(n_neighbors=n_neighbors, weights='uniform', metric='nan_euclidean')
    df_imputed = imputer.fit_transform(df)
    df_imputed = pd.DataFrame(df_imputed, index=df.index, columns=df.columns)
    return df_imputed

def _plot_embedding(df, labels, method='tSNE', cmap='tab20', figsize=(6, 6), markersize=50, show_legend=True,
                   return_emb=False, save=False, save_emb=False):
    """
    df: DataFrame nsamples x nfeatures
    labels: labels for each sample
    """
    df = df.fillna(0)
    if method == 'tSNE':
        from sklearn.manifold import TSNE
        X = TSNE(n_components=2, random_state=124).fit_transform(df)
    if method == 'UMAP':
        from umap import UMAP
        X = UMAP(n_neighbors=30, min_dist=0.1).fit_transform(df)
    if method == 'PCA':
        from sklearn.decomposition import PCA
        X = PCA(n_components=2, random_state=124).fit_transform(df)
    if method == 'PLS':
        from sklearn.cross_decomposition import PLSRegression
        from sklearn.preprocessing import LabelEncoder
        encode = LabelEncoder()
        ref = encode.fit_transform(labels)
        pls2 = PLSRegression(n_components=2)
        pls2.fit(df, ref)
        X = pls2.x_scores_
            
    plt.figure(figsize=figsize)

    if cmap is not None:
        cmap = cmap
    elif len(classes) <= 10:
        cmap = 'tab10'
    elif len(classes) <= 20:
        cmap = 'tab20'
    else:
        cmap = 'husl'
    palette = sns.color_palette(cmap, n_colors=len(np.unique(labels)))
        
    ax = sns.scatterplot(x=X[:,0],
                         y=X[:,1],
                         hue=labels,
                         palette=palette, 
                         marker="o",
                         legend='full',
                         s=markersize)

    ax.tick_params(axis='x', bottom=True, top=False, labeltop=False, labelbottom=True, labelsize=12, length=3, pad=3)
    ax.tick_params(axis='y', left=True, right=False, labelright=False, labelleft=True, labelsize=12, length=3, pad=3)

    ax.set_xlabel('{}_1'.format(method), fontsize=15, labelpad=10, va='center')
    ax.set_ylabel('{}_2'.format(method), rotation=90, fontsize=16, labelpad=10, va='center')

    if save:
        plt.savefig(save, format='pdf', bbox_inches='tight')
    else:
        plt.show()
        
    if save_emb:
        np.savetxt(save_emb, X)
    if return_emb:
        return X
    
def _split(df, labels):
    """
    """
    samples = []
    label_ = np.unique(labels)
    for label in label_:
        samples.append(df.iloc[np.where(np.array(labels)==label)[0],:].values)
    return samples

def _foldchange(df, labels):
    """
    """
    DE = {}
    label_ = np.unique(labels)
    for label in label_:
        df1 = df.iloc[np.where(np.array(labels)==label)[0],:].values
        df2 = df.iloc[np.where(np.array(labels)!=label)[0],:].values
        DE[label] = np.log2(df1.mean(axis=0)/df2.mean(axis=0))
        DE = pd.DataFrame(DE)
        DE.index = df.columns
    return DE

def _varAnalysis(df, labels):
    """
    """
    from scipy.stats import levene
    if df.shape[0] != len(labels):
        raise ValueError("The number of input samples is not equal to labels size")
        return 0
    
    label_ = np.unique(labels)
    groups = _split(df, labels)
    if len(label_) == 2:
        print('Performing t-test analysis...')
        from scipy.stats import ttest_ind
        F, P =[], []
        for i in range(df.shape[1]):
            sample = [item[:,i] for item in groups]
            stat, p = levene(*sample)
            if p < 0.05:
                f, p = ttest_ind(*sample, equal_var=False)
            else:
                f, p = ttest_ind(*sample, equal_var=True)
            F.append(f)
            P.append(p)

    elif len(label_) > 2:
        print('Performing anova analysis...')
        F, P =[], []
        for i in range(df.shape[1]):
            sample = [item[:,i] for item in groups]
            stat, p = levene(*sample)
            if p < 0.05:
                from pingouin import welch_anova
                meta = pd.DataFrame(df.iloc[:,i])
                meta.columns = ['feature']
                meta['labels'] = labels
                result = welch_anova(data=meta, dv='feature', between='labels')
                f = result['F'].values[0]
                p = result['p-unc'].values[0]
            else:
                from scipy.stats import f_oneway
                f, p = f_oneway(*sample)
            F.append(f)
            P.append(p)
    else:
        raise ValueError("Groups for comparison are less than 2!")
    F = pd.DataFrame(F)
    P = pd.DataFrame(P)
    F.index = df.columns
    P.index = df.columns
    return F, P

def _diff(df, labels, p_value=0.05, log2fold=0.5):
    """
    """
    label_ = np.unique(labels)
    de = _foldchange(df, labels)
    F, P = _varAnalysis(df, labels)
    ind1 = np.where(de>log2fold)[0]
    ind2 = np.where(P<p_value)[0]
    index = np.array(list(set(ind1).intersection(set(ind2))))
    index.sort()
    result = pd.concat([de, P], axis=1)
    result.columns = list(label_) + ['p_values']
    return result.iloc[index,:]