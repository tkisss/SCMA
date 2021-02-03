#!/usr/bin/env python
"""
# Author: Kang Tian
# File Name: run.py
# Description:
"""
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import os
import time
import matplotlib.pyplot as plt
from sklearn.preprocessing import MaxAbsScaler

from scma.function import _process, _filter3, _knn, _plot_embedding, _diff

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Single cell metabolomics analysis')
    parser.add_argument('--file', '-f', type=str, default='')
    parser.add_argument('--ppm_threshold_peak', type=int, default=10)
    parser.add_argument('--ppm_threshold_cell', type=int, default=20)
    parser.add_argument('--decrease', default=True)
    parser.add_argument('--peak', default=False)
    parser.add_argument('--filter', type=float, default=0.5)
    parser.add_argument('--knn', action='store_true')
    parser.add_argument('--n_neighbors', type=int, default=5)
    parser.add_argument('--method', '-m', type=str, default='PLS')
    parser.add_argument('--p_value', type=float, default=0.05)
    parser.add_argument('--log2fold', type=float, default=0.5)
    parser.add_argument('--seed', type=int, default=2020)
    parser.add_argument('--outdir', '-o', type=str, default='./')
    args = parser.parse_args()

    np.random.seed(args.seed) 
    
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    
    timestart = time.time()
    if os.path.isfile(args.file):
        df = pd.read_csv(args.file)
    else:
        raise ValueError("{} file does not exist!")
    print('Preprocessing...')
    meta = _process(df,
                    peak=args.peak,
                    ppm_threshold_peak=args.ppm_threshold_peak,
                    ppm_threshold_cell=args.ppm_threshold_cell,
                    decrease=args.decrease)
    meta.to_csv(outdir+'/processed.csv')
    
    plt.figure(figsize=(5,5))
    sns.violinplot((meta>0).sum(axis=1), orient='h')
    plt.savefig(outdir+"/violinplot.png")
    
    print('Filtering peaks...')
    meta = _filter3(meta, filter_peak=args.filter)
    meta.to_csv(outdir+'/filtered.csv')
    
    if not args.knn:
        print('Performing knn for imputation...')
        meta = _knn(meta,n_neighbors=args.n_neighbors)
        meta.to_csv(outdir+'/knn.csv')
    
    labels = [item.rsplit('-',1)[0] for item in meta.index]
    
    print('Performing {} for visualization...'.format(args.method))
    emb = _plot_embedding(meta, method=args.method, labels=labels, save=outdir+"/embedding.pdf", return_emb=True)
    emb = pd.DataFrame(emb, index=meta.index, columns=["{}_1".format(args.method),"{}_2".format(args.method)])
    emb.to_csv(outdir+"/{}.txt".format(args.method))
    
    print('Performing differential expression analysis...')
    de = _diff(meta, labels, p_value=args.p_value, log2fold=args.log2fold)
    de.to_csv(outdir+'/de.csv')
    
    meta = pd.DataFrame(MaxAbsScaler().fit_transform(meta), index=meta.index, columns=meta.columns)
    sns.clustermap(meta.loc[:,de.index].T,
                   col_cluster=False,cmap='RdBu_r',
                   figsize=(meta.shape[0]/15, de.shape[0]/8),
                   dendrogram_ratio=0.25)
    plt.savefig(outdir+"/heatmap.png")
    
    
    timeend = time.time()
    print('Running time: {}'.format(timeend-timestart))
    print('Done')
    