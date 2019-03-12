#! /usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from sklearn.decomposition import PCA

def centeroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    return sum_x/length, sum_y/length

parser = argparse.ArgumentParser(description='PCA of gene counts')
parser.add_argument('--projname',type=str,help='name of project')
parser.add_argument('--df',type=str,
                    help='sample gene counts dataframe, tab separated')
parser.add_argument('--md',type=str,
                    help='metadata tsv',default=None)
parser.add_argument('--review', action='store_true',
                help='retain outliers for analysis')

args = parser.parse_args()

data = pd.read_csv(args.df, sep='\t', index_col=0)

meta_data = pd.read_csv(args.md,sep='\t',index_col='Run')

if len(data.columns) > 2:
    data_normcounts = data/data.sum()
    log_normcounts = np.log2(data_normcounts+1)
    df = log_normcounts.loc[[i for i in list(log_normcounts.index) if '_' not in i]]

    pca_model = PCA(n_components=10)
    df = df.T
    fit = pca_model.fit(df)

    explained_var = np.cumsum(fit.explained_variance_ratio_)

    pca_two = PCA(2)
    projected = pca_two.fit_transform(df)

    centroid = centeroid(projected)
    dist = [np.linalg.norm(i-centroid) for i in projected]

    dist_std,dist_mean,dist = np.std(dist),np.mean(dist),np.array(dist)
    ingroup = np.nonzero(~(dist > (dist_mean + 2*dist_std)))
    outliers = np.nonzero((dist > (dist_mean + 2*dist_std)))
    outlier_removed_df = data[data.columns[ingroup]]

    meta_data['is_outlier'] = np.nan
    is_out = {}
    for row in meta_data.index:
        if row in data.columns[outliers]:
            is_out[row] = True
        else:
            is_out[row] = False
    meta_data['is_outlier'] = pd.Series(is_out)
    if args.review:
        meta_data.to_csv('./{}_outliersannot.csv'.format(args.projname))
        data.to_csv('./{}_alldatakept.csv'.format(args.projname))
    else:
        meta_data = meta_data[meta_data['is_outlier']!= True]
        outlier_removed_df.to_csv('./{}_outliersremoved.csv'.format(args.projname))
    print('Found {} outliers'.format(len(outliers[0])))
