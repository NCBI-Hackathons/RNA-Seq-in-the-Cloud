## CreateDataFrames.py - Create DataFrames of gene_counts, other operations 

## index: A list of operations and functions included in this function
'''
0. import libraries and initialize global variables
1. create dataframe
2. create expression matrix for set of runs

'''

## 0. import libraries and initialize global variables

import sys
import sklearn
import pandas as pd
import gcsfs
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math
import random
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

gene_counts_path = 'gs://ncbi_sra_rnaseq/genecounts/'


## 1. create dataframe
'''
df1 = pd.read_csv('gs://ncbi_sra_rnaseq/genecounts/ERR030875.genecounts', sep = '\t', index_col = 1, names = ['sample', 'gene', 'genecounts1'])
df2 = pd.read_csv('gs://ncbi_sra_rnaseq/genecounts/ERR315471.genecounts', sep = '\t', index_col = 1, names = ['sample', 'gene', 'genecounts2'])

df = pd.concat([df1, df2], axis = 1)
print(df)

t = sns.scatterplot(x = np.log(df['genecounts1']), y = np.log(df['genecounts2']))
print(df['genecounts1'].corr(df['genecounts2']))
print(df['genecounts1'].corr(df['genecounts2'], method = 'spearman'))

plt.tight_layout()
t.figure.savefig('../results/test.png')

'''

## 2. create expression matrix for set of runs
'''
run_file = '../data/ERP000546_runs.txt'


genecount_df = pd.DataFrame()
meta_df = pd.DataFrame()

for run in open(run_file, 'r'):
	try:
		run = run.strip()
		df = pd.read_csv('gs://ncbi_sra_rnaseq/genecounts/'+run+'.genecounts', sep = '\t', index_col = 1, names = ['sample', 'gene', run])
		df = df[run]
		df.columns = [run]
		genecount_df = pd.concat([genecount_df, df], axis = 1)
		
		meta_df.ix[run, 'type'] = random.randint(0,1) ## Make this A/B

	except FileNotFoundError:
		continue

genecount_df.to_csv('../data/ERP000546_genecounts.txt', sep = '\t')
meta_df.to_csv('../data/ERP000546_attributes.txt', sep = '\t')
'''



