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
project = 'ERP003613'
run_file = '../data/'+project+'_runs.txt'


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

genecount_df = genecount_df.drop(['__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned', '__too_low_aQual'], axis = 0)
genecount_df.to_csv('../data/'+project+'_genecounts.txt', sep = '\t')
meta_df.to_csv('../data/'+project+'_attributes.txt', sep = '\t')
'''
## 3. PCA for a given count matrix
'''
df = pd.read_csv('../data/'+project+'_genecounts.txt', sep = '\t', index_col = 0)

df = df.transpose()
df = StandardScaler().fit_transform(df)
pca = PCA(n_components = 5)
pc = pca.fit_transform(df)

pc_df = pd.DataFrame(pc, columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
print(pc_df)


t = sns.scatterplot(x = pc_df['PC1'], y = pc_df['PC2'])
t.figure.savefig('pca.png')
plt.close()
'''
## 4. PCA with runs from multiple projects

projects = ['ERP003613', 'ERP000546']

total_df = pd.DataFrame()
project_list = []

for project in projects:
	df = pd.read_csv('../data/'+project+'_genecounts.txt', sep = '\t', index_col = 0)
	total_df = pd.concat([total_df, df], axis = 1)
	for i in range(0,len(list(df))):
		project_list.append(project)

total_df = total_df.transpose()

total_df = StandardScaler().fit_transform(total_df)
pca = PCA(n_components = 5)
pc = pca.fit_transform(total_df)

pc_df = pd.DataFrame(pc, columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
pc_df['project'] = project_list

print(pc_df)

t = sns.scatterplot(x = pc_df['PC1'], y = pc_df['PC2'], hue = pc_df['project'])
t.figure.savefig('multi_project_pca.png')
plt.close()




