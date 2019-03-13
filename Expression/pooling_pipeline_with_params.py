#! /usr/bin/env python3

## pooling_pipeline.py - 

## index: A list of operations and functions included in this function
'''
0. import libraries and initialize global variables
1. create expression matrix for set of runs [Tommer]
2. single project PCA and elimination [Dan]
3. pool samples
4. multiproject PCA/CCA
5. elimination of run based on centroid distance
6. pass to DESeq2 for DGE analysis

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
import subprocess
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline


parser = argparse.ArgumentParser(description='Pooling pipeline')
parser.add_argument('--genecountspath',type=str,help='path to gene counts', default='gs://ncbi_sra_rnaseq/genecounts/')
parser.add_argument('--df',type=str,
                    help='sample gene counts dataframe, tab separated')
parser.add_argument('--metadatapath',type=str,
                    help='metadata path',default='/home/wagner/')
parser.add_argument('--outpath', type=str,
					help='output path', default ='../data')

parser.add_argument('--projects', nargs='2', default=['SRP062966_nblood', 'SRP071965_nblood'], help='The list of projects 2 required')

parser.add_argument('--review', action='store_true',
                help='retain outliers for analysis')

args = parser.parse_args()

out_path=args.outpath
gene_counts_path = args.genecountspath
#projects = ['ERP003613', 'ERP000546']
projects = args.projects
meta_path = args.metadatapath


## 1. create expression matrix for a project (set of runs)



for project in projects:
	df = pd.read_csv(meta_path + project+'.tsv', sep = '\t', index_col = 0)
	run_file1 = open(out_path + project+'_runs.txt', 'w')
	for run in list(df['Run']):
		run_file1.write(run+'\n')
	run_file1.close()
	genecount_df = pd.DataFrame()
	meta_df = pd.DataFrame()
	run_file = out_path + project + '_runs.txt'
	for run in open(run_file, 'r'):
		try:
			run = run.strip()
			df = pd.read_csv(gene_counts_path + run + '.genecounts', sep = '\t', index_col = 1, names = ['sample', 'gene', run])
			df = df[run]
			df.columns = [run]
			genecount_df = pd.concat([genecount_df, df], axis = 1)

		except FileNotFoundError:
			continue

	genecount_df = genecount_df.drop(['__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned', '__too_low_aQual'], axis = 0)
	genecount_df.to_csv(out_path + project + '_genecounts.txt', sep = '\t')


## 2. perform single project PCA

for project in projects:
	subprocess.call('python3 counts_pca.py --projname '+project+' --df '+ out_path + project + '_genecounts.txt --md ' + meta_path + project+'.tsv --review', shell = True)

## 3. pool data for runs that successfully pass through single sample pca

total_df = pd.DataFrame()
project_list = []

for project in projects:
	df = pd.read_csv(project+'_genecounts.txt', sep = '\t', index_col = 0)
	total_df = pd.concat([total_df, df], axis = 1)
	for i in range(0,len(list(df))):
		project_list.append(project)

total_df.to_csv(out_path + 'joint_genecounts.txt', sep = '\t')

## 3a. pool metadata

joint_df = pd.read_csv('joint_genecounts.txt', sep = '\t')

print(len(list(joint_df)))
meta_df = pd.DataFrame()

for project in projects:
	sub_meta_df = pd.read_csv(project + '_metadata_outlierannot.txt', sep = '\t', index_col = 0)
	sub_meta_df = sub_meta_df[['tissue', 'case','SRA_Study', 'is_outlier']]
	sub_meta_df = sub_meta_df.dropna()
	meta_df = pd.concat([sub_meta_df, meta_df], axis = 0)
	print(meta_df)


#meta_df = meta_df[['tissue', 'SRA_Study', 'is_outlier']]
meta_df.columns = ['type', 'condition', 'project', 'is_outlier']
meta_df['type'] = 'blood'

meta_df.to_csv(out_path + 'joint.tsv', sep = '\t')


## 4. multiproject PCA

for project in ['joint']:
	print('python3 counts_pca.py --projname ' + project + ' --df '+project+'_genecounts.txt --md '+project+'.tsv --review')
	subprocess.call('python3 counts_pca.py --projname ' + project + ' --df ' + project + '_genecounts.txt --md ' + project+'.tsv --review', shell = True)


## 5. provide metadata matrix to be passed to DESeq2

joint_df = pd.read_csv(out_path + 'joint_genecounts.txt', sep = '\t')

meta_df = pd.DataFrame()


meta_df = pd.read_csv('joint_metadata_outlierannot.txt', sep = '\t', index_col = 0)

joint_meta_df = pd.read_csv('joint.tsv', sep = '\t', index_col = 0)
joint_meta_df = joint_meta_df[joint_meta_df['is_outlier'] == True]

meta_df.ix[joint_meta_df.index, 'is_outlier'] = True


#meta_df = meta_df[['tissue', 'SRA_Study', 'is_outlier']]
meta_df.columns = ['type', 'condition', 'project', 'is_outlier']
meta_df['type'] = 'blood'

meta_df.to_csv(out_path + 'joint_metadata.txt', sep = '\t')

