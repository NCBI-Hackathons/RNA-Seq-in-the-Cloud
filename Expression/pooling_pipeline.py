#! /usr/bin/env python3

## pooling_pipeline.py - 

## index: A list of operations and functions included in this function
'''
0. import libraries and initialize global variables [Tommer]
1. parses input file [Tommer]
1a. create expression matrix for set of runs [Tommer]
2. single project PCA and elimination [Dan]
3. pool samples [Tommer]
3a. pool metadata [Tommer]
4. multiproject PCA/CCA [Dan]
4a. if (args.case): pool data and metadata again [Tommer]
5. elimination of run based on centroid distance [Tommer]
6. pass to DESeq2 for DGE analysis [Tommer]

'''

## 0. import libraries and initialize global variables, parse arguments

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
import argparse

gene_counts_path = 'gs://ncbi_sra_rnaseq/genecounts/'

print('parsing input file')

parser = argparse.ArgumentParser()
parser.add_argument('--case', action='store_true',
				help='for case control studies, perform PCA separately on cases and controls')
parser.add_argument('--input',type=str,
				help='metadata csv',default=None)
parser.add_argument('--name', type=str,
				help='name for whole analysis')
args = parser.parse_args()

joint_project = args.name
input_file_path = args.input
outdir = joint_project+'_output/'

subprocess.call('mkdir '+outdir, shell = True)

## 1. parses input file

total_df = pd.read_csv(input_file_path, sep = ',', index_col = 0)

projects = list(set(total_df['project']))


for project in projects:
	df = total_df[total_df['project'] == project]
	df = df.drop_duplicates(subset = 'Run')
	df.to_csv(outdir+project+'.tsv', sep = '\t')


## 1a. create expression matrix for a project (set of runs)



for project in projects:
	df = pd.read_csv(outdir+project+'.tsv', sep = '\t', index_col = 0)
	run_file1 = open(outdir+project+'_runs.txt', 'w')
	for run in list(df['Run']):
		run_file1.write(run+'\n')
	run_file1.close()
	genecount_df = pd.DataFrame()
	meta_df = pd.DataFrame()
	run_file = outdir+project+'_runs.txt'
	for run in open(run_file, 'r'):
		try:
			run = run.strip()
			df = pd.read_csv('gs://ncbi_sra_rnaseq/genecounts/'+run+'.genecounts', sep = '\t', index_col = 1, names = ['sample', 'gene', run])
			df = df[run]
			df.columns = [run]
			genecount_df = pd.concat([genecount_df, df], axis = 1)

		except FileNotFoundError:
			print(run)
			continue

	genecount_df = genecount_df.drop(['__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned', '__too_low_aQual'], axis = 0)
	genecount_df.to_csv(outdir+project+'_genecounts.txt', sep = '\t')


## 2. perform single project PCA

print('looking for outliers in individual projects')

if (args.case):
	print('separating cases from controls')
	new_projects = []
	for project in projects:
		meta_df = pd.read_csv(outdir+project+'.tsv', sep = '\t', index_col = 0)
		genecount_df = pd.read_csv(outdir+project+'_genecounts.txt', sep = '\t', index_col = 0)
		total_runs = list(set(list(meta_df['Run'])) & set(list(genecount_df)))
		meta_df = meta_df.loc[meta_df['Run'].isin(total_runs)]
		genecount_df = genecount_df[total_runs]
		case_meta_df = meta_df[meta_df['condition'] == 'case']
		control_meta_df = meta_df[meta_df['condition'] == 'control']
		case_genecounts_df = genecount_df[meta_df[meta_df['condition'] == 'case']['Run']]
		control_genecounts_df = genecount_df[meta_df[meta_df['condition'] == 'control']['Run']]
		print('project ID is: '+project)
		if (len(list(case_genecounts_df)) > 0):
			new_projects.append(project+'_case')
			case_meta_df.to_csv(outdir+project+'_case.tsv', sep = '\t')
			case_genecounts_df.to_csv(outdir+project+'_case_genecounts.txt', sep = '\t')
		if (len(list(control_genecounts_df)) > 0):
			new_projects.append(project+'_control')
			control_meta_df.to_csv(outdir+project+'_control.tsv', sep = '\t')
			control_genecounts_df.to_csv(outdir+project+'_control_genecounts.txt', sep = '\t')
	projects = new_projects


for project in projects:
	print(project)
	subprocess.call('python3 counts_pca-2.py --projname '+outdir+project+' --df '+outdir+project+'_genecounts.txt --md '+outdir+project+'.tsv --review', shell = True)

## 3. pool data for runs that successfully pass through single sample pca

print('pooling data from separate projects')

total_df = pd.DataFrame()
project_list = []

for project in projects:
	df = pd.read_csv(outdir+project+'_genecounts.txt', sep = '\t', index_col = 0)
	total_df = pd.concat([total_df, df], axis = 1)
	for i in range(0,len(list(df))):
		project_list.append(project)

total_df.to_csv(outdir+joint_project+'_genecounts.txt', sep = '\t')

## 3a. pool metadata

print('pooling metadata from separate projects')

joint_df = pd.read_csv(outdir+joint_project+'_genecounts.txt', sep = '\t')

print(len(list(joint_df)))
meta_df = pd.DataFrame()

for project in projects:
	try:
		sub_meta_df = pd.read_csv(outdir+project+'_metadata_outlierannot.txt', sep = '\t', index_col = 0)
		sub_meta_df = sub_meta_df[['type', 'condition', 'project', 'is_outlier']]
		sub_meta_df = sub_meta_df.dropna()
		meta_df = pd.concat([sub_meta_df, meta_df], axis = 0)
	except FileNotFoundError:
		projects.remove(project)


meta_df.columns = ['type', 'condition', 'project', 'is_outlier']

meta_df.to_csv(outdir+joint_project+'.tsv', sep = '\t')


## 4. multiproject PCA

print('looking for outliers in overall data')

projects = [joint_project]

if (args.case):
	new_projects = []
	for project in projects:
		meta_df = pd.read_csv(outdir+joint_project+'.tsv', sep = '\t', index_col = 0)
		genecount_df = pd.read_csv(outdir+joint_project+'_genecounts.txt', sep = '\t', index_col = 0)
		total_runs = list(set(list(meta_df.index)) & set(list(genecount_df)))
		meta_df = meta_df.ix[total_runs]
		genecount_df = genecount_df[total_runs]
		case_meta_df = meta_df[meta_df['condition'] == 'case']
		control_meta_df = meta_df[meta_df['condition'] == 'control']
		case_genecounts_df = genecount_df[meta_df[meta_df['condition'] == 'case'].index]
		control_genecounts_df = genecount_df[meta_df[meta_df['condition'] == 'control'].index]
		print(project)
		if (len(list(case_genecounts_df)) > 0):
			new_projects.append(project+'_case')
			case_meta_df.to_csv(outdir+project+'_case.tsv', sep = '\t')
			case_genecounts_df.to_csv(outdir+project+'_case_genecounts.txt', sep = '\t')
		if (len(list(control_genecounts_df)) > 0):
			new_projects.append(project+'_control')
			control_meta_df.to_csv(outdir+project+'_control.tsv', sep = '\t')
			control_genecounts_df.to_csv(outdir+project+'_control_genecounts.txt', sep = '\t')
	projects = new_projects

print(projects)


for project in projects:
	print('python3 counts_pca-2.py --projname '+outdir+project+' --df '+outdir+project+'_genecounts.txt --md '+outdir+project+'.tsv --review')
	subprocess.call('python3 counts_pca-2.py --projname '+outdir+project+' --df '+outdir+project+'_genecounts.txt --md '+outdir+project+'.tsv --review', shell = True)

## 4a. if (args.case): pool data and metadata again


if (args.case):
	print('pooling data and metadata')
	total_df = pd.DataFrame()
	total_meta_df = pd.DataFrame()
	for project in projects:
		condition_df = pd.read_csv(outdir+project+'_genecounts.txt', sep = '\t', index_col = 0)
		total_df = pd.concat([total_df, condition_df], axis = 1)
		meta_df = pd.read_csv(outdir+project+'_metadata_outlierannot.txt', sep = '\t', index_col = 0)
		total_meta_df = pd.concat([total_meta_df, meta_df], axis = 0)

	total_df.to_csv(outdir+joint_project+'_genecounts.txt', sep = '\t')
	total_meta_df.to_csv(outdir+joint_project+'_metadata_outlierannot.txt', sep = '\t')



## 5. provide metadata matrix to be passed to DESeq2

print('preparing output')

joint_df = pd.read_csv(outdir+joint_project+'_genecounts.txt', sep = '\t', index_col = 0)

meta_df = pd.DataFrame()


meta_df = pd.read_csv(outdir+joint_project+'_metadata_outlierannot.txt', sep = '\t', index_col = 0)

joint_meta_df = pd.read_csv(outdir+joint_project+'.tsv', sep = '\t', index_col = 0)
joint_meta_df = joint_meta_df[joint_meta_df['is_outlier'] == True]

meta_df.ix[joint_meta_df.index, 'is_outlier'] = True


meta_df.columns = ['type', 'condition', 'project', 'is_outlier']


total_runs = list(set(list(meta_df.index)) & set(list(joint_df)))

meta_df = meta_df.ix[total_runs]
meta_df = meta_df[~meta_df.index.duplicated(keep='first')]

joint_df = joint_df.loc[:,~joint_df.columns.duplicated(keep='first')]


outdir = outdir+'final_output/'
subprocess.call('mkdir '+outdir, shell = True)

meta_df.to_csv(outdir+joint_project+'_metadata.txt', sep = '\t')
joint_df.to_csv(outdir+joint_project+'_genecounts.txt', sep = '\t')


