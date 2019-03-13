# Team Expression

![Alt text](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Expression/progress_report.png "Title")


## PCA_counts
Example:<br/>
python3 counts_pca.py --projname *PROJNAME* --df *ERP000546_genecounts.txt* --md *meta_data.tsv* --review

### --projname (str)
Name of SRA project 

### --df (str)
Path to dataframe containing gene counts for all samples within the project

### --md (str)
Path to metadata file for all samples in project

### --review (flag)
Either keep outliers and update metadata file, or simply remove

### OUTPUT

1. Dataframe of counts
2. Dataframe of metadata

## Pooling
Intro

In order to maximize sample size to perform a differential expression analysis, we aim to pool runs together from different projects. However, this is only possible for runs that share similar properties between them. From general metadata, we are able to begin to gather clues as to which runs are able to be pooled together, but this does not provide conclusive support. We suggest performing several dimensionality reduction techniques in order to increase the number of runs that can be included in a differential gene expression analysis.

Methods

PCA → Eliminate 2+ sd away from centroid → CCA

Initially, we perform a principal component analysis (PCA) to identify linear combinations of gene measurements that maximize the variation between runs within a project. This enables us to identify and eliminate runs that are outliers within their own project, or identify and flag distinct biological variation in our data. For case/control studies, this is done by separating each of the case runs from control runs, and performing our PCA. For projects without a distinct case and control, there is no need to separate the data.

We then measure the euclidean distance on our principal component space to the centroid. Runs, that when plotted in this space, that have a distance > 2 standard deviations from our centroid are considered outliers, and eliminated from this project, and flagged for further analysis.

After eliminating outlier runs from multiple projects that we aim to pool, we perform a canonical correlation analysis (CCA), which captures shared signal across multiple datasets. We do this in order to distinguish biological differences from technical sources of variation


## Differential expression comparisons with DESeq2
### deseq2_contrast_groups.R
Example Usage:
Rscript deseq2_contrast_groups.R  geneCounts.txt --attributes metaData.txt --type blood --outliers --outdir output/ 

Users can compare expression across multiple tissue types and conditions by running Rscript deseq2_contrast_groups.R.
The DESEq model design includes include sample types (tissue, age, sex) and project IDs as factors (design = ~ projectID+sample type).
Users will define case vs. control comparisons to analyze with the --type flag by specifying a variable in the 'type' column of metaData.txt. This script will output 'control_vs_case.txt' containing log fold changes in expression between cases and controls matching matching the --type option. The output also contains p-values for each gene and depth normalized read counts for each gene for each sample. Only genes with at least two read counts across all samples are used for analysis.

Flags:
--counts
nxm table where n=genes and m=counts for each run

--attributes
nxm table where n=runs and m=meta data attributes (type, condition, project, outlier)
type is usually a tissue type (but can also be age, sex, etc.)
condition is either case or control (usually meaning healthy or disease samples)
project is the SRA project 
outlier is True or False and determined upstream in the pipeline

--type
MUST be defined!
MUST correspond to a variable in the 'type' column of the attribute table
For example, --type liver

--outliers
Include this flag if you want to perform analysis including outliers determined upstream in the pipeline
\n
--outdir
path to output directory that will be created
default is the working directory



