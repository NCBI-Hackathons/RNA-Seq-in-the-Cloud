# Team Expression

![Alt text](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Expression/progress_report.png "Title")


## PCA_counts
Example:<br/>
python counts_pca.py -projname PROJNAME -df ERP000546_genecounts.txt --md meta_data.tsv --review

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

## DESeq2
Rough DESeq2 methods: Users can compare expression across multiple tissue types and conditions by running "Rscript deseq_contrast_groups.R --counts $counts.txt --a $attributes.txt --i $comparisons.txt --o $output_path". The counts.txt and attributes.txt files are output by "Dan and Tommer scripts".  A 'count' data frame contains read counts for each gene for each individual run. An 'attribute' data frame contains metadata that associates each run with its project ID, sample type, and condition. Sample type can describe tissue of origin, sex, age, etc. Additionally, there can be several sample type columns in each attribute file. The condition of each run is a binary categorical (case or control). The comparisons.txt file is user generated based on which comparisons the user finds interesting. The comparisons.txt file should be a tab delimited text file with three columns. Column one contains the name of a sample type column, column two and three contain the condition names. For example, a row might be "kidney" "cancer" "healthy". The user can define comparisons of interest based on sample types and conditions of each run. The model design would include sample types and conditions as factors (design = ~ projectID+sample type) and would contrast all runs for healthy kidney samples vs. cancer kidney samples. Output depth normalized counts for each run, log fold-change and p-values for each genes


