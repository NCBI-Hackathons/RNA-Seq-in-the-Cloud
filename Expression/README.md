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



