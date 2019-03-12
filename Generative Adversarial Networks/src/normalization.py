#!/usr/bin/python3

#import libraries
import pandas as pd
import numpy as np
import os

#Get controls and add a column with 0.
os.chdir(/home/balan.ramesh/normal_skin)
my_file = pd.read_csv("counts.csv")
max_value = np.log(my_file.values.max())
normalize_df = my_file.applymap(lambda x: (np.log(x+1)/max_value))
normalize_df["case_control"] = 0
normalize_df.to_csv("control.csv",index=False)

#Get case and put 1 that column.
os.chdir(/home/balan.ramesh/melanoma)
my_file = pd.read_csv("counts.csv")
max_value = np.log(my_file.values.max())
normalize_df = my_file.applymap(lambda x: (np.log(x+1)/max_value))
normalize_df["case_control"] = 1
normalize_df.to_csv("case.csv",index=False)