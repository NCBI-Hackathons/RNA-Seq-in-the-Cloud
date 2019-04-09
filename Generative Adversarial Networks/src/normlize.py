import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--input_file", type=str,required=True)
parser.add_argument("--output_file", type=str,required=True)
parser.add_argument("--case_control", type=int, required=True)
args = parser.parse_args()
case_f = pd.read_csv(args.input_file, header=None)
case_f = case_f.applymap(lambda x: np.log(x + 1))
mean = np.mean(case_f.values)
std = np.std(case_f.values)
normalize_case_df = case_f.applymap(lambda x: (x-mean)/std)
normalize_case_df["case_control"] = args.case_control
normalize_case_df.to_csv(args.output_file,index=False)
