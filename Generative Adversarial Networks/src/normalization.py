#!/usr/bin/python3

#Import the libraries
import argparse
import pandas as pd
import numpy as np

'''
The case and the control file must be a csv (comma-separated) file.
Each gene has its own column and the rows contain each sample.
The value in each cell is a raw read count.
'''

def main():
	print("====================================================================\n")
	print("This program normalizes each count by the following formula:\n log(x+1)/log(max_value in the table + 1).\n The +1 adds 1 to all read counts so that log(0)=infinity is avoided.\n")
	print("====================================================================")
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument("--case", help="Specify the path/filename for the case file",type=str,required=True)
	parser.add_argument("--control", help="Specify the path/filename for the control file",type=str,required=True)
	parser.add_argument('-h', '--help',action='help',help="The case and the control file must be a csv (comma-separated) file.\nEach gene has its own column and the rows contain each sample.\nThe value in each cell is a raw read count.\n")
	args = parser.parse_args()
	if args.case:
		print("Your case file is ", args.case)
		case_f = pd.read_csv(args.case, header=None)
		max_value = np.log(case_f.values.max())
		normalize_case_df = case_f.applymap(lambda x: (np.log(x+1)/max_value))
		normalize_case_df["case_control"] = 0
		normalize_case_df.to_csv("normalized_case.csv",index=False)

	if args.control:
		print("Your control file is ", args.control)
		control_f = pd.read_csv(args.control, header=None)
		max_value = np.log(control_f.values.max())
		normalize_control_df = control_f.applymap(lambda x: (np.log(x+1)/max_value))
		normalize_control_df["case_control"] = 1
		normalize_control_df.to_csv("normalized_control.csv",index=False)


if __name__ == "__main__":
	main()
