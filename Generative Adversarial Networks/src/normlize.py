import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--case_file", type=str,required=True)
parser.add_argument("--control_file", type=str,required=True)
parser.add_argument("--output_file", type=str,required=True)
parser.add_argument("--balance", action="store_true")
args = parser.parse_args()

cases=np.loadtxt(args.case_file, delimiter=",")
print("cases:", cases.shape)
controls=np.loadtxt(args.control_file, delimiter=",")
print("controls:", controls.shape)
    
if args.balance:
    # balance the data: there are fewer controls
    n = min(cases.shape[0], controls.shape[0])
    cases=cases[0:n]
    controls=controls[0:n]

#normalize
X_df=pd.concat([pd.DataFrame(cases),pd.DataFrame(controls)])
print("X_df:", X_df.shape)

X_df = X_df.applymap(lambda x: np.log(x + 1))
mean = np.mean(X_df.values)
std = np.std(X_df.values)
normalize_df = X_df.applymap(lambda x: (x-mean)/std)

Xtarget=np.array([0]*cases.shape[0] + [1]*controls.shape[0])

normalize_df["case_control"] = Xtarget
print("normalize_df:", normalize_df.shape)
normalize_df.to_csv(args.output_file,index=False)


