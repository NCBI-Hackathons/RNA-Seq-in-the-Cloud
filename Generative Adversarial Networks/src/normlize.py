import argparse
import pandas as pd
import numpy as np
import json


def normalize(args):
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
    with open(args.stats_file, "w+") as f:
        f.write(json.dumps({"mean": mean, "std": std}))

    normalize_df = X_df.applymap(lambda x: (x-mean)/std)

    Xtarget=np.array([0]*cases.shape[0] + [1]*controls.shape[0])
    
    normalize_df["case_control"] = Xtarget
    print("normalize_df:", normalize_df.shape)
    normalize_df.to_csv(args.output_file,index=False)


def discretize(x):
    if x < 0.5:
        return 0
    else:
        return 1
    

def denormalize(args):
    df=np.loadtxt(args.input_file, delimiter=",", header=None)
    print("input:", df.shape)

    case_or_control_column = df.columns[-1]
    case_or_control_series0 = df.loc[case_or_control_column]
    Xtarget = case_or_control_series0.apply(discretize)
    X_df = df.drop(case_or_control_column)
    
    stats = json.load(args.stats_file)
    mean = stats["mean"]
    std = stats["std"]
    
    print("X_df:", X_df.shape)
    
    X_df = X_df.applymap(lambda x: x * std + mean)
    denormalize_df = X_df.applymap(lambda x: np.exp(x) - 1)
    
    denormalize_df["case_control"] = Xtarget
    print("denormalize_df:", denormalize_df.shape)
    denormalize_df.to_csv(args.output_file,index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers()

    parser_normalize = subparsers.add_parser("normalize")
    parser_normalize.add_argument("--case_file", type=str,required=True)
    parser_normalize.add_argument("--control_file", type=str,required=True)
    parser_normalize.add_argument("--output_file", type=str,required=True)
    parser_normalize.add_argument("--stats_file", type=str, required=True)
    parser_normalize.add_argument("--balance", action="store_true")
    parser_normalize.set_defaults(func=normalize)
    
    parser_denormalize = subparsers.add_parser("denormalize")
    parser_denormalize.add_argument("--input_file", type=str,required=True)
    parser_denormalize.add_argument("--output_file", type=str,required=True)
    parser_denormalize.add_argument("--stats_file", type=str, required=True)
    parser_denormalize.set_defaults(func=denormalize)
    
    args = parser.parse_args()
    
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
