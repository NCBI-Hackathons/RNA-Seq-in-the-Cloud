import pandas as pd
import sys

case_file = sys.argv[1]
control_file = sys.argv[2]
output_file = sys.argv[3]

case_df = pd.read_csv(case_file)
control_df = pd.read_csv(control_file)

merged_df = pd.concat([case_df, control_df])
merged_df.to_csv(output_file, index=False)

