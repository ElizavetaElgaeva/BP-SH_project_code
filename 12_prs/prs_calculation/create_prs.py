import pandas as pd
import sys

path_to_file = sys.argv[1]
file_name = sys.argv[2]

data = pd.DataFrame()
for i in range(1, 23):
    table = pd.read_csv(f"{path_to_file}/ukb_imp_chr_{i}_v3_PRS.sscore", delim_whitespace=True)
    data["IID"] = table["#IID"]
    data[f"chr{i}"] = table["SCORE1_SUM"]
data.set_index("IID", inplace=True)
data["PRS_values"] = data.sum(axis=1)
data[["PRS_values"]].to_csv(f"{path_to_file}/{file_name}.csv", sep="\t")
