import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, r2_score
from sklearn.metrics import classification_report
from scipy import stats
import sys
import csv
import os
import warnings
import pickle
warnings.filterwarnings('ignore')

prs_file = sys.argv[1]
phenotype_file = sys.argv[2]
trait = sys.argv[3]
path_to_val_list = "/home/anostaeva/Back_Pain_project/Back_Pain_2022/phenotype/test_list.txt"

ID = "ID"
IID = "IID"
PRS = "PRS_values"
PHENOTYPE = trait

PRS_COLUMNS = [ID, PRS]
prs_data = pd.read_csv(prs_file, delim_whitespace=True)
prs_data.columns = PRS_COLUMNS
prs_data.set_index(ID, inplace=True)
prs_data = prs_data.dropna(axis=0)

phenotype_data = pd.read_csv(phenotype_file, sep="\t", dtype=str)
phenotype_data[ID] = phenotype_data[IID].astype(int)
phenotype_data.set_index(ID, inplace=True)

val_ids = pd.read_csv(path_to_val_list)
val_ids[ID] = val_ids[IID].astype(int)
val_ids.set_index(ID, inplace=True)

phenotype_data = phenotype_data.loc[phenotype_data.index.isin(prs_data.index)]
phenotype_data = phenotype_data.loc[phenotype_data.index.isin(val_ids.index)]
prs_data = prs_data.loc[prs_data.index.isin(phenotype_data.index)]
phenotype_data.sort_index(inplace=True)
prs_data.sort_index(inplace=True)
if PHENOTYPE == "back":
    phenotype_data[PHENOTYPE] = phenotype_data[PHENOTYPE].astype(int)
else:
    phenotype_data[PHENOTYPE] = phenotype_data[PHENOTYPE].astype(float)
prs_data["Age"] = phenotype_data["Age"]
prs_data["Sex"] = phenotype_data["Sex"]

def perf_measure(y_actual, y_hat):
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    for i in range(len(y_hat)):
        if y_actual[i]==y_hat[i]==1:
           TP += 1
        if y_hat[i]==1 and y_actual[i]!=y_hat[i]:
           FP += 1
        if y_actual[i]==y_hat[i]==0:
           TN += 1
        if y_hat[i]==0 and y_actual[i]!=y_hat[i]:
           FN += 1

    return(TP, FP, TN, FN)


def run_log_reg(name, columns, df):
    data = prs_data[columns]
    X = data
    y = phenotype_data[PHENOTYPE]

    model=LogisticRegression()
    model.fit(X, y)
#    print("Intercept: ", model.intercept_)
#    print("Coefficient: ", model.coef_)
    score=accuracy_score(y, model.predict(X))
    roc=roc_auc_score(y, model.predict_proba(X)[:,1])

    df.loc[df["model"]==name, "AUC"] = roc

    q1 = roc/(2-roc)
    q2 = 2*roc**2/(1+roc)
    n1 = phenotype_data.loc[phenotype_data[PHENOTYPE]==1].shape[0]
    n2 = phenotype_data.loc[phenotype_data[PHENOTYPE]==0].shape[0]
    se_auc = np.sqrt((roc*(1-roc)+(n1-1)*(q1-roc**2)+(n2-1)*(q2-roc**2))/(n1*n2))
    df.loc[df["model"]==name, 'CI_left'] = roc-1.96*se_auc
    df.loc[df["model"]==name, 'CI_right'] = roc+1.96*se_auc

    df_y = pd.DataFrame({"y" : list(y), "probs" : list(model.predict_proba(X)[:,1])})
    df_y['procentile'] = pd.qcut(df_y.rank(method='first')["probs"], q = 20, labels = False)
    df_y["y_pred"] = 0
    df_y.loc[df_y["procentile"]==19, "y_pred"] = 1

    TP, FP, TN, FN = perf_measure(list(df_y["y"]), list(df_y["y_pred"]))
    df.loc[df["model"]==name, 'TP'] = TP
    df.loc[df["model"]==name, 'FP'] = FP
    df.loc[df["model"]==name, 'TN'] = TN
    df.loc[df["model"]==name, 'FN'] = FN
    df.loc[df["model"]==name, 'Sensitivity'] = TP/(TP+FN)
    df.loc[df["model"]==name, 'Specificity'] = TN/(TN+FP)


def run_linear_reg(name, columns, df):
    data = prs_data[columns]
    X = data
    y = phenotype_data[PHENOTYPE]
    N = len(X)
    p = len(columns) + 1

    model=LinearRegression()
    model.fit(X, y)

    y_predict = model.predict(X)
    r2 = r2_score(y, y_predict)
    rho, pval = stats.spearmanr(y, y_predict)

    df.loc[df["model"]==name, "R2"] = r2
    df.loc[df["model"]==name, "R_spearman"] = rho

    n = len(y)
    k = len(columns)
    se2 = (4*r2*((1-r2)**2)*((n-k-1)**2))/((n**2-1)*(n+3))
    se = se2 ** 0.5
    df.loc[df["model"]==name, "R2_lower_formula"] = r2 - 2 * se
    df.loc[df["model"]==name, "R2_upper_formula"] = r2 + 2 * se


if PHENOTYPE == "back":
    df_re_estimate = pd.DataFrame({'model':["PRS_age_sex", "PRS_age", "PRS", "age"]}, columns=["model", "AUC", "CI_left", "CI_right", "TP", "FP", "TN", "FN", "Sensitivity", "Specificity"])

    for name, columns in zip(["PRS_age_sex", "PRS_age", "PRS", "age"], [[PRS, "Age", "Sex"], [PRS, "Age"], [PRS], ["Age"]]):
        run_log_reg(name, columns, df_re_estimate)

    df_re_estimate.to_csv(f"./data/{trait}_test_results.csv", mode='a', index=False, header=False)
else:
    df_re_estimate = pd.DataFrame({'model':["PRS_age_sex", "PRS_age", "PRS", "age", "sex"]}, columns=["model", "R2", "R_spearman", "R2_lower_formula", "R2_upper_formula"])

    for name, columns in zip(["PRS_age_sex", "PRS_age", "PRS", "age", "sex"], [[PRS, "Age", "Sex"], [PRS, "Age"], [PRS], ["Age"], ["Sex"]]):
        run_linear_reg(name, columns, df_re_estimate)

    df_re_estimate.to_csv(f"./data/{trait}_test_results.csv", mode='a', index=False, header=False)

