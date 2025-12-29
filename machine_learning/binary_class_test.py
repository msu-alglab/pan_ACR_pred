# Partly created Chat-GPT

import pandas as pd
import numpy as np
import time
from naive_classifier import *

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import (
    GradientBoostingClassifier,
    AdaBoostClassifier,
    RandomForestClassifier,
    HistGradientBoostingClassifier,
    ExtraTreesClassifier
)

from sklearn.dummy import DummyClassifier
from sklearn.linear_model import RidgeClassifier, LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import (
    accuracy_score,
    roc_auc_score,
    precision_score,
    recall_score,
    f1_score,
    cohen_kappa_score,
    average_precision_score
)


# ============================================================
# Inputs
# ============================================================

train_file_path = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50/train.csv"
validate_file_path = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50/validate.csv"
test_file_path = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50/test.csv"

out_path_summary = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50/ml_results/summary.csv"
out_path_true_pred = "/home/mwarr/Data/paper_visualizing/upstream_true_pred_histgbc.txt"

rand_state = 29

save_true_predictions = False
save_summary = True

# ============================================================
# Load Data (remove FIRST COLUMN)
# ============================================================
train_df = pd.read_csv(train_file_path).iloc[:, 1:]

test_df  = pd.read_csv(test_file_path)
test_names = list(test_df.iloc[:, 0])
test_df = test_df.iloc[:, 1:]

# Split features and labels
# Assumes class assignment is the last column
X_train_full = train_df.iloc[:, :-1]
y_train_full = train_df.iloc[:, -1]

X_test_full  = test_df.iloc[:, :-1]
y_test_full  = test_df.iloc[:, -1]


# ============================================================
# Models
# ============================================================
models = {
    "GBC": GradientBoostingClassifier(),
    "HistGBC": HistGradientBoostingClassifier(),
    "AdaBoost": AdaBoostClassifier(),
    "RidgeClassifier": RidgeClassifier(),
    "LDA": LinearDiscriminantAnalysis(),
    "RandomForest": RandomForestClassifier(),
    "ExtraTrees": ExtraTreesClassifier(),
    "Dummy": DummyClassifier(),
    "Naive": NaiveClassifier()
}

cols = ["Model", "Accuracy", "AUC", "Precision", "Recall", "F1", "Kappa", "AUC_PR", "Train_Time", "Test_Time"]
data = {item : [] for item in cols}
true_predictions = []

# ============================================================
# MAIN LOOP
# ============================================================

for model_name, model in models.items():
    print(f"Starting model {model_name}", flush=True)

    # --- TRAINING ---
    start_train = time.time()
    model.fit(X_train_full, y_train_full)
    train_time = time.time() - start_train
    print(f"Training: {train_time} s", flush=True)


    # --- TESTING ---
    start_test = time.time()
    y_pred = model.predict(X_test_full)
    try:
        y_prob = model.predict_proba(X_test_full)[:, 1]
    except:
        y_prob = None
    test_time = time.time() - start_test

    print(f"Testing: {test_time} s", flush=True)

    # --- METRICS ---
    acc     = accuracy_score(y_test_full, y_pred)
    prec    = precision_score(y_test_full, y_pred, zero_division=0)
    rec     = recall_score(y_test_full, y_pred)
    f1      = f1_score(y_test_full, y_pred)
    kappa   = cohen_kappa_score(y_test_full, y_pred)
    auc     = roc_auc_score(y_test_full, y_prob) if y_prob is not None else np.nan
    auc_pr  = average_precision_score(y_test_full, y_prob) if y_prob is not None else np.nan    

    data["Model"].append(model_name)
    data["Accuracy"].append(acc)
    data["Precision"].append(prec)
    data["Recall"].append(rec)
    data["F1"].append(f1)
    data["Kappa"].append(kappa)
    data["AUC"].append(auc)
    data["AUC_PR"].append(auc_pr)
    data["Train_Time"].append(train_time)
    data["Test_Time"].append(test_time)

    if save_true_predictions:
        temp = [test_names[i] for i in range(len(y_pred)) if y_pred[i] == 1]
        true_predictions.extend(temp)

# ============================================================
# Save outputs
# ============================================================
if save_summary:
    pd.DataFrame(data).to_csv(out_path_summary, index=False)
if save_true_predictions:
    with open(out_path_true_pred, "w") as out:
        for reg in true_predictions:
            out.write(f"{reg}\n")
