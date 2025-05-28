import numpy as np
from sklearn.model_selection import StratifiedKFold,StratifiedGroupKFold,StratifiedShuffleSplit, ShuffleSplit
from sklearn.model_selection import train_test_split
import glob
import pandas as pd
import csv
import numpy as np
from collections import Counter

seed=7

labels=pd.read_csv("csv_files/lgggbm_surv_labels_v3.csv")
labels = labels.drop_duplicates(['case_id']).copy()

X=labels['case_id'].to_list()
y=labels['survival_months']

disc_labels, q_bins = pd.qcut(y, q=4, retbins=True, labels=False)

print(disc_labels)
y=disc_labels.to_list()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=seed, stratify=y)
with open("train_cases.txt", "w") as file:
    file.write("\n".join(X_train))
with open("val_cases.txt", "w") as file:
    file.write("\n".join(X_test))
X_train=np.array(X_train)

rs=StratifiedKFold(n_splits=10,shuffle=True,random_state=seed)


split=0
for i, (train_index, test_index) in enumerate(rs.split(X_train,y_train)):
  x_train=X_train[train_index]
  x_test=X_train[test_index]
  data={"train":list(set(x_train)),"val":list(set(x_test))}
  df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in data.items()]))
  df.to_csv(f"splits/splits_{split}.csv")
  split+=1
