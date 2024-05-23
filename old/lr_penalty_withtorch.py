"""
2023/10/19 version
For logistic regression model training by pytorch with elasticnet penalty feature selection method

@ Liuyuting

"""

import json
import os, sys
import pandas as pd
import numpy as np
from collections import Counter
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from torch.utlis.data import DataLoader

def read_merge_file(filepath):
    """Read expression file and check label counts,
    files whose label types less than 2 and each label counts less than 10 will not be read"""
    merged_df = pd.read_csv(filepath, header=0).fillna(0)
    labels = merged_df.loc[:, 'label']
    label_class = set(labels)

    # checkpoint1 有两种label
    if len(label_class) <= 1:
        return None, None
    # checkpoint2 每种label数量大于10
    for element, count in Counter(labels).items():
        if count <= 10:
            return None, None

    features = list(merged_df.columns)[1:-1]
    return merged_df, features


def label_dataset(merged_df):
    """encode old and young label as 0-1 label"""
    x = merged_df.iloc[:, 1:-1].values.astype(float)
    y = merged_df.iloc[:, -1].values
    le = preprocessing.LabelEncoder()
    y = le.fit_transform(y)

    return x, y


def prepare_dataset(x, y, features):
    """
    Do train/test dataset split
    Scaling

    """

    # split dataset
    x_train, x_test, y_train, y_test = train_test_split(x, y, stratify=y, random_state=22, test_size=0.3)

    # Scaling
    scaler = MinMaxScaler()
    x_train_scaled = scaler.fit_transform(x_train)
    x_test_scaled = scaler.transform(x_test)