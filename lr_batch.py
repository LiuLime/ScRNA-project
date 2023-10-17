"""
2023/10/17 revised version
Final batch processing version of logistic regression prediction revised by chatGPT4
The feature selection were conducted by penalty
Liu
"""
import json
import os, sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report, log_loss, confusion_matrix
from sklearn.feature_selection import SelectFromModel
from collections import Counter


def read_merge_file(filepath):
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


def label_data(merged_df):
    x = merged_df.iloc[:, 1:-1].values.astype(float)
    y = merged_df.iloc[:, -1].values
    le = preprocessing.LabelEncoder()
    y = le.fit_transform(y)

    return x, y


def train_and_evaluate(x, y, features):
    # split dataset
    x_train, x_test, y_train, y_test = train_test_split(x, y, stratify=y, random_state=22, test_size=0.3)

    # Scaling
    scaler = MinMaxScaler()
    x_train_scaled = scaler.fit_transform(x_train)
    x_test_scaled = scaler.transform(x_test)

    # parameter search
    param_grid = {'C': [0.001, 0.01, 0.1, 1, 5, 10, 15], 'penalty': ['elasticnet'], 'l1_ratio': [0.5],
                  'solver': ['saga'], 'max_iter': [5000]}
    grid_search = GridSearchCV(LogisticRegression(), param_grid, cv=5)
    grid_search.fit(x_train_scaled, y_train)
    print('grid search sucess')

    # train with searched parameter
    best_estimator = grid_search.best_estimator_
    print(best_estimator)
    sfm = SelectFromModel(best_estimator)
    x_train_sfm = sfm.fit_transform(x_train_scaled, y_train)
    x_test_sfm = sfm.transform(x_test_scaled)
    feature_mask = sfm.get_support()
    print('Before feature selection shape:', x_train_scaled.shape,
          'After feature selection shape:', x_train_sfm.shape)
    lr = LogisticRegression(max_iter=1000)
    lr.fit(x_train_sfm, y_train)

    # evaluation on validation set
    y_pred = lr.predict(x_test_sfm)
    y_pred_prob = lr.predict_proba(x_test_sfm)[:, 1]
    test_loss_ = log_loss(y_test, y_pred)
    conf_mat = confusion_matrix(y_test, y_pred)
    report = classification_report(y_test, y_pred, output_dict=True, zero_division=0)
    coef_repo = report_coef(lr, features, feature_mask)
    return test_loss_, conf_mat, report, coef_repo, best_estimator


def report_coef(model, features, feature_mask):
    coef = model.coef_

    # 对数组进行排序
    selected_features = np.array(features)[feature_mask]
    sorted_array = np.sort(coef, axis=None)[::1]
    sorted_indices = np.argsort(coef, axis=None)[::1].tolist()  # 得到降序排列的索引则用[::-1]
    feature_name = [selected_features[i] for i in sorted_indices]
    sorted_abs_array = abs(sorted_array)
    coef_repo = {'idx': sorted_indices,  # 根据coef排序idx
                 'gene': feature_name,
                 'coef': sorted_array.tolist(),
                 'abs_coef': sorted_abs_array.tolist()}
    return coef_repo


def main(argv):

    FOLDER_PATH = argv[0]
    SAVE_PATH = argv[1]
    FILE_END = argv[2]
    folder_name = FOLDER_PATH.split('/')[-1]

    cell_types = [filename for filename in os.listdir(FOLDER_PATH) if not filename.startswith('.') and filename.endswith(FILE_END)]

    reports = []
    coef_reports = {}
    missing_class_file = []

    for cell_type in cell_types:
        print(f"{cell_type} start process")
        filepath = os.path.join(FOLDER_PATH, cell_type)
        merged_df, features = read_merge_file(filepath)

        # check point
        if merged_df is None and features is None:
            missing_class_file.append(cell_type)
            print(f"Skipped processing for missing label class: {cell_type}")
            continue  # Continue with the next iteration of the loop

        x, y = label_data(merged_df)
        test_loss, conf_mat, report, coef_report, estimator = train_and_evaluate(x, y, features)

        coef_reports.update({cell_type: coef_report})

        reports.append({
            'cell_type': cell_type.rstrip('_merge.csv'),
            'num_samples': len(x),
            'num_genes': x.shape[1],
            'select_genes': len(coef_report['gene']),
            'accuracy': report['accuracy'],
            'test_loss': test_loss,
            'confusion_matrix': f'TN:{conf_mat[0][0]} FP:{conf_mat[0][1]} FN:{conf_mat[1][0]} TP:{conf_mat[1][1]}',
            'estimator': estimator,
            'classification_report': report
        })
        print(f"{cell_type} complete, accuracy={report['accuracy']}")

    # save info
    df = pd.DataFrame(reports)
    df.to_csv(f'{SAVE_PATH}/{folder_name}_report.csv', index=False)
    with open(f'{SAVE_PATH}/{folder_name}_coef.json', 'w') as j:
        json.dump(coef_reports, j)
    if len(missing_class_file) > 0:
        with open(f'{SAVE_PATH}/{folder_name}_missing_class_file.json', 'w') as j:
            json.dump(missing_class_file, j)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python ./main.py {DATA_FOLDER_PATH} {SAVE_PATH} {FILE_END}')
        print(
            'Example: DATA_FOLDER_PATH = /Users/human_TabulaSapiens/exp SAVE_PATH = /Users/human_TabulaSapiens/output/ FILE_END=_merge.csv' )
        sys.exit(1)
    main(sys.argv[1:])
