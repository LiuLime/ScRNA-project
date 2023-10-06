"""
2023/10/2
Final batch processing version of logistic regression prediction revised by chatGPT4
The feature selection were conducted by penalty

"""
import json
import os, sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.metrics import classification_report, log_loss, confusion_matrix
from sklearn.feature_selection import SelectFromModel

# 设置文件夹路径
# FOLDER_PATH = '/Users/liuyuting/WorkSpace/ScRNA_project/data_source/human_TabulaSapiens/exp/'
# SAVE_PATH = '/Users/liuyuting/WorkSpace/ScRNA_project/data_source/human_TabulaSapiens/output/'

"""这里的优化思路是使用temp_df读取完所有基因matrix并放入list中，在一起合并并删除重复cell_id，来取代迭代row的耗时做法"""

# 合并文件夹内的数据
def merge_data_in_folder(filepath):
    # 读取文件夹内所有txt文件
    filenames = [f for f in os.listdir(filepath) if f.endswith(".txt")]

    dfs = []
    label_dict = {}

    for filename in filenames:
        file_path = os.path.join(filepath, filename)
        gene_name = filename.replace(".txt", "")

        temp_df = pd.read_csv(file_path, sep="\t", header=None, names=['cell_id', gene_name, 'label'])  # names自定义列
        label_dict.update(temp_df.set_index('cell_id')['label'].to_dict())  # {cell_id:label}

        dfs.append(temp_df.set_index('cell_id').drop(columns='label'))

    merged_df = pd.concat(dfs, axis=1)  # 横向拼接前set_index，相同cell_id自动合并

    label_df = pd.DataFrame(list(label_dict.items()), columns=['cell_id', 'label']).set_index('cell_id')
    merged_df = pd.concat([merged_df, label_df], axis=1).reset_index().rename(columns={'index': 'cell_id'})
    features = list(merged_df.columns)[1:-1]
    print('merge success')
    return merged_df, features


# 对数据进行标记
def label_data(merged_df):
    x = merged_df.iloc[:, 1:-1].values.astype(float)
    y = merged_df.iloc[:, -1].values
    le = preprocessing.LabelEncoder()
    y = le.fit_transform(y)
    print('label success')
    return x, y


# 特征选择和训练模型
def train_and_evaluate(x, y, features):
    # split dataset
    x_train, x_test, y_train, y_test = train_test_split(x, y, stratify=y, random_state=22, test_size=0.3)
    print('split sucess')
    
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


# 主逻辑
def main(argv):
    FOLDER_PATH = argv[0]
    SAVE_PATH = argv[1]
    cell_types = [foldername for foldername in os.listdir(FOLDER_PATH) if not foldername.startswith('.')]
    reports = []
    coef_reports = {}
    for cell_type in cell_types:
        print(f"{cell_type} start process")
        filepath = os.path.join(FOLDER_PATH, cell_type)
        merged_df, features = merge_data_in_folder(filepath)

        x, y = label_data(merged_df)
        test_loss, conf_mat, report, coef_report, estimator = train_and_evaluate(x, y, features)

        # 将结果存入报告列表
        coef_reports.update({cell_type: coef_report})

        reports.append({
            'cell_type': cell_type,
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

    df = pd.DataFrame(reports)
    df.to_csv(f'{SAVE_PATH}/sum_report.csv', index=False)
    with open(f'{SAVE_PATH}/coef_reports.json', 'w') as j:
        json.dump(coef_reports, j)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python ./main.py {DATA_FOLDER_PATH} {SAVE_PATH}')
        print('Example: DATA_FOLDER_PATH = /Users/human_TabulaSapiens/exp/ SAVE_PATH = /Users/human_TabulaSapiens/output/')
        sys.exit(1)
    main(sys.argv[1:])

