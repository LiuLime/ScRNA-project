"""
2023/9/29
Initial version of regression coding for reading one files

"""


import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler, MinMaxScaler  # 特征缩放
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV  # 交叉验证
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import mean_squared_error, accuracy_score, f1_score, confusion_matrix, classification_report, \
    roc_curve, roc_auc_score, log_loss
import seaborn as sns

# %%  数据准备
# folder_path = '/Users/liuyuting/WorkSpace/ScRNA_project/data_source/human_TabulaSapiens/exp/Bladder.Smooth_muscle_cells'
folder_path = '/Users/liuyuting/WorkSpace/ScRNA_project/data_source/human_TabulaSapiens/exp/'


folder_path = '/Users/liuyuting/WorkSpace/ScRNA_project/data_source/human_TabulaSapiens/exp/Bladder.Smooth_muscle_cells'
count = 1
label = {}
for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        gene_name = filename.replace(".txt", "")
        filepath = os.path.join(folder_path, filename)
        if count == 1:
            init_df = pd.read_csv(filepath, sep="\t", header=None, names=['cell_id', gene_name, 'label'])
            merge_df = init_df.iloc[:, :-1]
            for idx, row in init_df.iterrows():
                label[row['cell_id']] = row['label']
        else:
            temp_df = pd.read_csv(filepath, sep="\t", header=None, names=['cell_id', gene_name, 'label'])
            for idx, row in temp_df.iterrows():
                label[row['cell_id']] = row['label']
            merge_df = pd.merge(merge_df, temp_df.iloc[:, :-1], how='outer', on='cell_id').fillna('0')
        count += 1
label_df = pd.DataFrame(list(label.items()), columns=['cell_id', 'label'])
merge_df = pd.merge(merge_df, label_df, how='right', on='cell_id')
merge_df.to_csv('merge_df.csv', index=False)


merge_df = pd.read_csv('merge_df.csv')
x = merge_df.iloc[:, 1:-1].values.astype(float)
y = merge_df.iloc[:, -1].values
# 标签one-hot转换
le = preprocessing.LabelEncoder()
le.fit(["Old", "Young"])
y_onehot = le.transform(y)
features = list(merge_df.columns)[1:-1]

# merge_df.to_csv('merge_df.csv', index=False)
# %% sklearn 线性回归
#
# x_train, x_test, y_train, y_test = train_test_split(x, y_onehot, random_state=22, test_size=0.2)
# # 最小二乘法线性回归
# model_sk_linear = LinearRegression()
# model_sk_linear.fit(x_train, y_train)
# score = model_sk_linear.score(x_test, y_test)
# print("R^2 score:", score)
# # 测试集预测精确度
# y_pred = model_sk_linear.predict(x_test)
# mse = mean_squared_error(y_test, y_pred)
# print('mse:',mse)
# print("Features weights (coefficients):", model_sk_linear.coef_)
# print("Intercept (bias):", model_sk_linear.intercept_)


# %% sklearn 逻辑回归模型
# 划分数据集
x_train, x_test, y_train, y_test = train_test_split(x, y, stratify=y, random_state=22,
                                                    test_size=0.2)  # stratify表示正反例子抽样的分布与y_onehot相同
# x数据缩放
scaler = MinMaxScaler()
x_train_scaled = scaler.fit_transform(x_train)
x_test_scaled = scaler.transform(x_test)  # 测试集缩放的参数分布是转移的训练集的，避免测试集数据泄漏
# %%
"""Feature selection based on L1 Penalty"""
# 基模型特征选择

solver = 'saga'
penalty = "elasticnet"
if penalty == 'elasticnet':
    clf = LogisticRegression(max_iter=1000, solver=solver, penalty=penalty, l1_ratio=0.5)
else:
    clf = LogisticRegression(max_iter=1000, solver=solver, penalty=penalty)
param_grid = {'C': [0.001, 0.01, 0.1, 1, 5, 10, 15]}

grid_search = GridSearchCV(clf, param_grid, cv=5)  # 分类问题默认score=accuracy
grid_search.fit(x_train_scaled, y_train)
best_params = grid_search.best_params_
best_C = best_params['C']
print("Best C:", best_C)
if penalty == 'elasticnet':
    clf_bestC = LogisticRegression(max_iter=1000, solver=solver, penalty=penalty, l1_ratio=0.5, C=best_C)
    print('l1_ratio set')
else:
    clf_bestC = LogisticRegression(max_iter=1000, solver=solver, penalty=penalty, C=best_C)



# %%

# sfm = SelectFromModel(clf)
sfm = SelectFromModel(clf_bestC)
x_train_sfm = sfm.fit_transform(x_train_scaled, y_train)
x_test_sfm = sfm.transform(x_test_scaled)
feature_mask = sfm.get_support()
print('Before feature selection shape:', x_train_scaled.shape,
      'After feature selection shape:', x_train_sfm.shape)

# 5折交叉验证->用来调整超参数的（optional）
scores = cross_val_score(clf_bestC, x_train_scaled, y_train, cv=5)
print("Cross-validation scores: ", scores)
print("Average cross-validation score: ", scores.mean())

lr = LogisticRegression(max_iter=1000)
lr.fit(x_train_sfm, y_train)

# 测试集进行评估
y_pred = lr.predict(x_test_sfm)
y_pred_prob = lr.predict_proba(x_test_sfm)[:, 1]

## loss
test_loss = log_loss(y_test, y_pred)
print(f"Test loss: {test_loss}")
## 准确性
acc = accuracy_score(y_test, y_pred)
print(f"Accuracy: {acc}")
## f1score
f1 = f1_score(y_test, y_pred)
print(f'f1_score:{f1}')
## 混淆矩阵
conf_mat = confusion_matrix(y_test, y_pred)
print(f"Confusion Matrix:\n {conf_mat}")
sns.heatmap(conf_mat, annot=True, fmt="d")
## ROC曲线和AUC
fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
plt.figure(figsize=(10, 6))
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc_score(y_test, y_pred_prob))
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC)')
plt.legend(loc="lower right")
plt.show()
## report
report = classification_report(y_test, y_pred, output_dict=True)
for label, metrics in report.items():
    if isinstance(metrics, dict):
        print(f"{label}:", {metric: round(value, 3) for metric, value in metrics.items()})
    else:
        print(f"{label}: {round(metrics, 3)}")

# %% 权重排序
coef = lr.coef_
# 对数组进行排序
selected_features = np.array(features)[feature_mask]
sorted_array = np.sort(coef, axis=None)[::1]
sorted_indices = np.argsort(coef, axis=None)[::1]  # 得到降序排列的索引则用[::-1]
feature_name = [selected_features[i] for i in sorted_indices]
sorted_abs_array = abs(sorted_array)
coef_repo = {'idx': sorted_indices,
             'gene': feature_name,
             'coef': sorted_array,
             'abs_coef': sorted_abs_array}
coef_repo = pd.DataFrame(coef_repo)
coef_repo.to_csv('Bladder.Smooth_muscle_cells_FS_L.csv', index=False, )