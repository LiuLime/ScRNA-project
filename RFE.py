import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import MinMaxScaler  # 特征缩放
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel, RFECV, RFE
from sklearn.metrics import mean_squared_error, accuracy_score, confusion_matrix, roc_curve, roc_auc_score




#%% 特征矩阵和标签
merge_df = pd.read_csv('merge_df.csv')
x = merge_df.iloc[:, 1:-1]
y = merge_df.iloc[:, -1]
# x = merge_df.iloc[:, 1:-1].values.astype(float)
# y = merge_df.iloc[:, -1].values

# 标签one-hot转换
le = preprocessing.LabelEncoder()
le.fit(["Old", "Young"])
y_onehot = le.transform(y)
features = list(merge_df.columns)[1:-1]

# 划分数据集
x_train, x_test, y_train, y_test = train_test_split(x, y_onehot, stratify=y_onehot, random_state=22,
                                                    test_size=0.2)  # stratify表示正反例子抽样的分布与y_onehot相同
# x数据缩放
scaler = MinMaxScaler()
x_train_scaled = scaler.fit_transform(x_train)
x_test_scaled = scaler.transform(x_test)  # 测试集缩放的参数分布是转移的训练集的，避免测试集数据泄漏

#%%
clf=LogisticRegression()
rfecv = RFECV(estimator=clf, cv=5, scoring='accuracy',n_jobs=10)
x_train_scaled_slt = rfecv.fit_transform(x_train_scaled, y_train)
x_test_scaled_slt = rfecv.transform(x_test_scaled)
# select_features = [features[idx] for idx, value in enumerate(rfecv.support_) if value]

clf.fit(x_train_scaled_slt, y_train)
# 打印最佳特征数量
print("最佳特征数量:", rfecv.n_features_)

# 打印最佳特征的排名
print("特征排名:", rfecv.ranking_)

#%%
# Plot number of features VS. cross-validation scores
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score of number of selected features")
plt.plot(range(1, len(rfecv.cv_results_['mean_test_score']) + 1), rfecv.cv_results_['mean_test_score'])
plt.show()

#%% selected features
# 测试集进行评估
y_pred = clf.predict(x_test_scaled_slt)
y_pred_prob = clf.predict_proba(x_test_scaled_slt)[:, 1]
## 准确性
acc = accuracy_score(y_test, y_pred)
print(f"Accuracy: {acc}")
## 混淆矩阵
conf_mat = confusion_matrix(y_test, y_pred)
print(f"Confusion Matrix:\n {conf_mat}")
sns.heatmap(conf_mat, annot=True, fmt="d")
# %% 权重排序
coef = clf.coef_
# 对数组进行排序
selected_features = np.array(features)[rfecv.get_support()]
sorted_array = np.sort(coef, axis=None)[::1]
sorted_indices = np.argsort(coef, axis=None)[::1]  # 得到降序排列的索引则用[::-1]
feature_name = [selected_features[i] for i in sorted_indices]
sorted_abs_array = abs(sorted_array)
coef_repo = {'coef': sorted_array,
             'idx': sorted_indices,
             'gene': feature_name,
             'abs_coef': sorted_abs_array}
coef_repo = pd.DataFrame(coef_repo)
coef_repo.to_csv('Bladder.Fibroblast_FSrfecv.csv', index=False, )



