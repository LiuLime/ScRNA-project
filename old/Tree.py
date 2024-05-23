import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import mean_squared_error, accuracy_score, confusion_matrix, roc_curve, roc_auc_score



merge_df = pd.read_csv('merge_df.csv')
# %% 特征矩阵和标签
x = merge_df.iloc[:, 1:-1].values.astype(float)
y = merge_df.iloc[:, -1].values

# 标签one-hot转换
le = preprocessing.LabelEncoder()
le.fit(["Old","Young"])
y_onehot = le.transform(y)
features = list(merge_df.columns)[1:-1]

# %% 树模型
# 划分数据集
x_train, x_test, y_train, y_test = train_test_split(x, y_onehot, stratify=y_onehot, random_state=22,
                                                    test_size=0.2)  # stratify表示正反例子抽样的分布与y_onehot相同
clf = RandomForestClassifier(n_estimators=100, random_state=42)
cv_scores = cross_val_score(clf, x_train, y_train, cv=5)
print("交叉验证准确率：", cv_scores)
print("平均交叉验证准确率：", np.mean(cv_scores))

clf.fit(x_train, y_train)
y_pred = clf.predict(x_test)
y_pred_prob = clf.predict_proba(x_test)[:, 1]
acc = accuracy_score(y_test, y_pred)
print("测试集准确率：", acc)
## 混淆矩阵
conf_mat = confusion_matrix(y_test, y_pred)
print(f"Confusion Matrix:\n {conf_mat}")
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


#%%
# 将特征重要性排序并绘图
# 特征重要性
importances = clf.feature_importances_
sorted_indices = np.argsort(importances)[::-1]
sorted_array = np.sort(importances, axis=None)[::-1]
feature_name = [features[i] for i in sorted_indices]

coef_repo = {'importance': sorted_array,
             'idx': sorted_indices,
             'gene': feature_name,
             }
coef_repo = pd.DataFrame(coef_repo)
coef_repo.to_csv('Bladder.Fibroblast_tree.csv', index=False, )
