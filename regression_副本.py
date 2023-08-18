import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error,accuracy_score, confusion_matrix, roc_curve, roc_auc_score
import statsmodels.api as sm
#%%
folder_path='/Users/liuyuting/ScRNA_project/data_source/human_TabulaSapiens/exp/Bladder.Fibroblast'
count=1
label={}
for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        gene_name = filename.replace(".txt", "")
        filepath = os.path.join(folder_path, filename)
        if count==1:
            init_df = pd.read_csv(filepath, sep="\t", header=None, names=['cell_id', gene_name, 'label'])
            merge_df = init_df.iloc[:,:-1]
            for idx,row in init_df.iterrows():
                label[row['cell_id']]=row['label']
        else:
            temp_df = pd.read_csv(filepath, sep="\t", header=None, names=['cell_id', gene_name, 'label'])
            for idx,row in temp_df.iterrows():
                label[row['cell_id']]=row['label']
            merge_df=pd.merge(merge_df, temp_df.iloc[:,:-1], how='outer', on='cell_id').fillna('0')

        count+=1

label_df=pd.DataFrame(list(label.items()), columns=['cell_id', 'label'])
merge_df=pd.merge(merge_df,label_df,how='right',on='cell_id')
#%% 特征矩阵和标签
x=merge_df.iloc[:,1:-1].values.astype(float)
y=merge_df.iloc[:,-1].values
le = preprocessing.LabelEncoder()
le.fit(["Young", "Old"])
y_onehot=le.transform(y)
features=list(merge_df.columns)[1:-1]
#%% sklearn 线性回归

x_train, x_test, y_train, y_test = train_test_split(x, y_onehot, random_state=22, test_size=0.2)
# 最小二乘法线性回归
model_sk_linear = LinearRegression()
model_sk_linear.fit(x_train, y_train)
score = model_sk_linear.score(x_test, y_test)
print("R^2 score:", score)
# 测试集预测精确度
y_pred = model_sk_linear.predict(x_test)
mse = mean_squared_error(y_test, y_pred)
print('mse:',mse)
print("Features weights (coefficients):", model_sk_linear.coef_)
print("Intercept (bias):", model_sk_linear.intercept_)


#%% sklearn 逻辑回归模型
x_train, x_test, y_train, y_test = train_test_split(x, y_onehot, random_state=22, test_size=0.2)

clf = LogisticRegression(max_iter=1000)
clf.fit(x_train, y_train)

y_pred = clf.predict(x_test)
y_pred_prob = clf.predict_proba(x_test)[:, 1]
# 准确性
acc = accuracy_score(y_test, y_pred)
print(f"Accuracy: {acc}")
# 混淆矩阵
conf_mat = confusion_matrix(y_test, y_pred)
print(f"Confusion Matrix:\n {conf_mat}")
# ROC曲线和AUC
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
coef=clf.coef_
# 对数组进行排序
sorted_array = np.sort(coef, axis=None)[::1]
sorted_indices = np.argsort(coef, axis=None)[::1]  # 得到降序排列的索引则用[::-1]
feature_name=[features[i] for i in sorted_indices]
sorted_abs_array = abs(sorted_array)
coef_repo={'coef':sorted_array,
        'idx': sorted_indices,
        'gene': feature_name,
        'abs_coef': sorted_abs_array}
coef_repo=pd.DataFrame(coef_repo)
coef_repo.to_csv('Bladder.Fibroblast.csv',index=False,)
