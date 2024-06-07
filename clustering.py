"""
Define marker pattern by clustering


@ 2024/5/31 Liu
"""

from sql import *
import pandas as pd
import matplotlib.pyplot as plt
from utils import common
from collections import Counter
import json


def replace_nonzero(x):
    return 1 if x != 0 else 0


def is_marker(x, markers: list):
    return True if x in markers else False


with open("./config.json") as j:
    config = json.load(j)
    markers = config["markers"]
    loadPath = config["loadPath"]

db = "scrna_mt"
save_path = loadPath[db]["save_path_clus"]

cut_number = 21
save_figure_hier = "hierarchy_ward"
save_cluster_hier = f"hierarchy_ward_clus{cut_number}"
# %%
# 取出log10p top 2000
s = sql(db=db)
with s:
    top_p = pd.DataFrame(s.filter_data_by_rank(), columns=["id", "gene", "log10p", "abbr_id"])

# 转换为gene*tissue 0-1 形式
top_p.loc[:, "is_in"] = top_p["log10p"].apply(replace_nonzero)
pivot_df = top_p.pivot_table(values="is_in", index="gene", columns="abbr_id", fill_value=0)

X = pivot_df.values
y = pivot_df.index

# %% 层次聚类
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

linked = linkage(X, method='ward')  # 'single','complete', 'average', 'ward'

# scipy绘制树状图
plt.figure(figsize=(100, 70))
dendrogram(linked,
           orientation='top',
            
           distance_sort='descending',
           show_leaf_counts=False)
plt.savefig(f"{save_path}/{save_figure_hier}.pdf")
plt.show()
plt.close()
# %% fcluster标记cluster label

f1 = fcluster(linked, cut_number, criterion='maxclust')
f1_count = Counter(f1)
f1_count = pd.DataFrame(f1_count.items(), columns=["cluster_idx", "gene_nums"])
hier_cluster_ = pd.DataFrame({"gene": y, "cluster": f1})
hier_cluster_.loc[:, "is_marker"] = hier_cluster_["gene"].apply(is_marker, args=(markers,))
common.save_csv(hier_cluster_, f"{save_path}/{save_cluster_hier}.csv")
common.save_csv(f1_count, f"{save_path}/{save_cluster_hier}_counts.csv")
print(f"Clusters number: {f1_count}")

# %%
# kmeans
# from sklearn.cluster import KMeans

# krange = 2001
# save_title = "kmeans_cluster_2000"
#
# inertia = []
# for k in range(1, krange, 10):
#     kmeans = KMeans(n_clusters=k, random_state=0, n_init="auto")
#     kmeans.fit(X)
#     inertia.append(kmeans.inertia_)
#
# cluster_res = pd.DataFrame({"gene": y, "labels": kmeans.labels_})
# common.save_csv(cluster_res, f"{save_path}/{save_title}.csv")
#
# plt.plot(range(1, krange, 10), inertia, marker='o')
# plt.xlabel('Number of clusters')
# plt.ylabel('Inertia')
# plt.title('Elbow Method')
# plt.savefig(f"{save_path}/{save_title}.pdf")
# plt.show()
# plt.close()

# %%
# # 使用scikit-learn确定聚类
# from sklearn.cluster import AgglomerativeClustering
# cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
# cluster.fit_predict(X)
#
# # 可视化聚类结果
# plt.scatter(X[:, 0], X[:, 1], c=cluster.labels_, cmap='rainbow')
# plt.show()
# %% PCA
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt
# import mglearn
# from sklearn.preprocessing import StandardScaler, RobustScaler
#
#
# def pca_fig(X1):
#     plt.figure(figsize=(8, 8))
#     mglearn.discrete_scatter(X1[:, 0], X1[:, 1], X1)
#     plt.legend(y, loc="best")
#     plt.gca().set_aspect("equal")
#     plt.xlabel("First principal component")
#     plt.ylabel("Second principal component")
#
#
# pca = PCA(n_components=2, random_state=22)
#
# # 没有scaler
# X_pca = pca.fit_transform(X)
# plt.scatter(X_pca[:, 0], X_pca[:, 1])
# # scaler
# # scaler = RobustScaler()
# # X_scaler = scaler.fit_transform(X)
# # X_scaler = pca.fit_transform(X_scaler)
# # plt.scatter(X_scaler[:, 0], X_scaler[:, 1])
