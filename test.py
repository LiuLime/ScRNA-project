import pandas as pd

data = {"abbr_id": ["ta", "tb", "tb", "tc"],
        "gene": ["g1", "g2", "g1", "g3"],
        "pvalue": [True, True, False, True]}
df = pd.DataFrame(data)
new_df = df.pivot_table(values="pvalue", index="gene", columns="abbr_id", fill_value=0)

#%%
import numpy as np
from sklearn.cluster import AgglomerativeClustering

# 示例数据同上
# data = np.array([[1, 2], [1, 4], [1, 0],
#                  [10, 2], [10, 4], [10, 0]])

# 创建AgglomerativeClustering实例
clustering = AgglomerativeClustering(n_clusters=2, metric='euclidean', linkage='ward')

# 拟合模型
clustering.fit(new_df)

# 获取聚类标签
labels = clustering.labels_

print("Cluster labels:", labels)

#%%
from sql import *

s =sql(db = "scrna_mt")
with s:
        results = s.filter_data_by_rank()

#%%
# Import the dendrogram function and the ward clustering function from SciPy
from scipy.cluster.hierarchy import dendrogram, ward
from sklearn.datasets import make_blobs

X, y = make_blobs(random_state=0, n_samples=12)
# Apply the ward clustering to the data array X
# The SciPy ward function returns an array that specifies the distances
# bridged when performing agglomerative clustering
linkage_array = ward(X)
# Now we plot the dendrogram for the linkage_array containing the distances
# between clusters
dendrogram(linkage_array)

# mark the cuts in the tree that signify two or three clusters
ax = plt.gca()
bounds = ax.get_xbound()
ax.plot(bounds, [7.25, 7.25], '--', c='k')
ax.plot(bounds, [4, 4], '--', c='k')

ax.text(bounds[1], 7.25, ' two clusters', va='center', fontdict={'size': 15})
ax.text(bounds[1], 4, ' three clusters', va='center', fontdict={'size': 15})
plt.xlabel("Sample index")
plt.ylabel("Cluster distance")


