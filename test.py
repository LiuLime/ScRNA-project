import pandas as pd

data = {"abbr_id": ["ta", "tb", "tb", "tc"],
        "gene": ["g1", "g2", "g1", "g3"],
        "pvalue": [True, True, False, True]}
df = pd.DataFrame(data)
new_df = df.pivot_table(values="pvalue", index="gene", columns="abbr_id", fill_value=0)

# %%
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

# %%
from sql import *

s = sql(db="scrna_mt")
with s:
    results = s.filter_data_by_rank()

# %%
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

# %%
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette
import matplotlib.cm as cm
import matplotlib as mpl

ytdist = np.array([662., 877., 255., 412., 996., 295., 468., 268.,
                   400., 754., 564., 138., 219., 869., 669.])
# colors_cool = cm.get_cmap('jet', 10)
cmap = cm.rainbow(np.linspace(0, 1, 10))
# set_link_color_palette(colors_cool)
set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])
Z = linkage(ytdist, 'single')
dn = dendrogram(Z, color_threshold=230)

dn['color_list']
# %%
import matplotlib as mpl
from matplotlib.pyplot import cm
from scipy.cluster import hierarchy

cmap = cm.rainbow(np.linspace(0, 1, 10))
# hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

print([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

# %%
import matplotlib.pyplot as plt
import numpy as np


def sample_colormap(cmap_name, n_colors):
    """从指定的颜色映射中取样 n 个颜色"""
    cmap = plt.cm.get_cmap(cmap_name)  # 获取颜色映射
    colors = [cmap(i) for i in np.linspace(0, 1, n_colors)]  # 均匀地从颜色映射中取样
    return colors


# 示例：从 'viridis' 颜色映射中取样 10 个颜色
sampled_colors = sample_colormap('viridis', 10)
print(sampled_colors)

# %%
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
from collections import Counter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# a = np.array([[1.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]])
# b = np.array([1.0, 0.0, 0.0, 1.0])
# X = pdist(a, metric="jaccard")

df = pd.DataFrame({"o1": [1.0, 0.0, 1.0, 0.0], "o2": [1.0, 0.0, 0.0, 1.0], "o3": [1.0, 0.0, 0.0, 1.0]},
                  index=["idx1", "idx2", "idx3", "idx4"])
X2 = linkage(df, method='ward')
dn = dendrogram(X2, no_plot=True)
f = fcluster(X2, 3, criterion='maxclust')
f_count = Counter(f)
fig, ax = plt.subplots()
hmp=sns.heatmap(df, ax=ax)
ax.axhline(1,color='orange',linestyle='--')
ax.axvline(2,color='green',linestyle='--')


# %%
import pandas as pd

data = {"abbr_id": ["ta", "tb", "tb", "tc"],
        "gene": ["g1", "g2", "g1", "g3"],
        "value": [1, 1, 3, 4]
        }
df = pd.DataFrame(data)
sort_index = [0, 3, 2, 1]
df_row = df.iloc[sort_index]
df_row2 = df[[True, False, True, False]]
df_col = df[["gene", "abbr_id", "value"]]

#%%
import networkx as nx
import pandas as pd
edges=pd.DataFrame({"source":["a","b","c"],"target":["b","a","b"],"notice":[True,True,False]})
G = nx.from_pandas_edgelist(edges, edge_attr=True)
