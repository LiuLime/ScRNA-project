"""
Define marker pattern by clustering

@ 2024/5/31 Liu
"""

from sql import *
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import json

from utils import common
from DrawTool import draw_dendrogram, _dendro, draw_dendro_with_heatmap


def replace_nonzero(x):
    return 1 if x != 0 else 0


def is_marker(x, markers: list):
    return True if x in markers else False


with open("./config.json") as j:
    config = json.load(j)
    markers = config["markers"]
    loadPath = config["loadPath"]

db = "scrna_mc"
save_path = loadPath[db]["save_path_clus"]

cut_threshold = 25
linkage_method = "ward"  # "ward"得用euclidean metric
pdist_metric = "euclidean"  # distance metric for binary data: "jaccard","hamming","yule","dice"
save_figure_hier = f"hier_{linkage_method}_{pdist_metric}_hmp_dist{cut_threshold}_top1000_old"  # top1000还是top100，old还是young，在pickTopP.sql文件里最后一行改。
save_cluster_hier = f"hier_{linkage_method}_{pdist_metric}_hmp_dist{cut_threshold}_top1000_old"
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
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette, cophenet
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

# Perform hierarchical/agglomerative clustering.
Z = linkage(X, method=linkage_method, metric=pdist_metric)  # 'single','complete', 'average', 'ward'


def check_cluster_effect(Z, X):
    """check the Cophenetic Correlation Coefficient. This compares (correlates) the actual pairwise distances
    of all samples to those implied by the hierarchical clustering."""
    c, coph_dists = cophenet(Z, pdist(X))
    print(c)


check_cluster_effect(Z, X)
# %%
# draw figure-upper clustering; lower heatmap
draw_dendro_with_heatmap(Z,
                         y,
                         pivot_df,
                         draw_hline=cut_threshold,
                         font_size=1,
                         dendro_show_leaf=True,
                         hmp_xticklabels=True,
                         save_figure=True,
                         save_path=save_path,
                         save_figure_hier=save_figure_hier)

# %% fcluster标记cluster label
# from scipy.cluster.hierarchy import fcluster
#
# f1 = fcluster(Z, cut_threshold, criterion='distance')
# f1_count = Counter(f1)
#
# f1_count = pd.DataFrame(f1_count.items(), columns=["cluster_idx", "gene_nums"])
#
# hier_cluster_ = pd.DataFrame({"gene": y, "cluster": f1, "color": dn['leaves_color_list']})
#
# hier_cluster_.loc[:, "is_marker"] = hier_cluster_["gene"].apply(is_marker, args=(markers,))
# common.save_csv(hier_cluster_, f"{save_path}/{save_cluster_hier}.csv")
# common.save_csv(f1_count, f"{save_path}/{save_cluster_hier}_counts.csv")
# print(f"Clusters number: {f1_count}")
