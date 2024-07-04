"""
Define gene and tissue pattern by clustering
output:
    - heatmap with dendrogram figure .pdf
    - cluster results .csv

@ 2024/5/31 Liu Yuting | Niigata University
"""

import json

import numpy as np
import pandas as pd
from collections import Counter
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import cophenet, linkage, fcluster, ward, complete, average, centroid, single

from sql import *
from utils import common
from DrawTool import draw_dendro_with_heatmap


def replace_nonzero(x):
    return 1 if x != 0 else 0


def is_marker(x, markers: list):
    """check gene is known marker or not"""
    return True if x in markers else False


def check_cluster_effect(X, metric_list):
    """check the Cophenetic Correlation Coefficient.
    This compares (correlates) the actual pairwise distances of all samples
    to those implied by the hierarchical clustering.
    ---------------------------------------------------
    checked distance function,`ward`, `complete`, `average`,`conetroid`,`single`"""

    result = {}
    for m in metric_list:
        X_pdist = pdist(X, metric=m)
        Z1 = ward(X_pdist)
        Z2 = complete(X_pdist)
        Z3 = average(X_pdist)
        Z4 = centroid(X_pdist)
        Z5 = single(X_pdist)
        c1, coph_dists1 = cophenet(Z1, X_pdist)
        c2, coph_dists2 = cophenet(Z2, X_pdist)
        c3, coph_dists3 = cophenet(Z3, X_pdist)
        c4, coph_dists4 = cophenet(Z4, X_pdist)
        c5, coph_dists5 = cophenet(Z5, X_pdist)
        result[m] = {"ward": c1, "complete": c2, "average": c3, "centroid": c4, "single": c5}
    print(f"Cophenetic Correlation Coefficient:\n{pd.DataFrame(result)}")


with open("./config.json") as j:
    config = json.load(j)
    markers = config["markers"]
    loadPath = config["loadPath"]

db = "scrna_mc"
save_path = loadPath[db]["save_path_clus"]

# 取出log10p top 2000
s = sql(db=db)
with s:
    top_p = pd.DataFrame(s.filter_data_by_rank(), columns=["id", "gene", "log10p", "abbr_id"])

# 转换为gene*tissue 0-1 形式
top_p.loc[:, "is_in"] = top_p["log10p"].apply(replace_nonzero)
pivot_df = top_p.pivot_table(values="is_in", index="gene", columns="abbr_id", fill_value=0)

# remove genes which present less than 10% of  total number of tissue/cell types
total_num_organs = len(pivot_df.columns)
min_num_organs = total_num_organs // 10
# 直接使用sum函数计算每行1的数量并筛选行
pivot_df = pivot_df[pivot_df.sum(axis=1) >= min_num_organs]
print(f"present over {min_num_organs} organ's genes were chosen:", pivot_df.shape)

# %% 测试metric 和method组合
will_check = str(input(f"check cluster effect-yes or no:\n"))
if will_check == "yes":
    metric_list = ["euclidean", "jaccard", "hamming", "yule", "dice"]
    check_cluster_effect(pivot_df.values, metric_list)
    hier_method = str(input(f"chose pdist_metric from `ward`, `complete`, `average`, `centroid`, `single`:"))
    pdist_metric = str(input(f"chose hier_method from {metric_list}:"))
else:
    hier_method = "average"
    pdist_metric = "hamming"
    print(f"clustering default setting were chosen: {hier_method} method, {pdist_metric} metric")

max_pdist_height = 0.8
max_pdist_vertical = 0.8
print(f"""max pdistance for dendrogram : 
                    height line: {max_pdist_height}
                    vertical line:{max_pdist_vertical}""")
save_figure_hier = f"hier_{hier_method}_{pdist_metric}_hmp_hdist{max_pdist_height}_vdist{max_pdist_vertical}_top1000_old"  # top1000还是top100，old还是young，在pickTopP.sql文件里最后一行改。
save_cluster_hier = f"hier_{hier_method}_{pdist_metric}_hmp_dist{max_pdist_vertical}_vdist{max_pdist_vertical}_top1000_old"

# %% 层次聚类
# Perform agglomerative hierarchical clustering
# draw dendrogram + heatmap
Z_col, dn_col, Z_row, dn_row = draw_dendro_with_heatmap(pivot_df,
                                                        hier_method=hier_method,
                                                        pdist_metric=pdist_metric,
                                                        draw_hline=max_pdist_height,
                                                        draw_vline=max_pdist_vertical,
                                                        font_size=1,
                                                        dendro_show_leaf=True,
                                                        hmp_xticklabels=False,
                                                        hmp_yticklabels=False,
                                                        save_figure=True,
                                                        save_path=save_path,
                                                        save_figure_hier=save_figure_hier,
                                                        column_cluster=True,
                                                        row_cluster=True)

# %% fcluster标记cluster label
from scipy.cluster.hierarchy import fcluster


def divide_cluster(linkage_matrix, max_pdist) -> tuple[np.ndarray, pd.DataFrame]:
    f = fcluster(linkage_matrix, max_pdist, criterion='distance')
    f_count = Counter(f)
    f_count = pd.DataFrame(f_count.items(), columns=["cluster_idx", "nums"])
    return f, f_count


will_divide = str(input(f"divide cluster by fcluster - yes or no:\n"))
if will_divide == "yes":

    max_pdist_height = float(input("fcluster max pdistance (float)-height level:"))
    max_pdist_vertical = float(input("fcluster max pdistance (float)-vertical level:"))

    save_fcluster_hier = f"hier_{hier_method}_{pdist_metric}_hdist{max_pdist_height}_vdist{max_pdist_vertical}_top1000_old"

    if Z_col is not None:
        f1, f1_count = divide_cluster(Z_col, max_pdist_height)
        hier_cluster_col = pd.DataFrame({"gene": pivot_df.index, "cluster": f1})  # "color": dn['leaves_color_list']

        hier_cluster_col.loc[:, "is_marker"] = hier_cluster_col["gene"].apply(is_marker, args=(markers,))
        common.save_csv(hier_cluster_col, f"{save_path}/{save_fcluster_hier}_gene.csv")
        common.save_csv(f1_count, f"{save_path}/{save_fcluster_hier}_counts_gene.csv")

    if Z_row is not None:
        f2, f2_count = divide_cluster(Z_row, max_pdist_vertical)
        hier_cluster_row = pd.DataFrame(
            {"organism": pivot_df.columns, "cluster": f2})  # "color": dn['leaves_color_list']

        common.save_csv(hier_cluster_row, f"{save_path}/{save_fcluster_hier}_organ.csv")
        common.save_csv(f2_count, f"{save_path}/{save_fcluster_hier}_counts_organ.csv")
