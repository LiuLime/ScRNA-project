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
from scipy.cluster.hierarchy import cophenet, fcluster, ward, complete, average, centroid, single

from sql import *
from utils import common
from DrawTool import draw_dendro_with_heatmap, calculate_pdist


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


def reshape(top_p: pd.DataFrame) -> pd.DataFrame:
    """reshape dataframe as gene*tissue format, fill with binary value 0-1
     (if gene is screen out in that tissue (top1000), it is represented as 1 in that tissue)"""

    top_p.loc[:, "is_in"] = top_p["log10p"].apply(replace_nonzero)
    pivot_df = top_p.pivot_table(values="is_in", index="gene", columns="abbr_id", fill_value=0)
    return pivot_df


def chose_universe_genes(pivot_df: pd.DataFrame, over_percent: float = 0.05) -> pd.DataFrame:
    """remove genes which present less than `over_percent` of  total number of tissue/cell types"""
    print(f"top pvalue genes:", pivot_df.shape)
    total_num_organs = len(pivot_df.columns)
    min_num_organs = int(total_num_organs * over_percent)
    # 直接使用sum函数计算每行1的数量并筛选行
    pivot_df = pivot_df[pivot_df.sum(axis=1) >= min_num_organs]
    print(f"present over {min_num_organs} organ's genes were chosen:", pivot_df.shape)
    return pivot_df


def set_max_pdist(max_pdist_height: float, max_pdist_vertical: float, pdist_metric: str, hier_method: str) -> tuple[
    str, str]:

    save_figure_hier = f"hier_{hier_method}_{pdist_metric}_hmp_hdist{max_pdist_height}_vdist{max_pdist_vertical}_top1000_old"  # top1000还是top100，old还是young，在pickTopP.sql文件里最后一行改。
    save_cluster_hier = f"hier_{hier_method}_{pdist_metric}_hmp_hdist{max_pdist_height}_vdist{max_pdist_vertical}_top1000_old"
    return save_figure_hier, save_cluster_hier


# %% 层次聚类
# Perform agglomerative hierarchical clustering
# draw dendrogram + heatmap
#
# Z_col, dn_col, Z_row, dn_row = draw_dendro_with_heatmap(pivot_df,
#                                                         hier_method=hier_method,
#                                                         pdist_metric=pdist_metric,
#                                                         draw_hline=max_pdist_height,
#                                                         draw_vline=max_pdist_vertical,
#                                                         font_size=1,
#                                                         dendro_show_leaf=True,
#                                                         hmp_xticklabels=False,
#                                                         hmp_yticklabels=False,
#                                                         save_figure=True,
#                                                         save_path=save_path,
#                                                         save_figure_hier=save_figure_hier,
#                                                         column_cluster=True,
#                                                         row_cluster=True)

# %% fcluster标记cluster label

def divide_cluster(linkage_matrix, max_pdist, label, search_marker=False) -> tuple[
    pd.DataFrame, pd.DataFrame]:
    f = fcluster(linkage_matrix, max_pdist, criterion='distance')
    f_count = Counter(f)
    f_count = pd.DataFrame(f_count.items(), columns=["cluster_idx", "nums"]).sort_values(by=["cluster_idx"])
    hier_cluster = pd.DataFrame({"leaves": label, "cluster": f})
    if search_marker:
        hier_cluster.loc[:, "is_marker"] = hier_cluster["leaves"].apply(is_marker, args=(markers,))

    return hier_cluster, f_count


# will_divide = str(input(f"divide cluster by fcluster - yes or no:\n"))
# if will_divide == "yes":
#
#     max_pdist_height = float(input("fcluster max pdistance (float)-height level:"))
#     max_pdist_vertical = float(input("fcluster max pdistance (float)-vertical level:"))
#
#     save_fcluster_hier = f"hier_{hier_method}_{pdist_metric}_hdist{max_pdist_height}_vdist{max_pdist_vertical}_top1000_old"
#
#     if Z_col is not None:
#         f1, f1_count = divide_cluster(Z_col, max_pdist_height)
#         hier_cluster_col = pd.DataFrame({"gene": pivot_df.index, "cluster": f1})  # "color": dn['leaves_color_list']
#
#         hier_cluster_col.loc[:, "is_marker"] = hier_cluster_col["gene"].apply(is_marker, args=(markers,))
#         common.save_csv(hier_cluster_col, f"{save_path}/{save_fcluster_hier}_gene.csv")
#         common.save_csv(f1_count, f"{save_path}/{save_fcluster_hier}_counts_gene.csv")
#
#     if Z_row is not None:
#         f2, f2_count = divide_cluster(Z_row, max_pdist_vertical)
#         hier_cluster_row = pd.DataFrame(
#             {"organism": pivot_df.columns, "cluster": f2})  # "color": dn['leaves_color_list']
#
#         common.save_csv(hier_cluster_row, f"{save_path}/{save_fcluster_hier}_organ.csv")
#         common.save_csv(f2_count, f"{save_path}/{save_fcluster_hier}_counts_organ.csv")


def parse_input(instruction) -> tuple[list, list]:
    """parse input indicate"""
    instruction_list = instruction.split(" ")
    instruction1 = []
    instruction2 = []
    for i in instruction_list:
        i1, i2 = i.split("+")
        instruction1.append(i1)
        instruction2.append(i2)
    return instruction1, instruction2


if __name__ == "__main__":
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

    pivot_df = reshape(top_p)
    # remove genes which present less than `over_percent` of  total number of tissue/cell types
    pick_percent = float(input(f"input percent for picking genes over *% organisms (example 0.05 for picking genes over 5% organisms):\n"))
    pivot_df = chose_universe_genes(pivot_df, over_percent=pick_percent)

    # will_check = str(input(f"check cluster effect-yes or no:\n"))
    # if will_check == "yes":
    # 测试metric 和method组合
    pdist_metric_list = ["euclidean", "jaccard", "hamming", "yule", "dice"]
    hier_method_list = ["ward", "complete", "average", "centroid", "single"]
    check_cluster_effect(pivot_df.values, pdist_metric_list)

    # set max distance
    instruct1 = str(input(f"""
    chose pdistance metric from {pdist_metric_list};
    chose hierarchy method from {hier_method_list};
    input combination-> `pdis1+hie1 pdist2+hie2`;
    multiple combo is supported by spliting with space->`euclidean+ward jaccard+average`:\n
    """))
    pdist_metric, hier_method = parse_input(instruct1)
    instruct2 = str(input(f"""
    set max_pdist_height and max_pdist_vertical;
    input combination-> `pdist_height1+pdist_vertical1 pdist_height2+pdist_vertical2`;
    multiple combo is supported by spliting with space->`0.8+0.8 0.6+0.5`:\n
    """))
    max_pdist_height, max_pdist_vertical = parse_input(instruct2)
    max_pdist_height = [float(i) for i in max_pdist_height]
    max_pdist_vertical = [float(i) for i in max_pdist_vertical]

    for pm, hm in zip(pdist_metric, hier_method):
        for disth in max_pdist_height:
            for distv in max_pdist_vertical:
                save_figure_hier, save_cluster_hier = set_max_pdist(disth, distv, pm, hm)
                hmp_col_cluster_sizes, hmp_row_cluster_sizes = None, None
                Z_col, Z_row = calculate_pdist(pivot_df,
                                               pdist_metric=pm,
                                               hier_method=hm,
                                               column_cluster=True,
                                               row_cluster=True)

                if Z_col is not None:
                    hier_cluster_col, f_count_col = divide_cluster(Z_col, disth, label=pivot_df.index,
                                                                   search_marker=True)
                    hmp_col_cluster_sizes = f_count_col["nums"].tolist()
                    common.save_csv(hier_cluster_col, f"{save_path}/{save_cluster_hier}_gene.csv")
                    common.save_csv(f_count_col, f"{save_path}/{save_cluster_hier}_counts_gene.csv")
                if Z_row is not None:
                    hier_cluster_row, f_count_row = divide_cluster(Z_row, distv, label=pivot_df.columns,
                                                                   search_marker=False)
                    c = common.converter()
                    hier_cluster_row.loc[:, "full_id"] = hier_cluster_row["leaves"].apply(c.convert_abbrid_to_fullname)
                    hmp_row_cluster_sizes = f_count_row["nums"].tolist()
                    common.save_csv(hier_cluster_row, f"{save_path}/{save_cluster_hier}_organ.csv")
                    common.save_csv(f_count_row, f"{save_path}/{save_cluster_hier}_counts_organ.csv")

                dn_col, dn_row = draw_dendro_with_heatmap(pivot_df,
                                                          Z_col=Z_col,
                                                          Z_row=Z_row,
                                                          dendro_hline=disth,
                                                          dendro_vline=distv,
                                                          hmp_col_cluster_sizes=hmp_col_cluster_sizes,
                                                          hmp_row_cluster_sizes=hmp_row_cluster_sizes,
                                                          font_size=1,
                                                          dendro_show_leaf=True,
                                                          hmp_xticklabels=False,
                                                          hmp_yticklabels=False,
                                                          save_figure=True,
                                                          save_path=save_path,
                                                          save_figure_hier=save_figure_hier,
                                                          column_cluster=True,
                                                          row_cluster=True)
