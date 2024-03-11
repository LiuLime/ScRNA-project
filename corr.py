"""
Step 1:
根据相关性筛出network
根据表达量列出每个细胞中marker相关基因随着marker基因表达量排序的表达量变化

Step 2:
3-D网络图
"""
import os

import networkx as nx
import matplotlib.pyplot as plt
import json
import pandas as pd
import numpy as np


def create_network(nodes, highlight_nodes=None, save_name=None, save_fig=True, save_object=True):
    G = nx.Graph(name=save_name)
    G.add_edges_from(nodes)

    if highlight_nodes:
        node_colors = ["orange" if node in highlight_nodes else "blue" for node in G.nodes()]
        node_labels = {node: node for node in highlight_nodes if node in G.nodes()}
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos=pos, node_color=node_colors, node_size=15, alpha=0.7)
        nx.draw_networkx_edges(G, pos, edge_color="grey")
        nx.draw_networkx_labels(G, pos, labels=node_labels)
    else:
        nx.draw(G, with_labels=True, node_size=15, width=0.5, alpha=0.7)

    if save_fig:
        if save_name:
            save_graph_as_figure(save_name)
        else:
            raise TypeError("missing argument `graph_name`")
    if save_object:
        if save_name:
            save_graph_object(G, save_name)
        else:
            raise TypeError("missing argument `graph_name`")
    return G


def save_graph_as_figure(name, format="pdf"):
    plt.savefig(f"{name}.{format}")


def save_graph_object(G, name):
    nx.write_graphml(G, f"{name}.graphml")


# def graph_laplacian_matrix(G):
#     """calculate graph feature by laplacian matrix"""
#     # 计算拉普拉斯矩阵
#     L = nx.laplacian_matrix(G).toarray()
#     # 计算特征值
#     eigvals = np.linalg.eigvals(L)
#     # 排序特征值
#     eigvals_sorted = np.sort(eigvals)
#     return eigvals_sorted
#
#
# def graph_distance(G1, G2):
#     """计算两个图的距离 graph similarity（distance）"""
#     eigvals1_sorted = graph_laplacian_matrix(G1)
#     eigvals2_sorted = graph_laplacian_matrix(G2)
#     # 比较特征值（这里使用欧式距离计算谱距离）
#     distance = np.linalg.norm(eigvals1_sorted - eigvals2_sorted)
#     return distance


def create_edge_tuple(x):
    return tuple([x["gene1"], x["gene2"], {"weight": x["cor_pearson.binMean"], "cell_n": x["cell_n"]}])


def create_df_network(df, group_level, fig_save_path, markers):
    """
    Create network from dataframe
    param:
        * df: dataframe
        * group_level: "tissue","cell"
        * fig_save_path
        * markers: aging marker list
    """
    match group_level:
        case "tissue":
            groupby = ["study", "group", "tissue"]

        case "cell":
            groupby = ["study", "group", "tissue", "cellType"]
        case _:
            raise TypeError("invalid group_level argument, use 'tissue' 'cell'")

    df.loc[:, "edge_attr"] = df.apply(lambda x: create_edge_tuple(x), axis=1)
    result = df.groupby(by=groupby)["edge_attr"].apply(lambda x: x.tolist()).reset_index()
    result.to_csv(f"{fig_save_path}/temp_result.csv", index=False)
    graph_list = []
    for _, data in result.iterrows():
        if len(groupby) == 3:
            save_name = f"{fig_save_path}/{data['study']}|{data['group']}|{data['tissue']}"
        else:
            save_name = f"{fig_save_path}/{data['study']}|{data['group']}|{data['tissue']}|{data['cellType']}"
        G = create_network(nodes=data["edge_attr"], highlight_nodes=markers, save_name=save_name)
        graph_list.append(G)
    return graph_list


# %%
with open("./config.json") as j:
    config = json.load(j)

fig_save_path = config["fig_save_path"]
corrDatafilebymedian = config["corrDataFolderBymedian"]
corrDatafileby4060 = config["corrDataFolderBy40-60"]
corr_shrefold = config["corr_shrefold"]
markers = config["markers"]
TissueMarkerCorrFile = "tissue_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt"
CellMarkerCorrFile = "tissue-cellType_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt"

"""
Aging markers-related correlation network <- median age
"""

# Tissue level
# marker_tissue_corr_data_path = os.path.join(corrDatafilebymedian, TissueMarkerCorrFile)
# marker_tissue_corr_data = pd.read_csv(marker_tissue_corr_data_path, sep="\t", header=0)
# # correlation shrefold
# marker_tissue_corr_slice = marker_tissue_corr_data[marker_tissue_corr_data["cor_pearson.binMean"] > corr_shrefold]
# graph_list = create_df_network(marker_tissue_corr_slice, group_level="tissue", fig_save_path=fig_save_path, markers=markers)

# Cell level
# marker_cell_corr_data_path = os.path.join(corrDatafilebymedian, CellMarkerCorrFile)
# marker_cell_corr_data = pd.read_csv(marker_cell_corr_data_path, sep="\t", header=0)
#
# marker_cell_corr_slice = marker_cell_corr_data[marker_cell_corr_data["cor_pearson.binMean"] > corr_shrefold]
# graph_list = create_df_network(marker_cell_corr_slice, group_level="cell", fig_save_path=fig_save_path, markers=markers)

"""
Aging markers-related correlation network <- 40-60 age
"""
# Tissue level
# marker_tissue_corr_data_path = os.path.join(corrDatafileby4060, TissueMarkerCorrFile)
# marker_tissue_corr_data = pd.read_csv(marker_tissue_corr_data_path, sep="\t", header=0)
# # correlation shrefold
# marker_tissue_corr_slice = marker_tissue_corr_data[marker_tissue_corr_data["cor_pearson.binMean"] > corr_shrefold]
# graph_list = create_df_network(marker_tissue_corr_slice, group_level="tissue", fig_save_path=fig_save_path, markers=markers)

# Cell level
# marker_cell_corr_data_path = os.path.join(corrDatafileby4060, CellMarkerCorrFile)
# marker_cell_corr_data = pd.read_csv(marker_cell_corr_data_path, sep="\t", header=0)
#
# marker_cell_corr_slice = marker_cell_corr_data[marker_cell_corr_data["cor_pearson.binMean"] > corr_shrefold]
# graph_list = create_df_network(marker_cell_corr_slice, group_level="cell", fig_save_path=fig_save_path, markers=markers)

# %%
import ast
import MySQL

mysql = MySQL.sql()


def query_mysql(query_tissue, marker, target):
    with mysql as sql:
        cell_id = sql.search(table=query_tissue, define_column="gene_name", query_column="cell_id", query=marker)[0]
        marker_exp = \
            sql.search(table=query_tissue, define_column="gene_name", query_column="gene_expression", query=marker)[0]

        target_exp = sql.search(table=query_tissue, define_column=["cell_id", "gene_name"],
                                query_column="gene_expression",
                                query=[cell_id, target])
        if target_exp:
            target_exp = target_exp[0]
    return cell_id, target_exp, marker_exp


test_file = "/Users/liuyuting/WorkSpace/ScRNA_project/03figure/cell_level_medianAge_corr0.7/temp_result.csv"
# test 17.human_TabulaSapiens normal Bladder Fibroblast
with open(test_file, "r") as t:
    data = pd.read_csv(t, header=0)
idx = 11
edges = data.loc[idx, "edge_attr"]
study=data.loc[idx, 'study']
tissue=data.loc[idx, 'tissue']
cellType=data.loc[idx, 'cellType']
print(f"testing {study,tissue , cellType}")
edges = ast.literal_eval(edges)
project_list = []
for tuples in edges:
    target = tuples[0]
    marker = tuples[1]
    cell_id, gene_exp, marker_exp = query_mysql(tissue+"_cell", marker, target)
    projection = {"marker": marker, "target": target, "cell_id": cell_id, "marker_exp": marker_exp,
                  "target_exp": gene_exp}
    project_list.append(projection)
    # print(project_list)
project_df = pd.DataFrame(project_list)
project_df.to_csv(f"{fig_save_path}/exp_result_bladder_fibroblast.csv", index=False)
# %%
