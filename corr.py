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


# G = nx.Graph()
# G.add_edges_from([("A", "B", {"weight": 0.8}), ("B", "C", {"weight": 0.5})])
# G.add_edge("C", "A", weight=0.4, color="red")
# nx.draw(G, with_labels=True)
# plt.savefig("testG2.png")

# %%
def create_network(nodes, save_name=None, save_fig=True, save_object=True):
    G = nx.Graph()
    G.add_edges_from(nodes)
    nx.draw(G, with_labels=True)
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


def save_graph_as_figure(name):
    plt.savefig(f"{name}.png")


def save_graph_object(G, name):
    nx.write_graphml(G, f"{name}.graphml")


def laplacian_matrix():
    # 创建两个图
    G1 = nx.gnp_random_graph(5, 0.5)
    G2 = nx.gnp_random_graph(5, 0.5)

    # 计算拉普拉斯矩阵
    L1 = nx.laplacian_matrix(G1).toarray()
    L2 = nx.laplacian_matrix(G2).toarray()

    # 计算特征值
    eigvals1 = np.linalg.eigvals(L1)
    eigvals2 = np.linalg.eigvals(L2)

    # 排序特征值
    eigvals1_sorted = np.sort(eigvals1)
    eigvals2_sorted = np.sort(eigvals2)

    # 比较特征值（这里简单使用欧式距离作为示例）
    distance = np.linalg.norm(eigvals1_sorted - eigvals2_sorted)

    print(f"特征值1: {eigvals1_sorted}")
    print(f"特征值2: {eigvals2_sorted}")
    print(f"谱距离: {distance}")


def create_edge_tuple(x):
    return tuple([x["gene1"], x["gene2"], {"weight": x["cor_pearson.binMean"]}])


with open("./config.json") as j:
    config = json.load(j)
fig_save_path = config["fig_save_path"]
corrDatafilebymedian = config["corrDataFolderBymedian"]
corr_shrefold = config["corr_shrefold"]
"""
Aging markers-related network using correlation data divided by median age
"""
# Tissue level
marker_tissue_corr_data_path = os.path.join(corrDatafilebymedian,
                                            "tissue_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt")
marker_tissue_corr_data = pd.read_csv(marker_tissue_corr_data_path, sep="\t", header=0)
# correlation shrefold
marker_tissue_corr_slice = marker_tissue_corr_data[marker_tissue_corr_data["cor_pearson.binMean"] > corr_shrefold]
marker_tissue_corr_slice.loc[:, "edge_attr"] = marker_tissue_corr_slice.apply(lambda x: create_edge_tuple(x), axis=1)
result = marker_tissue_corr_slice.groupby(by=["study", "group", "tissue"])["edge_attr"].apply(
    lambda x: x.tolist()).reset_index()

for idx, data in result.iterrows():
    save_name = f"{fig_save_path}/{data['study']}|{data['group']}|{data['tissue']}"
    create_network(nodes=data["edge_attr"], save_name=save_name)

