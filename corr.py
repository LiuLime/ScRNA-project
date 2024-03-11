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


def generate_edge_attributes(x):
    return tuple([x["gene1"], x["gene2"], {"weight": x["cor_pearson.binMean"], "cell_n": x["cell_n"]}])


def generate_network_from_dataframe(df, group_level, fig_save_path, markers):
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

    df.loc[:, "edge_attr"] = df.apply(lambda x: generate_edge_attributes(x), axis=1)
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


@staticmethod
def is_marker(x, markers):
    if x in markers:
        return True
    else:
        return False


def generate_count_from_dataframe(df, group_level, fig_save_path, markers):
    match group_level:
        case "tissue":
            groupby = ["study", "group", "tissue", "gene1"]
        case "cell":
            groupby = ["study", "group", "tissue", "cellType", "gene1"]
        case _:
            raise TypeError("invalid group_level argument, use 'tissue' 'cell'")

    result = df.groupby(by=groupby)["gene2"].count().reset_index()
    result["is_marker"] = result["gene1"].apply(is_marker, markers=markers)
    result.to_csv(f"{fig_save_path}/marker_degree.csv", index=False)
    return result


def create_graphs_from_csv_data(folder_path, file, group_level, fig_save_path, **kwargs):
    """
    param:
        * kwargs: "corr_shrefold", "p_shrefold"
    """
    marker_tissue_corr_data_path = os.path.join(folder_path, file)
    with open(marker_tissue_corr_data_path, "r") as p:
        marker_tissue_corr_data = pd.read_csv(p, sep="\t", header=0)

    # correlation shrefold & p-value shrefold
    marker_tissue_corr_slice = marker_tissue_corr_data[
        marker_tissue_corr_data["cor_pearson.binMean"] > kwargs["corr_shrefold"]]
    marker_tissue_corr_slice = marker_tissue_corr_slice[marker_tissue_corr_slice["p1"] < kwargs["p_shrefold"]]
    marker_tissue_corr_slice = marker_tissue_corr_slice[marker_tissue_corr_slice["p2"] < kwargs["p_shrefold"]]
    graph_list = generate_network_from_dataframe(marker_tissue_corr_slice, group_level=group_level,
                                                 fig_save_path=fig_save_path,
                                                 markers=markers)
    marker_degree = generate_count_from_dataframe(marker_tissue_corr_slice, group_level=group_level,
                                                  fig_save_path=fig_save_path,
                                                  markers=markers)
    return graph_list, marker_degree


# %%
with open("./config.json") as j:
    config = json.load(j)

fig_save_path = config["fig_save_path"]
corrDatafilebymedian = config["corrDataFolderBymedian"]
corrDatafileby4060 = config["corrDataFolderBy40-60"]
corr_shrefold = config["corr_shrefold"]
p_shrefold = config["p_shrefold"]
markers = config["markers"]
TissueMarkerCorrFile = "tissue_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt"
CellMarkerCorrFile = "tissue-cellType_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt"

"""
Aging markers-related correlation&p-value network <- median age
p-value使用的双边检的值
"""

# Tissue level
create_graphs_from_csv_data(corrDatafilebymedian,
                            TissueMarkerCorrFile,
                            group_level="tissue",
                            fig_save_path=fig_save_path + "tissue_level_medianAge_corr0.7_p0.001/",
                            corr_shrefold=corr_shrefold,
                            p_shrefold=p_shrefold,
                            )
# marker_tissue_corr_data_path = os.path.join(corrDatafilebymedian, TissueMarkerCorrFile)
# marker_tissue_corr_data = pd.read_csv(marker_tissue_corr_data_path, sep="\t", header=0)
# # correlation shrefold & p-value shrefold
# marker_tissue_corr_slice = marker_tissue_corr_data[marker_tissue_corr_data["cor_pearson.binMean"] > corr_shrefold]
# marker_tissue_corr_slice = marker_tissue_corr_slice[marker_tissue_corr_slice["p1"] < p_shrefold]
# graph_list = create_df_network(marker_tissue_corr_slice, group_level="tissue", fig_save_path=fig_save_path, markers=markers)

# Cell level
create_graphs_from_csv_data(corrDatafilebymedian,
                            CellMarkerCorrFile,
                            group_level="cell",
                            fig_save_path=fig_save_path + "cell_level_medianAge_corr0.7_p0.001/",
                            corr_shrefold=corr_shrefold,
                            p_shrefold=p_shrefold)

# marker_cell_corr_data_path = os.path.join(corrDatafilebymedian, CellMarkerCorrFile)
# marker_cell_corr_data = pd.read_csv(marker_cell_corr_data_path, sep="\t", header=0)
#
# marker_cell_corr_slice = marker_cell_corr_data[marker_cell_corr_data["cor_pearson.binMean"] > corr_shrefold]
# marker_cell_corr_slice = marker_cell_corr_slice[marker_cell_corr_slice["p1"] < p_shrefold]
# graph_list = create_df_network(marker_cell_corr_slice, group_level="cell", fig_save_path=fig_save_path, markers=markers)

"""
Aging markers-related correlation network <- 40-60 age
"""
# Tissue level
create_graphs_from_csv_data(corrDatafileby4060,
                            TissueMarkerCorrFile,
                            group_level="tissue",
                            fig_save_path=fig_save_path + "tissue_level_40-60Age_corr0.7_p0.001/",
                            corr_shrefold=corr_shrefold,
                            p_shrefold=p_shrefold,
                            )

# marker_tissue_corr_data_path = os.path.join(corrDatafileby4060, TissueMarkerCorrFile)
# marker_tissue_corr_data = pd.read_csv(marker_tissue_corr_data_path, sep="\t", header=0)
# # correlation shrefold
# marker_tissue_corr_slice = marker_tissue_corr_data[marker_tissue_corr_data["cor_pearson.binMean"] > corr_shrefold]
# marker_tissue_corr_slice = marker_tissue_corr_slice[marker_tissue_corr_slice["p1"] < p_shrefold]
# graph_list = create_df_network(marker_tissue_corr_slice, group_level="tissue", fig_save_path=fig_save_path, markers=markers)

# Cell level
create_graphs_from_csv_data(corrDatafileby4060,
                            CellMarkerCorrFile,
                            group_level="cell",
                            fig_save_path=fig_save_path + "cell_level_40-60Age_corr0.7_p0.001/",
                            corr_shrefold=corr_shrefold,
                            p_shrefold=p_shrefold)

# marker_cell_corr_data_path = os.path.join(corrDatafileby4060, CellMarkerCorrFile)
# marker_cell_corr_data = pd.read_csv(marker_cell_corr_data_path, sep="\t", header=0)
#
# marker_cell_corr_slice = marker_cell_corr_data[marker_cell_corr_data["cor_pearson.binMean"] > corr_shrefold]
# marker_cell_corr_slice = marker_cell_corr_slice[marker_cell_corr_slice["p1"] < p_shrefold]
# graph_list = generate_network_from_dataframe(marker_cell_corr_slice, group_level="cell", fig_save_path=fig_save_path,
#                                              markers=markers)

# %%
