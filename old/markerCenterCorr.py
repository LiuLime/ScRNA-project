"""
Step 1:
根据相关性和p-value筛出network

"""

import os

import networkx as nx
import matplotlib.pyplot as plt
import json
import pandas as pd


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
    plt.close()
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


def is_marker(x, markers):
    if x in markers:
        return "yes"
    else:
        return "no"


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


def create_graphs_from_csv_data(folder_path, file, group_level, fig_save_path, markers, **kwargs):
    """
    param:
        * kwargs: "corr_shrefold", "p_shrefold"
    """
    marker_tissue_corr_data_path = os.path.join(folder_path, file)
    with open(marker_tissue_corr_data_path, "r") as p:
        marker_tissue_corr_data = pd.read_csv(p, sep="\t", header=0)

    # correlation shrefold & p-value shrefold
    marker_tissue_corr_slice = marker_tissue_corr_data[
        abs(marker_tissue_corr_data["cor_pearson.binMean"]) > kwargs["corr_shrefold"]]
    marker_tissue_corr_slice = marker_tissue_corr_slice[marker_tissue_corr_slice["p1"] < kwargs["p_shrefold"]]
    marker_tissue_corr_slice = marker_tissue_corr_slice[marker_tissue_corr_slice["p2"] < kwargs["p_shrefold"]]
    graph_list = generate_network_from_dataframe(marker_tissue_corr_slice, group_level=group_level,
                                                 fig_save_path=fig_save_path,
                                                 markers=markers)
    marker_degree = generate_count_from_dataframe(marker_tissue_corr_slice, group_level=group_level,
                                                  fig_save_path=fig_save_path,
                                                  markers=markers)
    return graph_list, marker_degree


def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)


# %%


def main(config, corr_shrefold, p_shrefold):
    fig_save_path = config["fig_save_path"]
    folder_dict = config["folder_dict"]
    corrDatafilebymedian = config["corrDataFolderBymedian"]
    corrDatafileby4060 = config["corrDataFolderBy40-60"]
    TissueMarkerCorrFile = "tissue_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt"
    CellMarkerCorrFile = "tissue-cellType_level/agingMarkers-related_genePairCor_pLT0.2_1ageGrpGE50genes.txt"
    markers = config["markers"]
    """
    Aging markers-related correlation&p-value network <- median age
    p-value使用的双边检的值
    """
    # Tissue level
    path1 = fig_save_path + f"{folder_dict['tm']}corr{corr_shrefold}_p{p_shrefold}/"
    create_folder(path1)
    create_graphs_from_csv_data(corrDatafilebymedian,
                                TissueMarkerCorrFile,
                                group_level="tissue",
                                fig_save_path=path1,
                                markers=markers,
                                corr_shrefold=corr_shrefold,
                                p_shrefold=p_shrefold,
                                )

    # Cell level
    path2 = fig_save_path + f"{folder_dict['cm']}corr{corr_shrefold}_p{p_shrefold}/"
    create_folder(path2)
    create_graphs_from_csv_data(corrDatafilebymedian,
                                CellMarkerCorrFile,
                                group_level="cell",
                                fig_save_path=path2,
                                markers=markers,
                                corr_shrefold=corr_shrefold,
                                p_shrefold=p_shrefold)

    """
    Aging markers-related correlation network <- 40-60 age
    """
    # Tissue level
    path3 = fig_save_path + f"{folder_dict['t46']}corr{corr_shrefold}_p{p_shrefold}/"
    create_folder(path3)
    create_graphs_from_csv_data(corrDatafileby4060,
                                TissueMarkerCorrFile,
                                group_level="tissue",
                                fig_save_path=path3,
                                markers=markers,
                                corr_shrefold=corr_shrefold,
                                p_shrefold=p_shrefold,
                                )

    # Cell level
    path4 = fig_save_path + f"{folder_dict['c46']}corr{corr_shrefold}_p{p_shrefold}/"
    create_folder(path4)
    create_graphs_from_csv_data(corrDatafileby4060,
                                CellMarkerCorrFile,
                                group_level="cell",
                                fig_save_path=path4,
                                markers=markers,
                                corr_shrefold=corr_shrefold,
                                p_shrefold=p_shrefold)


# %%
if __name__ == "__main__":
    with open("./config.json") as j:
        config = json.load(j)
    # single output
    # corr_shrefold = config["corr_shrefold"]
    # p_shrefold = config["p_shrefold"]
    # main(config, corr_shrefold, p_shrefold)
    # multiple output
    corr_shrefold_list = config["corr_shrefold_list"]
    p_shrefold_list = config["p_shrefold_list"]
    for c in corr_shrefold_list:
        for p in p_shrefold_list:
            main(config, c, p)
