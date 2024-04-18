"""
Create network for top 50 genes owning most connections in all gene pairs
output:
* degree table according to each study+organ
* network.graphml for top 50 genes in each organ
* network.pdf for top 50 genes in each organ
* heatmap for highly multi-connected organs-genes

2024/4/17
@Liu
"""

import json
import os
import pandas as pd
import networkx as nx
from utils import common, log
import re
import math
import matplotlib.pyplot as plt
import numpy as np

logger = log.logger()


def read_file(path) -> pd.DataFrame:
    with open(path) as p:
        delimiter = common.detect_delimiter(path)
        data = pd.read_csv(p, header=0, sep=delimiter)
    return data


def save_csv(df: pd.DataFrame, path):
    df.to_csv(path, index=False, sep=",")
    logger.debug(f"save dataframe {path} 🍺")


def save_net_as_figure(name, format="pdf"):
    plt.savefig(f"{name}.{format}")
    logger.debug(f"save network {name}.{format} 🍺🍺")


def save_net_as_graphml(G, name):
    nx.write_graphml(G, f"{name}.graphml")
    logger.debug(f"save network {name}.graphml 🍺🍺🍺")


# def parse_string(input_string):
#     pattern = r'(s\d+)(g\d+)(t\d+)(c\d+)?'
#     match = re.match(pattern, input_string)
#
#     if match:
#         # 获取所有匹配的组，这会返回一个包含所有组的元组
#         return match.groups()  # 这将返回一个包含结果的元组，如 ('s1', 'g2', 't2', 'c3') 或 ('s1', 'g2', 't2', None)
#     else:
#         return "No match found."
#

class preprocess:
    def __init__(self, filepath):
        self.path = filepath

    @staticmethod
    def map_non_coding_gene(x):
        """删除RNA gene和non-coding gene"""
        pattern1 = "\."
        pattern2 = "^LINC\d+"
        pattern3 = "^(RPS|RPL)[0-9]{3,}"
        pattern4 = "^RF[0-9]{3,}"
        flag = False
        if re.search(pattern1, x):
            flag = True
        elif re.search(pattern2, x):
            flag = True
        elif re.search(pattern3, x):
            flag = True
        elif re.search(pattern4, x):
            flag = True
        return flag

    def remove_non_coding_gene(self, df):
        non_coding_gene = df[df["gene1"].map(preprocess.map_non_coding_gene)]
        save_csv(non_coding_gene, os.path.join(self.path, "non_coding_gene.csv"))
        coding_gene = df.drop(non_coding_gene.index)
        return coding_gene

    @staticmethod
    def top(df: pd.DataFrame) -> pd.DataFrame:
        """return top 50 records according to connection degree"""
        return df.sort_values(by="gene2", ascending=False).head(50)

    def process_single_file(self, df: pd.DataFrame) -> pd.DataFrame:
        # 删除non-coding gene
        df = self.remove_non_coding_gene(df)

        # 取每个study的top50
        coding_gene = df.sort_values(by=["gene2"], ascending=False)
        coding_gene_slice = coding_gene.groupby(by=["abbr_id"]).apply(preprocess.top).reset_index(drop=True)

        # 分割study和organ
        coding_gene_slice["study"] = coding_gene_slice["abbr_id"].apply(lambda x: re.search("s\d+", x)[0])
        coding_gene_slice["organ"] = coding_gene_slice["abbr_id"].apply(lambda x: re.search("g\d+t\d+(c\d+)?", x)[0])

        # sum same study
        edges = coding_gene_slice.groupby(by=["organ", "gene1", "is_marker"])["gene2"].sum().reset_index()

        # 连接数normalization，作为network的edge weight
        max_conn = 2 ** np.ceil(np.log(np.abs(edges["gene2"]).max()) / np.log(2))
        edges["normConn"] = edges["gene2"].map(lambda x: np.sqrt(np.abs(x) / max_conn))

        # rename column
        edges.columns = ["source", "target", "target_is_marker", "connect", "normConn"]

        # save top 50 degree files
        save_csv(edges, os.path.join(self.path, "degree_top50.csv"))

        return edges


class network:
    def __init__(self, save_path):
        self.path = save_path

    def save_node_degree(self, node_degree: dict) -> pd.DataFrame:
        """generate node degree file"""
        node_degree_df = (pd.DataFrame.from_dict(node_degree, orient="index")
        .reset_index()
        .set_axis(
            ["node", "node_degree"], axis=1))
        save_csv(node_degree_df,
                 os.path.join(self.path, "node_degree.csv"))
        return node_degree_df

    def draw_network(self, edges: pd.DataFrame, net_name):
        """draw network entry, save network.graphml, network.pdf, node_degree.csv"""
        fig = plt.figure(figsize=(40, 50))

        G = nx.from_pandas_edgelist(edges, edge_attr=True)
        # test to add single nodes
        # blank_nodes = ["test1","test2","test3"]
        # G.add_nodes_from(blank_nodes)

        node_degree = dict(G.degree())
        self.save_node_degree(node_degree)

        nx.set_node_attributes(G, node_degree, 'degree_')

        pos = nx.spring_layout(G, seed=10)
        # pos = nx.shell_layout(G)
        # pos = nx.circular_layout(G)
        # pos = nx.kamada_kawai_layout(G)

        # draw nodes
        tissue_nodes = edges["source"]
        target_nodes = edges["target"]

        nx.draw_networkx_nodes(G, pos=pos, nodelist=tissue_nodes, alpha=0.7, node_color="orange",
                               node_shape="o", node_size=500)
        nx.draw_networkx_nodes(G, pos=pos, nodelist=target_nodes, alpha=0.7,
                               node_color=[node_degree[n] for n in target_nodes],
                               node_shape="o", node_size=[node_degree[n] * 100 for n in target_nodes], cmap="YlGn")
        # nx.draw_networkx_nodes(G, pos=pos, nodelist=blank_nodes, alpha=0.7, node_color="red", node_shape="v", node_size=500)

        # draw edges
        nx.draw_networkx_edges(G, pos, edge_color="grey", width=edges["normConn"])

        # draw labels
        node_labels = {node: node for node in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=node_labels, )
        save_net_as_graphml(G, os.path.join(self.path, net_name))
        save_net_as_figure(os.path.join(self.path, net_name))
        plt.close()

class heatmap:
    def __init__(self):
        pass
    def
logger.debug("开启客服大门💁🚪🚪")
with open("./config.json") as j:
    config = json.load(j)

figConfig = config["figConfig"]
loadPath = config["loadPath"]

abbr_id = read_file("./01datasource/abbr_id.csv")
study_id = read_file("./01datasource/study_id.csv")
group_id = read_file("./01datasource/group_id.csv")
tissue_id = read_file("./01datasource/tissue_id.csv")
cell_id = read_file("./01datasource/cell_id.csv")
abbr_id_dict = {
    i["abbr_id"]: {"study": i["study"], "group": i["health_group"], "tissue": i["tissue"], "cellType": i["cellType"]}
    for _, i in abbr_id.iterrows()}

for group in loadPath.keys():  # group-> 'scrna_mt'...
    logger.debug(f"Start process group {group}")
    group_path = loadPath[group]["save_path"]

    folder_list = [i for i in os.listdir(group_path) if not i.startswith(".")]
    for folder in folder_list:
        save_path = os.path.join(group_path, folder)
        file = read_file(os.path.join(save_path, "degree.csv"))
        p = preprocess(save_path)
        d = network(save_path)
        edges = p.process_single_file(file)
        d.draw_network(edges, folder)
        logger.debug(f"Finish process {group} {folder}")

logger.debug("关闭客服大门😏👋👋")
# test_file = read_file("./03figure/ageGrp_byMedian/tissue_level/corr0.9_log10p3/degree.csv")
# path = "./03figure/ageGrp_byMedian/tissue_level/corr0.9_log10p3/"
# p = preprocess(path)
# d = network(path)
# edges = p.process_single_file(test_file)
# d.draw_network(edges, "corr0.9_log10p3")
