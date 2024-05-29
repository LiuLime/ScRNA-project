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
from utils import common, log
import re
import numpy as np
import DrawTool

logger = log.logger()


class preprocess:
    def __init__(self, filepath,
                 non_coding_gene_savetitle="non_coding_gene.csv",
                 fetch_top50_savetitle="degree_top50.csv", ):
        self.path = filepath
        self.non_coding_gene_savetitle = non_coding_gene_savetitle
        self.fetch_top50_savetitle = fetch_top50_savetitle
        self.fetch_marker_savetitle = "degree_marker.csv"
        self.log = log.logger()

    def fetch_marker(self, degree_df: pd.DataFrame) -> pd.DataFrame:
        """fetch marker genes"""
        coding_gene_slice = degree_df[degree_df["is_marker"] == "yes"]
        coding_gene_slice["study"] = coding_gene_slice["abbr_id"].apply(lambda x: re.search("s\d+", x)[0])
        coding_gene_slice["organ"] = coding_gene_slice["abbr_id"].apply(lambda x: re.search("g\d+t\d+(c\d+)?", x)[0])
        top_df = coding_gene_slice.groupby(by=["organ", "gene1", "is_marker", "jaccard_index"])[
            "gene2"].sum().reset_index()
        top_df.columns = ["source", "target", "is_marker", "jaccard_index", "connect"]
        common.save_csv(top_df, os.path.join(self.path, self.fetch_marker_savetitle))
        return top_df

    def fetch_top50(self, degree_df: pd.DataFrame, rank_by="gene2", fetch_by="gene2") -> pd.DataFrame:
        """fetch top 50 genes in each study-organ
        :param: rank_by: 指降序排列的列
        :param: fetch_by: 指后续画图gene1的connect列，是'gene1'的study和或者是'jaccard index'的study和

        input dataframe(columns=['gene1','gene2','jaccard_index','is_marker','abbr_id'])
        return dataframe(columns=["source", "target", "is_marker", "connect"]), save as 'degree_top50.csv'"""
        # 删除RNA gene, non-coding gene and ribosome gene(RPS|RPL series)
        df = self.remove_non_coding_gene(degree_df)

        # 取每个study的top50
        coding_gene = df.sort_values(by=[rank_by], ascending=False)
        coding_gene_slice = coding_gene.groupby(by=["abbr_id"]).apply(preprocess.top, rank_by=rank_by).reset_index(
            drop=True)

        # 分割study和organ
        coding_gene_slice["study"] = coding_gene_slice["abbr_id"].apply(lambda x: re.search("s\d+", x)[0])
        coding_gene_slice["organ"] = coding_gene_slice["abbr_id"].apply(lambda x: re.search("g\d+t\d+(c\d+)?", x)[0])
        top_df = self.fetch_top50_by_(fetch_by, coding_gene_slice)
        return top_df

    @staticmethod
    def map_non_coding_gene(x):
        """删除RNA gene,non-coding gene,ribosome gene和mitochondria"""
        pattern1 = "\."
        pattern2 = "^LINC\d+"
        pattern3 = "^RP(S|L)"  # 删除ribosome gene
        pattern4 = "^RF[0-9]{3,}"
        pattern5 = "^MT-"  # 删除mitochondria gene
        flag = False
        if re.search(pattern1, x):
            flag = True
        elif re.search(pattern2, x):
            flag = True
        elif re.search(pattern3, x):
            flag = True
        elif re.search(pattern4, x):
            flag = True
        elif re.search(pattern5, x):
            flag = True
        return flag

    @staticmethod
    def map_abbr_id(x):
        """match abbreviation id"""
        pattern1 = "^(s\d+)?g\d+t\d+(c\d+)?"
        flag = False
        if re.search(pattern1, x):
            flag = True
        return flag

    @staticmethod
    def top(df: pd.DataFrame, rank_by) -> pd.DataFrame:
        """return top 50 records where the value of `rank_by` is greater than zero.
        If there are fewer than 50 such entries, return all entries where the value is greater than zero."""
        filter_df = df[df[rank_by] > 0]
        return filter_df.sort_values(by=rank_by, ascending=False).head(50)

    def remove_non_coding_gene(self, df):
        non_coding_gene = df[df["gene1"].map(preprocess.map_non_coding_gene)]
        common.save_csv(non_coding_gene, os.path.join(self.path, self.non_coding_gene_savetitle))
        coding_gene = df.drop(non_coding_gene.index)
        return coding_gene

    def fetch_top50_by_(self, by, coding_gene_slice):
        # self.log.debug("fetch top 50 by")

        match by:
            case "gene2":
                # sum same study's gene2
                top_df = coding_gene_slice.groupby(by=["organ", "gene1", "is_marker", "jaccard_index"])[
                    "gene2"].sum().reset_index()
                top_df.columns = ["source", "target", "is_marker", "jaccard_index", "connect"]

            case "jaccard_index":
                # sum same study's jaccard_index
                top_df = coding_gene_slice.groupby(by=["organ", "gene1", "is_marker"])[
                    "jaccard_index"].sum().reset_index()
                top_df.columns = ["source", "target", "is_marker", "connect"]
        common.save_csv(top_df, os.path.join(self.path, self.fetch_top50_savetitle))

        return top_df

    @staticmethod
    def cutoff_by_maxrows(df: pd.DataFrame, sort_columns: str, max_rows=50) -> pd.DataFrame:
        """
        DataFrame 按照指定列降序排列。对于每个 int 值，其组内所有行要么完全包含在结果中，要么完全不包含，以避免截断。
        最终选择的行数接近但不超过max_rows行。
        """
        value_counts = df[sort_columns].value_counts().sort_index(ascending=False)

        # 确定哪些组可以完全包含在输出中
        cumulative_count = 0
        included_groups = []

        for value, count in value_counts.items():
            if cumulative_count + count <= max_rows:
                cumulative_count += count
                included_groups.append(value)
            else:
                break

        # 筛选出可以完全包含的组
        selected_df = df[df[sort_columns].isin(included_groups)]
        return selected_df

    def prepare_file_for_net(self, edges: pd.DataFrame) -> pd.DataFrame:
        """ input dataframe(columns=["source", "target", "target_is_marker", "connect"])
        return dataframe(columns=["source", "target", "target_is_marker", "connect", "normConn"])"""
        # 连接数normalization，作为network的edge weight
        # print(edges.columns)
        max_conn = 2 ** np.ceil(np.log(np.abs(edges["connect"]).max()) / np.log(2))
        edges["normConn"] = edges["connect"].map(lambda x: np.sqrt(np.abs(x) / max_conn))

        return edges

    def prepare_file_for_heatmap(self, node_degree_df: pd.DataFrame,
                                 degree_top50_df: pd.DataFrame,
                                 degree_threshold: int = 5,
                                 max_node_degree_row: int = 30) -> pd.DataFrame:

        """
        node_degree_df: columns=['node','node_degree']
        degree_top50_df: columns=["source", "target", "target_is_marker", "connect", "normConn"]
        """

        # 去掉organ行
        gene_degree_df = node_degree_df.drop(node_degree_df[
                                                 node_degree_df["node"].map(preprocess.map_abbr_id)].index)

        # 根据degree_threshold筛选候选gene list,值越高表示在network中连接的organ越多
        gene_degree_df = gene_degree_df[gene_degree_df["node_degree"] >= degree_threshold]

        # 取某列值累加不超过max行的df切片，避免溢出heatmap边界
        gene_degree_df = preprocess.cutoff_by_maxrows(gene_degree_df, 'node_degree', max_node_degree_row)
        gene_list = gene_degree_df["node"].unique()
        self.log.debug(f"{len(gene_list)},{gene_list}")
        # 根据画heatmap的候选gene list得到degree df的切片
        new_degree_top50_df = degree_top50_df[degree_top50_df["target"].isin(gene_list)]

        # transform abbr_id to full_id
        new_degree_top50_df.loc[:, "source"] = new_degree_top50_df["source"].apply(map_id)
        # print(new_degree_top50_df.head(3))
        # reshape df
        pivot_df = self.reshape_df(new_degree_top50_df)

        return pivot_df

    def reshape_df(self, df) -> pd.DataFrame:
        """dataframe reshaped for drawing heatmap
         re-arraged by **summing** counts with same group, tissue and cellType across studies.

        :param: df: marker-degree dataframe
        :return: new df for heatmap drawing
        """
        new_df = df.pivot_table(index="source", columns="target", values="connect", aggfunc="sum")

        return new_df


def parse_string(input_string):
    pattern = r'^(s\d+)?(g\d+)(t\d+)(c\d+)?$'
    match = re.match(pattern, input_string)

    if match:
        s, g, t, c = match.groups()  # 这将返回一个包含结果的元组，如 ('s1', 'g2', 't2', 'c3') 或 ('s1', 'g2', 't2', None)
        return s, g, t, c
    else:
        return "No match found."


def map_id(abbr_id: str):
    """map abbr id to full id"""
    # study_id = common.read_file("./01datasource/study_id.csv")
    # group_id = common.read_file("./01datasource/group_id.csv")
    # tissue_id = common.read_file("./01datasource/tissue_id.csv")
    # cell_id = common.read_file("./01datasource/cell_id.csv")
    #
    # study_id_dict = study_id.set_index("study_id")["study"].to_dict()
    # group_id_dict = group_id.set_index("group_id")["health_group"].to_dict()
    # tissue_id_dict = tissue_id.set_index("tissue_id")["tissue"].to_dict()
    # cell_id_dict = cell_id.set_index("cell_id")["cell"].to_dict()

    study, group, tissue, cell = parse_string(abbr_id)
    study_full_id = study_id_dict.get(study)
    group_full_id = group_id_dict.get(group)
    tissue_full_id = tissue_id_dict.get(tissue)
    cell_full_id = cell_id_dict.get(cell)
    full_id = "|".join([i for i in [study_full_id, group_full_id, tissue_full_id, cell_full_id] if i is not None])
    return full_id


def main(pre=True, draw_network=True, draw_heatmap=True):
    for group in loadPath.keys():  # group-> 'scrna_mt'...
        logger.debug(f"Start process group {group}")
        group_path = loadPath[group]["save_path"]

        folder_list = [i for i in os.listdir(group_path) if not i.startswith(".")]
        for folder in folder_list:
            save_path = os.path.join(group_path, folder)
            file = common.read_file(
                os.path.join(save_path, "degree.csv"))  # degree file是 gene pair 连接数的直接输出文件，包含所有符合corr和p筛选条件的基因连接数

            # prepare table
            p = preprocess(save_path)
            if pre:
                # 取degree file中每个study的前50个（首先去掉了non-coding gene和ribosome gene）
                top_df = p.fetch_top50(file, rank_by="gene2", fetch_by="gene2")

            # draw network
            if draw_network:
                top_df = common.read_file(os.path.join(save_path, "degree_top50.csv"))
                edges = p.prepare_file_for_net(top_df)
                d = DrawTool.network(save_path)
                d.draw_network(edges, folder)
            if draw_heatmap:
                node_degree_df = common.read_file(os.path.join(save_path, "node_degree.csv"))
                degree_top50_df = common.read_file(os.path.join(save_path, "degree_top50.csv"))
                # print(node_degree_df.head())
                # print(degree_top50_df.head())
                degree_threshold = 3
                max_row = 50
                heatmap_table = p.prepare_file_for_heatmap(node_degree_df,
                                                           degree_top50_df,
                                                           degree_threshold=degree_threshold,
                                                           max_node_degree_row=max_row
                                                           )
                DrawTool.draw_heatmap(heatmap_table,
                                      os.path.join(save_path, "heatmap_top50.png"),
                                      figConfig,
                                      color_marker=True,
                                      marker_list=markers,
                                      fig_title=f"Gene pair connection number(>={degree_threshold} degree top{max_row})")
            logger.debug(f"Finish process {group} {folder}")


if __name__ == "__main__":
    logger.debug("开启客服大门💁🚪🚪")
    with open("./config.json") as j:
        config = json.load(j)

    figConfig = config["figConfig"]
    loadPath = config["loadPath"]
    markers = config["markers"]

    study_id = common.read_file("./01datasource/study_id.csv")
    group_id = common.read_file("./01datasource/group_id.csv")
    tissue_id = common.read_file("./01datasource/tissue_id.csv")
    cell_id = common.read_file("./01datasource/cell_id.csv")

    study_id_dict = study_id.set_index("study_id")["study"].to_dict()
    group_id_dict = group_id.set_index("group_id")["health_group"].to_dict()
    tissue_id_dict = tissue_id.set_index("tissue_id")["tissue"].to_dict()
    cell_id_dict = cell_id.set_index("cell_id")["cell"].to_dict()

    main(pre=False,
         draw_network=True,
         draw_heatmap=False)
    logger.debug("关闭客服大门😏👋👋")

# test_file = read_file("./03figure/ageGrp_byMedian/tissue_level/corr0.9_log10p3/degree.csv")
# path = "./03figure/ageGrp_byMedian/tissue_level/corr0.9_log10p3/"
# p = preprocess(path)
# d = network(path)
# edges = p.process_single_file(test_file)
# d.draw_network(edges, "corr0.9_log10p3")
