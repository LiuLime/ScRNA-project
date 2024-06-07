"""
read arrow file
根据相关性和p-value筛出dataframe切片
save `degree.csv`

@ Liu
"""

import os
import json
import pandas as pd
from utils import common, log
import pyarrow as pa
import MySQL

logger = log.logger()
with open("./config.json") as c:
    config = json.load(c)

markers = config["markers"]
stringdb = config["stringdb_filepath"]


def read_arrow_to_pd(path):
    """'./scrna_mc/s12g2t14c6.arrow'"""
    schema = pa.schema([
        ('gene1', pa.string()),
        ('log10p1', pa.float64()),
        ('gene2', pa.string()),
        ('log10p2', pa.float64()),
        ('cor_pearson.binMean', pa.float64())
    ])

    with pa.ipc.open_file(path) as reader:
        df = reader.read_pandas()
    # print(df[:5])
    return df


def read_arrow_to_pa(path):
    with pa.memory_map(path, 'r') as source:
        loaded_arrays = pa.ipc.open_file(source).read_all()
    return loaded_arrays


class dataFilter:
    def __init__(self):
        self.corr_threshold = None
        self.p_threshold = None
        self.m = matchRealLink(stringdb)

    def filter_data_by_criteria(self, df, corr_threshold, p_threshold, save_path):
        """Filter dataframe by correalation value and pvalue, return dataframe slice"""
        self.corr_threshold = corr_threshold
        self.p_threshold = p_threshold
        df_slice = df[
            (abs(df["cor_pearson.binMean"]) >= self.corr_threshold) & (
                    abs(df["log10p1"]) >= self.p_threshold) & (
                    abs(df["log10p2"]) >= self.p_threshold)]
        # common.save_csv(df_slice, os.path.join(save_path, "data_p_corr_filter.csv"))
        return df_slice


class matchRealLink:
    def __init__(self, stringdb, real_link_threshold=None):
        """integrate stringdb protein real link info with calculate link info"""
        if real_link_threshold is None:
            self.real_link_threshold = 900
        else:
            self.real_link_threshold = real_link_threshold
        self.filter_link = self.filter_stringdb_combine_score(self.real_link_threshold)
        self.stringdb = stringdb

    def is_realpath(self, df_slice):
        """create column `is_real_path`"""
        df_slice = df_slice.set_index(["gene1", "gene2"])
        filter_link = self.filter_link.set_index(["protein1", "protein2"])
        is_real_path = pd.Series(df_slice.index.isin(filter_link.index))
        df_slice = df_slice.reset_index()
        df_slice["is_real_path"] = is_real_path
        return df_slice

    def generate_realpath_similarity(self, df_slice):
        """evaluate similarity between calculate protein interaction link set with stringdb human protein link set by `jaccard index`

        """
        # filter stringdb protein link by combine score
        filter_link_set = self.filter_link.groupby(by="protein1")["protein2"].apply(set)
        calculate_link_set = df_slice.groupby(by=["gene1"])["gene2"].apply(set)

        # 计算Jaccard指数
        jaccard_scores = {}
        for int_val in calculate_link_set.index:
            set1 = filter_link_set.get(int_val, set())
            set2 = calculate_link_set.get(int_val, set())
            jaccard_scores[int_val] = matchRealLink.jaccard_index(set1, set2)

        # 将Jaccard指数添加到df
        df_slice.loc[:, 'jaccard_index'] = df_slice["gene1"].map(jaccard_scores)

        return df_slice

    def filter_stringdb_combine_score(self, real_link_threshold: int) -> pd.DataFrame:
        """match network link with stringdb real world protein interaction link
        :param: protein_link_threshold: int,range[1,999], means 1%~99.9%
        """

        link = common.read_file(self.stringdb["link"])
        info = common.read_file((self.stringdb["info"]))
        info_map = info.set_index("#string_protein_id")["preferred_name"].to_dict()

        filter_link = link[link["combined_score"] >= real_link_threshold]
        filter_link.loc[:, "protein1"] = filter_link["protein1"].map(info_map)
        filter_link.loc[:, "protein2"] = filter_link["protein2"].map(info_map)
        return filter_link

    @staticmethod
    def jaccard_index(set1: set, set2: set) -> float:
        """计算两个集合的Jaccard 相似度指数"""
        if not set1 or not set2:
            return 0.0

        intersection = set1.intersection(set2)
        union = set1.union(set2)
        # 计算Jaccard指数 交集/并集
        jaccard_index = len(intersection) / len(union)
        return jaccard_index


def is_marker(x, marker_list: list) -> str:
    if x in marker_list:
        return "yes"
    else:
        return "no"


def generate_degree(df_slice, abbr_id, stringdb_, marker_list=None):
    """ Prepare for heatmap.Return connection degree df from dataframe slice which filtered by criteria, save csv as 'degree.csv'
        :param: df_slice
        :return: degree dataframe
        """
    m = matchRealLink(stringdb_)
    df_slice = m.is_realpath(df_slice)  # mark stringdb interaction or not
    df_slice = m.generate_realpath_similarity(df_slice)  # calculate jaccard_index
    df_slice = df_slice[df_slice["is_real_path"]]  # only keep stringdb interaction

    degree = df_slice.groupby(by=["gene1", "jaccard_index"])[
        "gene2"].count().reset_index()
    degree["is_marker"] = degree["gene1"].apply(is_marker, marker_list=marker_list)
    degree["abbr_id"] = abbr_id
    # degree.to_csv(f"{self.save_folder}/degree.csv", index=False)
    return degree


def process_group(load_folder: str,
                  corr_threshold: float,
                  p_threshold: float,

                  save_path: str = None,
                  arrow_list: list = None) -> pd.DataFrame:
    if arrow_list is None:
        arrow_list = os.listdir(load_folder)
        arrow_list = [i for i in arrow_list if i.endswith('.arrow')]

    dfs = []
    d = dataFilter()

    for s in arrow_list:
        df = read_arrow_to_pd(os.path.join(load_folder, s))
        new_df = d.filter_data_by_criteria(df, corr_threshold, p_threshold, save_path)
        new_df = generate_degree(new_df,
                                 abbr_id=s.replace(".arrow", ""),
                                 stringdb_=stringdb,
                                 marker_list=markers)
        dfs.append(new_df)
        logger.debug(f"(process_group) {load_folder.split('/')[-2]} {s} count degree ✅")
    total_dfs = pd.concat(dfs, axis=0, ignore_index=True)

    if save_path:
        common.create_folder(save_path)
        common.save_csv(total_dfs, os.path.join(save_path, "degree.csv"))
        logger.debug(
            f"(process_group) {load_folder.split('/')[-2]}:corr{corr_threshold}_log10p{p_threshold} save degree.csv ✅")

    return total_dfs


def process_group_batch(load_folder: str,
                        save_folder: str,

                        corr_threshold_list: list,
                        p_threshold_list: list,

                        study_list=None) -> None:
    for cor in corr_threshold_list:
        for log10p in p_threshold_list:
            sub_save_folder = save_folder + f"corr{cor}_log10p{log10p}/"
            process_group(load_folder,
                          corr_threshold=cor,
                          p_threshold=log10p,

                          arrow_list=study_list,
                          save_path=sub_save_folder)
            logger.debug(
                f"(process_group_batch) {load_folder.split('/')[-2]}:corr{cor}_log10p{log10p} count finish 🎉🎉🎉～～～")
