"""
根据相关性和p-value筛出dataframe切片，并进行可视化

"""

import os
import json
import pandas as pd
import utils as us
import DrawTool
from MySQL import sql
from utils import common
from utils.log import logger
import pyarrow.compute as pc
import pyarrow as pa
import time


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
    def __init__(self, markers):
        self.corr_threshold = None
        self.p_threshold = None
        self.markers = markers
        self.log = logger()

    def filter_data_by_criteria(self, df, corr_threshold, p_threshold, ):
        """Filter dataframe by correalation value and pvalue, return dataframe slice"""
        self.corr_threshold = corr_threshold
        self.p_threshold = p_threshold
        df_slice = df[
            (abs(df["cor_pearson.binMean"]) >= self.corr_threshold) & (
                    abs(df["log10p1"]) >= self.p_threshold) & (
                    abs(df["log10p2"]) >= self.p_threshold)]

        return df_slice

    def generate_degree(self, df_slice, abbr_id):
        """ Prepare for heatmap.Return connection degree df from dataframe slice which filtered by criteria, save csv as 'degree.csv'
        :param: df_slice
        :return: degree dataframe
        """
        degree = df_slice.groupby(by=["gene1"])[
            "gene2"].count().reset_index()
        degree["is_marker"] = degree["gene1"].apply(dataFilter.is_marker, markers=self.markers)
        degree["abbr_id"] = abbr_id
        # degree.to_csv(f"{self.save_folder}/degree.csv", index=False)
        return degree

    @staticmethod
    def is_marker(x, markers: list) -> str:
        if x in markers:
            return "yes"
        else:
            return "no"


def save_csv(save_path: str, df: pd.DataFrame) -> None:
    df.to_csv(f"{save_path}/degree.csv", index=False)


def process_group(load_folder: str,
                  corr_threshold: float,
                  p_threshold: float,
                  markers: list,
                  save_path: str = None,
                  arrow_list: list = None) -> pd.DataFrame:

    if arrow_list is None:
        arrow_list = os.listdir(load_folder)

    dfs = []
    d = dataFilter(markers)
    for s in arrow_list:
        df = read_arrow_to_pd(os.path.join(load_folder, s))
        new_df = d.filter_data_by_criteria(df, corr_threshold, p_threshold, )
        new_df = d.generate_degree(new_df, abbr_id=s.replace(".arrow", ""))
        dfs.append(new_df)
        logger().debug(f"(process_group) {load_folder.split('/')[-2]} {s} count degree ✅")
    total_dfs = pd.concat(dfs, axis=0, ignore_index=True)

    if save_path:
        common.create_folder(save_path)
        save_csv(save_path, total_dfs)
        logger().debug(
            f"(process_group) {load_folder.split('/')[-2]}:corr{corr_threshold}_log10p{p_threshold} save degree.csv ✅")

    return total_dfs


def process_group_batch(load_folder: str,
                        save_folder: str,
                        corr_threshold_list: list,
                        p_threshold_list: list,
                        markers: list,
                        study_list=None) -> None:

    for cor in corr_threshold_list:
        for log10p in p_threshold_list:
            sub_save_folder = save_folder + f"corr{cor}_log10p{log10p}/"
            process_group(load_folder,
                          corr_threshold=cor,
                          p_threshold=log10p,
                          markers=markers,
                          arrow_list=study_list,
                          save_path=sub_save_folder)
            logger().debug(
                f"(process_group_batch) {load_folder.split('/')[-2]}:corr{cor}_log10p{log10p} count finish 🎉🎉🎉～～～")


if __name__ == "__main__":
    logger().debug("开启客服大门🙄🧨🚪")

    with open("./config.json") as c:
        config = json.load(c)

    loadPath = config["loadPath"]
    corr_cutoffs = config["corr_cutoffs"]
    log10p_abs_cutoffs = config["log10p_abs_cutoffs"]
    markers = config["markers"]

    for group in loadPath.keys():
        process_group_batch(load_folder=loadPath[group]["joined"],
                            save_folder=loadPath[group]["save_path"],
                            corr_threshold_list=corr_cutoffs,
                            p_threshold_list=log10p_abs_cutoffs,
                            markers=markers)

    # test_folder = "./01datasource/joined_table/scrna_4060t/"
    # process_group_batch(load_folder=test_folder,
    #                         save_folder="./02result/ageGrp_by40-60/",
    #                         corr_threshold_list=[0.8],
    #                         p_threshold_list=[3],
    #                         markers=markers)
    logger().debug("关闭客服大门😊🧑‍🤝‍")
