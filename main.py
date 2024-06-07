import os
import json
import pandas as pd
from utils import common, log
import pyarrow as pa
from Filter import process_group_batch

logger = log.logger()


def excute_filter():
    """根据相关性和p-value筛出dataframe切片"""
    for group in loadPath.keys():
        process_group_batch(load_folder=loadPath[group]["joined"],
                            save_folder=loadPath[group]["save_path"],

                            corr_threshold_list=corr_cutoffs,
                            p_threshold_list=log10p_abs_cutoffs,
                            )


if __name__ == "__main__":
    logger.debug("开启客服大门🙄🧨🚪")

    with open("./config.json") as c:
        config = json.load(c)

    loadPath = config["loadPath"]
    corr_cutoffs = config["corr_cutoffs"]
    log10p_abs_cutoffs = config["log10p_abs_cutoffs"]

    # Excute filtering and matching stringdb
    excute_filter()

    # test_folder = "./01datasource/joined_table/scrna_4060t/"
    # process_group_batch(load_folder=test_folder,
    #                         save_folder="./02result/ageGrp_by40-60/",
    #                         corr_threshold_list=[0.8],
    #                         p_threshold_list=[3],
    #                        )
    logger.debug("关闭客服大门😊🧑‍🤝‍")
