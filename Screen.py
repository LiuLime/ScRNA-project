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

class dataFilter:

    def __init__(self, df, level, corr_shrefold, p_shrefold, save_folder, markers: list = None, cor_col=None,
                 p1_col=None, p2_col=None):
        self.df = df
        self.level = level
        self.corr_shrefold = corr_shrefold
        self.p_shrefold = p_shrefold
        self.save_folder = save_folder
        self.markers = markers
        self.cor_col = "cor_pearson.binMean" if cor_col is None else cor_col
        self.p1_col = "log10p1" if p1_col is None else p1_col
        self.p2_col = "log10p2" if p2_col is None else p2_col
        self.groupby_levels = {"tissue": ["study", "group", "tissue", "gene1"],
                               "cell": ["study", "group", "tissue", "cellType", "gene1"]}
        self.groupby = self.groupby_levels.get(self.level, None)

    # def filter_data_by_criteria(self):
    #     """Filter dataframe by correalation value and pvalue, return dataframe slice"""
    #     df_slice = self.df[abs(self.df[self.cor_col]) > self.corr_shrefold]
    #     df_slice = df_slice[abs(df_slice[self.p1_col]) > self.p_shrefold]
    #     df_slice = df_slice[abs(df_slice[self.p2_col]) > self.p_shrefold]
    #     return df_slice

    def generate_degree(self, df_slice):
        """ Prepare for heatmap.Return connection degree df from dataframe slice which filtered by criteria, save csv as 'degree.csv'
        :param: df_slice
        :return: degree dataframe
        """
        degree = df_slice.groupby(by=self.groupby)["gene2"].count().reset_index()
        degree["is_marker"] = degree["gene1"].apply(dataFilter.is_marker, markers=self.markers)
        degree.to_csv(f"{self.save_folder}/degree.csv", index=False)
        return degree

    @staticmethod
    def is_marker(x, markers):
        if x in markers:
            return "yes"
        else:
            return "no"


def main(path, level, corr, p, save_path, markers):
    with open(path, "r") as c:
        data = pd.read_csv(c, sep=",", header=0)
    F = dataFilter(data, level, corr, p, save_path, markers)
    DN = DrawTool.drawNetwork(level, save_path)

    df_slice = F.filter_data_by_criteria()
    DN.draw_network_from_df(df_slice, markers)
    df_degree = F.generate_degree(df_slice)



def create_graphs_from_csv_data(folder_path, file, group_level, fig_save_path, markers, **kwargs):
    """
    param:
        * kwargs: "corr_shrefold", "p_shrefold"
    """
    marker_tissue_corr_data_path = os.path.join(folder_path, file)

    with open(marker_tissue_corr_data_path, "r") as p:
        marker_tissue_corr_data = pd.read_csv(p, sep=",", header=0)

    graph_list = generate_network_from_dataframe(marker_tissue_corr_slice, group_level=group_level,
                                                 fig_save_path=fig_save_path,
                                                 markers=markers)
    print("graph list sucess")
    marker_degree = generate_count_from_dataframe(marker_tissue_corr_slice, group_level=group_level,
                                                  fig_save_path=fig_save_path,
                                                  markers=markers)
    print("marker degree sucess")
    return graph_list, marker_degree


def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)


def convert_corrMatrix(file):
    "convert single correlation Matrix into long matrix"
    with open(file, "r") as f:
        data = pd.read_csv(f, sep="\t", header=0).set_index("Unnamed: 0")
    df_long = data.stack().reset_index()
    df_long.columns = ['gene1', 'gene2', 'cor_pearson.binMean']

    return df_long


def convert_pvalue(p_file, df):
    with open(p_file, "r") as f:
        data = pd.read_csv(f, sep="\t", header=0).set_index("gene")
    p1 = []
    p2 = []
    for idx, row in df.iterrows():
        # query = ":".join([row["study"], row["group"], row["tissue"], row["cellType"]])
        query = ":".join([row["study"], row["group"], row["tissue"]])

        p1.append(data.loc[row["gene1"], query])
        p2.append(data.loc[row["gene2"], query])
    df.loc[:, "log10p1"] = p1
    df.loc[:, "log10p2"] = p2
    return df


def convert_matrix(file, pvalue_file, level):
    title = file.split("/")[-1]
    print("titleXX", title)
    df_long = convert_corrMatrix(file)
    print("df_longXX", len(df_long))
    match level:
        case "tissue":

            study, group, tissue = us.process_title(title)

            df_long.loc[:, "study"] = study
            df_long.loc[:, "group"] = group
            df_long.loc[:, "tissue"] = tissue
            df_long = convert_pvalue(pvalue_file, df_long)
        case "cell":
            study, group, tissue, cellType = convert_title(title, level="cell")
            df_long.loc[:, "study"] = study
            df_long.loc[:, "group"] = group
            df_long.loc[:, "tissue"] = tissue
            df_long.loc[:, "cellType"] = cellType
            df_long = convert_pvalue(pvalue_file, df_long)
        case _:
            raise TypeError
    return df_long


def main(corr_shrefold, p_shrefold):
    """
    Aging markers-related correlation&p-value network <- median age
    p-value使用的双边检的值
    """
    # Tissue level
    # path1 = fig_save_path + f"{folder_dict['tissue_median_age']}corr{corr_shrefold}_p{p_shrefold}/"
    # create_folder(path1)
    # create_graphs_from_csv_data(corrDatafilebymedian,
    #                             genePairFile,
    #                             group_level="tissue",
    #                             fig_save_path=path1,
    #                             markers=markers,
    #                             corr_shrefold=corr_shrefold,
    #                             p_shrefold=p_shrefold,
    #                             )

    # Cell level
    # folder read path
    # folder_path2 = corrDatafilebymedian + cellPairCorrFolder  # file folder path waiting to process
    # # figure save path
    # save_path2 = fig_save_path + f"{folder_dict['cm']}corr{corr_shrefold}_log10p{p_shrefold}/"
    # print("create folder save_path2")
    # create_folder(save_path2)
    # print("after create folder save_path2")
    # # gene pair merged matrix save path
    # save_pair_path = corrDatafilebymedian + cellTypeLevel

    # file_lists = os.listdir(folder_path2)
    # count = 0
    # for file in file_lists:
    #     file_path = folder_path2 + file
    #     converted_matrix = convert_matrix(file_path, pvalueFile2, level="cell")
    #     new_matrix_list.append(converted_matrix)
    #     count += 1
    #     if count % 20 == 0:
    #         temp_merge_matrix = pd.concat(new_matrix_list)
    #         temp_merge_matrix.to_csv(f"{save_pair_path}/temp_merge_matrix.csv", index=False)
    #         print(f"save at {count},file:{file}")
    # merge_matrix = pd.concat(new_matrix_list)
    # merge_matrix.to_csv(f"{save_pair_path}/merge_matrix.csv", sep=",", index=False)

    # print("in main create_graphs_from_csv_data")
    # create_graphs_from_csv_data(save_pair_path,
    #                             genePairFile,
    #                             group_level="cell",
    #                             fig_save_path=save_path2,
    #                             markers=markers,
    #                             corr_shrefold=corr_shrefold,
    #                             p_shrefold=p_shrefold)
    # print("after create_graphs_from_csv_data")
    """
    Aging markers-related correlation network <- 40-60 age
    """
    # Tissue level
    # folder read path
    folder_path3 = corrDatafileby4060 + tissuePairCorrFolder  # file folder path waiting to process
    # figure save path
    save_path3 = fig_save_path + f"{folder_dict['tissue_40_60_age']}corr{corr_shrefold}_log10p{p_shrefold}/"
    print("create folder save_path3")
    create_folder(save_path3)
    print("after create folder save_path3")
    # gene pair merged matrix save path
    save_pair_path = corrDatafileby4060 + tissueLevel
    # new_matrix_list = []
    # file_lists = os.listdir(folder_path3)
    # count = 0
    # for file in file_lists:
    #     print(f"process {count}")
    #     file_path = folder_path3 + file
    #     print("file path:", file_path)
    #     converted_matrix = convert_matrix(file_path, pvalueFile3, level="tissue")
    #     new_matrix_list.append(converted_matrix)
    #     print(f"matrix_list:{len(new_matrix_list)}")
    #     count += 1
    #     if count % 20 == 0:
    #         temp_merge_matrix = pd.concat(new_matrix_list)
    #         temp_merge_matrix.to_csv(f"{save_pair_path}/temp_merge_matrix.csv", index=False)
    #         print(f"save at {count},file:{file}")
    # merge_matrix = pd.concat(new_matrix_list)

    # merge_matrix.to_csv(f"{save_pair_path}/merge_matrix.csv", sep=",", index=False)

    print("in main create_graphs_from_csv_data")
    create_graphs_from_csv_data(save_pair_path,
                                genePairFile,
                                group_level="tissue",
                                fig_save_path=save_path3,
                                markers=markers,
                                corr_shrefold=corr_shrefold,
                                p_shrefold=p_shrefold)
    print("after create_graphs_from_csv_data")

    # Cell level
    # folder_path4 = corrDatafileby4060 + cellPairCorrFolder  # file folder path waiting to process
    # # figure save path
    # save_path4 = fig_save_path + f"{folder_dict['c46']}corr{corr_shrefold}_log10p{p_shrefold}/"
    # create_folder(save_path4)
    # # gene pair merged matrix save path
    # save_pair_path = corrDatafileby4060 + cellTypeLevel

    # file_lists = os.listdir(folder_path4)
    # new_matrix_list = []
    # count = 0
    # for file in file_lists:
    #     file_path = folder_path4 + file
    #     converted_matrix = convert_matrix(file_path, pvalueFile4, level="cell")
    #     new_matrix_list.append(converted_matrix)
    #     count += 1
    #     if count % 20 == 0:
    #         temp_merge_matrix = pd.concat(new_matrix_list)
    #         temp_merge_matrix.to_csv(f"{save_pair_path}/temp_merge_matrix.csv", index=False)
    #         print(f"save at {count},file:{file}")
    # merge_matrix = pd.concat(new_matrix_list)
    # merge_matrix.to_csv(f"{save_pair_path}/merge_matrix.csv", sep=",", index=False)
    # create_graphs_from_csv_data(save_pair_path,
    #                             genePairFile,
    #                             group_level="cell",
    #                             fig_save_path=save_path4,
    #                             markers=markers,
    #                             corr_shrefold=corr_shrefold,
    #                             p_shrefold=p_shrefold)


# %%
if __name__ == "__main__":
    print("start")
    with open("./config.json") as j:
        config = json.load(j)
    print("open yyy ./config.json")
    fig_save_path = config["fig_save_path"]
    folder_dict = config["folder_dict"]
    corrDatafilebymedian = config["corrDataFolderBymedian"]
    corrDatafileby4060 = config["corrDataFolderBy40-60"]
    tissueLevel = config["tissueLevel"]
    cellTypeLevel = config["cellTypeLevel"]
    tissuePairCorrFolder = config["tissuePairCorrFolder"]
    cellPairCorrFolder = config["cellPairCorrFolder"]
    cellPairPFile = config["cellPairPFile"]
    tissuePairPFile = config["tissuePairPFile"]
    genePairFile = config["savedGenePairFile"]
    markers = config["markers"]
    pvalueFile2 = corrDatafilebymedian + cellPairPFile
    pvalueFile3 = corrDatafileby4060 + tissuePairPFile
    pvalueFile4 = corrDatafileby4060 + cellPairPFile
    # single output
    corr = config["corr_shrefold"]
    p = config["log10p_abs_shrefold"]
    main(corr, p)
    # multiple output
    # corr_shrefold_list = config["corr_shrefold_list"]
    # p_shrefold_list = config["p_shrefold_list"]
    # for c in corr_shrefold_list:
    #     for p in p_shrefold_list:
    #         main(config, c, p)
    print("finish")
