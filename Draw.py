"""Tools for Visualization
 * heatmap

2024/3/8 @Liu
"""

import pandas as pd
import os
import json
import DrawTool

global fig_format


def read_file(path) -> pd.DataFrame:
    with open(path) as p:
        data = pd.read_csv(p, header=0, sep="\t")
    return data


def check_folder(path):
    if not os.path.exists(path):
        raise FileNotFoundError


def add_blank_marker(df, markers) -> pd.DataFrame:
    left_markers = [m for m in markers if m not in df.columns]
    add_df = pd.DataFrame(index=df.index, columns=left_markers)
    agg_df = pd.concat([df, add_df], axis=1)
    agg_df = agg_df[markers]
    return agg_df


def add_blank_tissue(df, organs) -> pd.DataFrame:
    left_organs = [o for o in organs if o not in df.index]
    add_df = pd.DataFrame(index=left_organs, columns=df.columns)
    agg_df = pd.concat([df, add_df], axis=0)
    return agg_df


def df_arrange_by_tissue(df) -> pd.DataFrame:
    """dataframe reshaped for drawing heatmap
     re-arraged by **summing** counts with same group and tissue across studies.

    :param: df: marker-degree dataframe
    :return: new df for heatmap drawing
    """
    data_tissue = df.groupby(by=["group", "tissue", "gene1"])["gene2"].sum().reset_index()
    data_tissue.loc[:, "group_tissue"] = data_tissue["tissue"] + "_" + data_tissue["group"]
    new_df = data_tissue.pivot_table(index="group_tissue", columns="gene1", values="gene2", aggfunc="sum")

    return new_df


def df_arrange_by_cellType(df) -> pd.DataFrame:
    """dataframe reshaped for drawing heatmap
     re-arraged by **averaging** counts with same group and cellType across studies.

    :param: df: marker-degree dataframe
    :return: new df for heatmap drawing
    """
    data_cell = df.groupby(by=["group", "cellType", "gene1"])["gene2"].mean().reset_index()
    data_cell.loc[:, "group_cell"] = data_cell["cellType"] + "_" + data_cell["group"]
    new_df = data_cell.pivot_table(index="group_cell", columns="gene1", values="gene2", aggfunc="mean")

    return new_df


def df_arrange_by_tissueCell(df) -> pd.DataFrame:
    """dataframe reshaped for drawing heatmap
     re-arraged by **summing** counts with same group, tissue and cellType across studies.

    :param: df: marker-degree dataframe
    :return: new df for heatmap drawing
    """
    df.loc[:, "group_tissue_cell"] = df["tissue"] + "_" + df["cellType"] + "_" + df["group"]
    new_df = df.pivot_table(index="group_tissue_cell", columns="gene1", values="gene2", aggfunc="sum")

    return new_df


def df_arrange_organ_list(tissue_df, level) -> list:
    """ return group-tissue/cellType list in correlation analysis, used for presenting full organs for
    overview."study" was ignored.
    :param: tissue_df: study_tissue_cellType_list dataframe
    :param: level: "tissue" or "cell"

    :return: tissue or cellType list
    """
    match level:
        case "tissue":
            tissue_df.loc[:, "ID"] = tissue_df["tissue"] + "_" + tissue_df["group"]
            tissue_list = tissue_df["ID"].unique().tolist()
        case "cell":
            tissue_df.loc[:, "ID"] = tissue_df["tissue"] + "_" + tissue_df["cellType"] + "_" + tissue_df["group"]
            tissue_list = tissue_df["ID"].unique().tolist()
        case _:
            raise TypeError
    return tissue_list


def _process_tissue_level_markers(file):
    data = pd.read_csv(file, sep=",", header=0)
    data_slice = data[data["is_marker"] == "yes"]
    tissue_df = df_arrange_by_tissue(data_slice)
    return tissue_df


def _process_cell_level_markers(file):
    data = pd.read_csv(file, sep=",", header=0)
    data_slice = data[data["is_marker"] == "yes"]

    tissue_df = df_arrange_by_tissue(data_slice)
    cell_df = df_arrange_by_cellType(data_slice)
    tissue_cell_df = df_arrange_by_tissueCell(data_slice)

    return tissue_df, cell_df, tissue_cell_df


def _process_tissue_level_candidates(file):
    """return dataframe with marker candidates which has top connections in not known marker groups, dealing with tissue level data"""
    data = pd.read_csv(file, sep=",", header=0)
    # median_value = data[data["is_marker"] == "yes"]["gene2"].max()
    # data_slice = data[data["is_marker"] == "no"].sort_values(by="gene2", ascending=False)
    # data_slice_max = data_slice[data["gene2"] >= median_value]
    data_slice_max = data[data["is_marker"] == "no"].sort_values(by="gene2", ascending=False).head(20)
    tissue_df = df_arrange_by_tissue(data_slice_max)
    return tissue_df


def _process_cell_level_candidates(file):
    """return dataframe with marker candidates which has top connections in not known marker groups, dealing with cell level data"""

    data = pd.read_csv(file, sep=",", header=0)
    # median_value = data[data["is_marker"] == "yes"]["gene2"].max()
    # data_slice = data[data["is_marker"] == "no"].sort_values(by="gene2", ascending=False)
    # data_slice_max = data_slice[data["gene2"] >= median_value]
    data_slice_max = data[data["is_marker"] == "no"].sort_values(by="gene2", ascending=False).head(20)

    tissue_df = df_arrange_by_tissue(data_slice_max)
    cell_df = df_arrange_by_cellType(data_slice_max)
    tissue_cell_df = df_arrange_by_tissueCell(data_slice_max)

    return tissue_df, cell_df, tissue_cell_df


def excute(path, level, title: str | dict,
           process_markers=True,
           add_blank_markers=False,
           add_blank_organs=False,
           process_marker_candidates=False,
           **kwargs):
    check_folder(path)
    file = path + "marker_degree.csv"
    tissue_df = pd.DataFrame()
    cell_df = pd.DataFrame()
    tissue_cell_df = pd.DataFrame()
    match level:
        case "tissue":
            if process_markers:
                tissue_df = _process_tissue_level_markers(file)
            if add_blank_markers:
                tissue_df = add_blank_marker(tissue_df, kwargs["markers"])
            if add_blank_organs:
                tissue_df = add_blank_tissue(tissue_df, kwargs["organs"])
            if process_marker_candidates:
                tissue_df_candidate = _process_tissue_level_candidates(file)
                tissue_df = pd.concat([tissue_df, tissue_df_candidate], axis=1)

            tissue_df = tissue_df.astype(float).sort_index(axis=0)
            DrawTool.draw_heatmap(tissue_df, path + title, fig_format)

        case "cell":
            if process_markers:
                tissue_df, cell_df, tissue_cell_df = _process_cell_level_markers(file)

            if add_blank_markers:
                tissue_df = add_blank_marker(tissue_df, kwargs["markers"])
                cell_df = add_blank_marker(cell_df, kwargs["markers"])
                tissue_cell_df = add_blank_marker(tissue_cell_df, kwargs["markers"])
            if add_blank_organs:
                tissue_df = add_blank_tissue(tissue_df, kwargs["organs"])
                cell_df = add_blank_tissue(cell_df, kwargs["organs"])
                tissue_cell_df = add_blank_tissue(tissue_cell_df, kwargs["organs"])
            if process_marker_candidates:
                tissue_df_candidate, cell_df_candidate, tissue_cell_df_candidate = _process_cell_level_candidates(file)
                tissue_df = pd.concat([tissue_df, tissue_df_candidate], axis=1)
                cell_df = pd.concat([cell_df, cell_df_candidate], axis=1)
                tissue_cell_df = pd.concat([tissue_cell_df, tissue_cell_df_candidate], axis=1)

            tissue_df = tissue_df.astype(float).sort_index(axis=0)
            cell_df = cell_df.astype(float).sort_index(axis=0)
            tissue_cell_df = tissue_cell_df.astype(float).sort_index(axis=0)
            DrawTool.draw_heatmap(tissue_df, path + title["tissue"], fig_format)
            DrawTool.draw_heatmap(cell_df, path + title["cell"], fig_format)
            DrawTool.draw_heatmap(tissue_cell_df, path + title["tissue_cell"], fig_format)

        case _:
            raise TypeError("missing argument `level`=`tissue` or `cell`")


def main(process, config, c, p):
    corr_shrefold = c
    p_shrefold = p

    
    # Aging markers-related gene number heatmap <- median age
    # path1 = fig_save_path + f"{folder_dict['tm']}corr{corr_shrefold}_log10p{p_shrefold}/"
    # path2 = fig_save_path + f"{folder_dict['cm']}corr{corr_shrefold}_log10p{p_shrefold}/"

    # Aging markers-related gene number heatmap <- 40-60 age
    path3 = fig_save_path + f"{folder_dict['t46']}corr{corr_shrefold}_log10p{p_shrefold}/"
    # path4 = fig_save_path + f"{folder_dict['c46']}corr{corr_shrefold}_log10p{p_shrefold}/"

    # organ1 = read_file(corrDataFolderBymedian + tissueFile)
    # organ1_list = df_arrange_organ_list(organ1, level="tissue")

    # organ2 = read_file(corrDataFolderBymedian + cellTypeFile)
    # organ2_list = df_arrange_organ_list(organ2, level="cell")

    organ3 = read_file(corrDataFolderBy4060 + tissueFile)
    organ3_list = df_arrange_organ_list(organ3, level="tissue")

    # organ4 = read_file(corrDataFolderBy4060 + cellTypeFile)
    # organ4_list = df_arrange_organ_list(organ4, level="cell")

    match process:
        case "marker":
            """heatmap with only full markers"""
            # excute(path1, level="tissue", markers=markers, title=fig_save_title_tissue_level, add_blank_markers=True)
            # excute(path2, level="cell", markers=markers, title=fig_save_title_cellType_level, add_blank_markers=True)
            excute(path3, level="tissue", markers=markers, title=fig_save_title_tissue_level, add_blank_markers=True)
            # excute(path4, level="cell", markers=markers, title=fig_save_title_cellType_level, add_blank_markers=True)
        case "organ":
            """heatmap with full markers and full organs"""
            # excute(path1, level="tissue", organs=organ1_list, markers=markers, title=fig_save_title_tissue_level,
            #        add_blank_markers=True, add_blank_organs=True)
            # excute(path2, level="cell", organs=organ2_list, markers=markers, title=fig_save_title_cellType_level,
            #        add_blank_markers=True, add_blank_organs=True)
            excute(path3, level="tissue", organs=organ3_list, markers=markers, title=fig_save_title_tissue_level,
                   add_blank_markers=True, add_blank_organs=True)
            # excute(path4, level="cell", organs=organ4_list, markers=markers, title=fig_save_title_cellType_level,
            #        add_blank_markers=True, add_blank_organs=True)
        case "candidate":
            """heatmap with full markers and full organs and candidate markers"""
            # excute(path1, level="tissue", organs=organ1_list, markers=markers, title=fig_save_title_tissue_level,
            #        add_blank_markers=True, add_blank_organs=True, process_marker_candidates=True)
            # excute(path2, level="cell", organs=organ2_list, markers=markers, title=fig_save_title_cellType_level,
            #        add_blank_markers=True, add_blank_organs=True, process_marker_candidates=True)
            excute(path3, level="tissue", organs=organ3_list, markers=markers, title=fig_save_title_tissue_level,
                   add_blank_markers=True, add_blank_organs=True, process_marker_candidates=True)
            # excute(path4, level="cell", organs=organ4_list, markers=markers, title=fig_save_title_cellType_level,
            #        add_blank_markers=True, add_blank_organs=True, process_marker_candidates=True)
        case "simple_candidate":
            """heatmap with full markers and candidate markers, not full organs"""
            # excute(path1, level="tissue", organs=organ1_list, markers=markers, title=fig_save_title_tissue_level,
            #        add_blank_markers=True, add_blank_organs=False, process_marker_candidates=True)
            # excute(path2, level="cell", organs=organ2_list, markers=markers, title=fig_save_title_cellType_level,
            #        add_blank_markers=True, add_blank_organs=False, process_marker_candidates=True)
            excute(path3, level="tissue", organs=organ3_list, markers=markers, title=fig_save_title_tissue_level,
                   add_blank_markers=True, add_blank_organs=False, process_marker_candidates=True)
            # excute(path4, level="cell", organs=organ4_list, markers=markers, title=fig_save_title_cellType_level,
            #        add_blank_markers=True, add_blank_organs=False, process_marker_candidates=True)
        case _:
            raise TypeError("missing argument, please indicate 'marker','organ','candidate'")


if __name__ == "__main__":
    with open("./config.json") as j:
        config = json.load(j)
    fig_format = config["fig_config"]
    fig_save_path = config["fig_save_path"]
    folder_dict = config["folder_dict"]
    markers = config["markers"]
    corrDataFolderBymedian = config["corrDataFolderBymedian"]
    corrDataFolderBy4060 = config["corrDataFolderBy40-60"]
    tissueFile = config["tissueFile"]
    cellTypeFile = config["cellTypeFile"]

    fig_save_title_tissue_level = config["fig_save_title_tissue_level"]
    fig_save_title_cellType_level = config["fig_save_title_cellType_level"]

    # corr_shrefold_list = config["corr_shrefold_list"]
    # p_shrefold_list = config["p_shrefold_list"]
    # for c in corr_shrefold_list:
    #     for p in p_shrefold_list:
            # main("marker", config, c, p)  # output: {name}_m.png
            # main("organ", config, c, p)  # output {name}_f.png
            # main("candidate", config, c, p)  # output {name}_c.png
            # main("simple_candidate", config, c, p)  # output {name}_mc.png
    corr = config["corr_shrefold"]
    p = config["log10p_abs_shrefold"]
    main("marker", config, corr, p) # output: {name}_m.png
    print("marker sucess---------")
    # main("organ", config, corr, p) # output {name}_o.png
    # print("organ sucess--------")
    # main("simple_candidate", config, corr, p) # output {name}_c.png
    # print("simple candidate sucess-------")
    # main("candidate", config, corr, p)   # output {name}_mc.png
    # print("candidate sucess--------")
