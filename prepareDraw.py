"""Tools for Visualization
 * heatmap

2024/3/8 @Liu
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import json
import seaborn as sns

global fig_format


def check_folder(path):
    if not os.path.exists(path):
        raise FileNotFoundError


def add_blank_marker(df, markers):
    left_markers = [m for m in markers if m not in df.columns]
    add_df = pd.DataFrame(index=df.index, columns=left_markers)
    agg_df = pd.concat([df, add_df], axis=1)
    agg_df = agg_df[markers]
    return agg_df


def add_blank_tissue(df, organs):
    left_organs = [o for o in organs if o not in df.index]
    add_df = pd.DataFrame(index=left_organs, columns=df.columns)
    agg_df = pd.concat([df, add_df], axis=0)
    return agg_df


def df_arrange_by_tissue(df):
    """same group and cellType would record sum value of them"""

    data_tissue = df.groupby(by=["group", "tissue", "gene1"])["gene2"].sum().reset_index()
    data_tissue["group_tissue"] = data_tissue["tissue"] + "_" + data_tissue["group"]
    new_df = data_tissue.pivot_table(index="group_tissue", columns="gene1", values="gene2", aggfunc="sum")

    return new_df


def df_arrange_by_cellType(df):
    """same group and cellType would record mean value of them"""
    data_cell = df.groupby(by=["group", "cellType", "gene1"])["gene2"].mean().reset_index()
    data_cell["group_cell"] = data_cell["cellType"] + "_" + data_cell["group"]
    new_df = data_cell.pivot_table(index="group_cell", columns="gene1", values="gene2", aggfunc="mean")

    return new_df


def df_arrange_by_tissueCell(df):
    """dataframe reshaped for drawing heatmap
     re-arraged by summing counts with same group, tissue and cellType across studies.

    :param: df: marker-degree dataframe
    :return: new df for heatmap drawing
    """
    df["group_tissue_cell"] = df["tissue"] + "_" + df["cellType"] + "_" + df["group"]
    new_df = df.pivot_table(index="group_tissue_cell", columns="gene1", values="gene2", aggfunc="sum")

    return new_df


def df_arrange_organ_list(tissue_df, level):
    """ return group-tissue/cellType list in correlation analysis, used for presenting full organs for overview.
        "study" was ignored.
    :param: tissue_df: study_tissue_cellType_list dataframe
    :param: level: "tissue" or "cell"

    :return: tissue or cellType list
    """
    match level:
        case "tissue":
            tissue_df["ID"] = tissue_df["tissue"] + "_" + tissue_df["group"]
            tissue_list = tissue_df["ID"].unique().tolist()
        case "cell":
            tissue_df["ID"] = tissue_df["tissue"] + "_" + tissue_df["cellType"] + "_" + tissue_df["group"]
            tissue_list = tissue_df["ID"].unique().tolist()
        case _:
            raise TypeError
    return tissue_list


def draw_heatmap(df, title):
    # 设置图形大小
    height_per_row = 0.5
    width_per_col = 0.5
    fig_height = df.shape[0] * height_per_row
    fig_width = df.shape[1] * width_per_col if df.shape[1] > 3 else 4
    cbar_kws = {"shrink": 0.5}  # color bar缩小到原来的一半
    # 使用plt.subplots来创建图像和轴
    fig, ax = plt.subplots(figsize=(fig_width + 1, fig_height + 1), constrained_layout=True)

    # 绘制heatmap
    sns.heatmap(df, annot=True, cmap="YlGnBu", fmt=".0f", square=True, ax=ax, cbar_kws=cbar_kws)

    # 设置字体
    plt.rcParams["font.family"] = fig_format["font_family"]

    # 设置标题、轴标签和其他参数
    ax.set_title("Heatmap of marker-related gene number", fontsize=fig_format["title_size"], fontweight="bold")
    ax.set_ylabel("Group", fontsize=fig_format["font_size"])
    ax.set_xlabel("Marker", fontsize=fig_format["font_size"])
    ax.tick_params(axis='both', labelsize=fig_format["label_size"])

    # 设置轴标签角度
    plt.setp(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.get_yticklabels(), rotation=0)

    # 保存图像时使用bbox_inches='tight'来自动调整边距
    plt.savefig(f"{title}", dpi=300, bbox_inches='tight')
    plt.close()


def _excute_tissue(file):
    data = pd.read_csv(file, sep=",", header=0)
    data_slice = data[data["is_marker"] == "yes"]
    tissue_df = df_arrange_by_tissue(data_slice).sort_index(axis=0)
    return tissue_df


def _excute_cell(file):
    data = pd.read_csv(file, sep=",", header=0)
    data_slice = data[data["is_marker"] == "yes"]

    tissue_df = df_arrange_by_tissue(data_slice).sort_index(axis=0)
    cell_df = df_arrange_by_cellType(data_slice).sort_index(axis=0)
    tissue_cell_df = df_arrange_by_tissueCell(data_slice).sort_index(axis=0)

    return tissue_df, cell_df, tissue_cell_df


def excute(path, level, title: str | dict, add_blank_markers=False, add_blank_organs=False, **kwargs):
    check_folder(path)
    file = path + "marker_degree.csv"
    match level:
        case "tissue":
            tissue_df = _excute_tissue(file)
            if add_blank_markers:
                tissue_df = add_blank_marker(tissue_df, kwargs["markers"]).astype(float)
            if add_blank_organs:
                tissue_df = add_blank_tissue(tissue_df, kwargs["organs"]).astype(float).sort_index(axis=0)
            draw_heatmap(tissue_df, path + title)

        case "cell":
            tissue_df, cell_df, tissue_cell_df = _excute_cell(file)
            if add_blank_markers:
                tissue_df = add_blank_marker(tissue_df, kwargs["markers"]).astype(float)
                cell_df = add_blank_marker(cell_df, kwargs["markers"]).astype(float)
                tissue_cell_df = add_blank_marker(tissue_cell_df, kwargs["markers"]).astype(float)
            if add_blank_organs:
                tissue_df = add_blank_tissue(tissue_df, kwargs["organs"]).astype(float).sort_index(axis=0)
                cell_df = add_blank_tissue(cell_df, kwargs["organs"]).astype(float).sort_index(axis=0)
                tissue_cell_df = add_blank_tissue(tissue_cell_df, kwargs["organs"]).astype(float).sort_index(axis=0)
            draw_heatmap(tissue_df, path + title["tissue"])
            draw_heatmap(cell_df, path + title["cell"])
            draw_heatmap(tissue_cell_df, path + title["tissue_cell"])

        case _:
            raise TypeError("missing argument `level`=`tissue` or `cell`")


def read_file(path):
    with open(path) as p:
        data = pd.read_csv(p, header=0, sep="\t")
    return data


with open("./config.json") as j:
    config = json.load(j)

fig_save_path = config["fig_save_path"]
folder_dict = config["folder_dict"]
markers = config["markers"]
corrDataFolderBymedian = config["corrDataFolderBymedian"]
corrDataFolderBy4060 = config["corrDataFolderBy40-60"]
tissueFile = config["tissueFile"]
cellTypeFile = config["cellTypeFile"]
corr_shrefold = config["corr_shrefold"]
p_shrefold = config["p_shrefold"]
fig_format = config["fig_config"]
fig_save_title_tissue_level = config["fig_save_title_tissue_level"]
fig_save_title_cellType_level = config["fig_save_title_cellType_level"]

# Aging markers-related gene number heatmap <- median age
path1 = fig_save_path + f"{folder_dict['tm']}corr{corr_shrefold}_p{p_shrefold}/"
path2 = fig_save_path + f"{folder_dict['cm']}corr{corr_shrefold}_p{p_shrefold}/"

# Aging markers-related gene number heatmap <- 40-60 age
path3 = fig_save_path + f"{folder_dict['t46']}corr{corr_shrefold}_p{p_shrefold}/"
path4 = fig_save_path + f"{folder_dict['c46']}corr{corr_shrefold}_p{p_shrefold}/"

"""heatmap with only full markers"""
# excute(path1, level="tissue", markers=markers, title=fig_save_title_tissue_level, add_blank_markers=True)
# excute(path2, level="cell", markers=markers, title=fig_save_title_cellType_level, add_blank_markers=True)
# excute(path3, level="tissue", markers=markers, title=fig_save_title_tissue_level, add_blank_markers=True)
# excute(path4, level="cell", markers=markers, title=fig_save_title_cellType_level, add_blank_markers=True)

"""heatmap with full markers and full organs"""
organ1 = read_file(corrDataFolderBymedian + tissueFile)
organ1_list = df_arrange_organ_list(organ1, level="tissue")

organ2 = read_file(corrDataFolderBymedian + cellTypeFile)
organ2_list = df_arrange_organ_list(organ2, level="cell")

organ3 = read_file(corrDataFolderBy4060 + tissueFile)
organ3_list = df_arrange_organ_list(organ3, level="tissue")

organ4 = read_file(corrDataFolderBy4060 + cellTypeFile)
organ4_list = df_arrange_organ_list(organ4, level="cell")

excute(path1, level="tissue", organs=organ1_list, markers=markers, title=fig_save_title_tissue_level,
       add_blank_markers=True, add_blank_organs=True)
excute(path2, level="cell", organs=organ2_list, markers=markers, title=fig_save_title_cellType_level,
       add_blank_markers=True, add_blank_organs=True)
excute(path3, level="tissue", organs=organ3_list, markers=markers, title=fig_save_title_tissue_level,
       add_blank_markers=True, add_blank_organs=True)
excute(path4, level="cell", organs=organ4_list, markers=markers, title=fig_save_title_cellType_level,
       add_blank_markers=True, add_blank_organs=True)
