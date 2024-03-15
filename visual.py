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


def df_arrange_by_tissue(df):
    data_tissue = df.groupby(by=["group", "tissue", "gene1"])["gene2"].sum().reset_index()
    data_tissue["group_tissue"] = data_tissue["group"] + "_" + data_tissue["tissue"]
    new_df = data_tissue.pivot_table(index="group_tissue", columns="gene1", values="gene2", aggfunc="first")

    return new_df


def df_arrange_by_cellType(df):
    data_cell = df.groupby(by=["group", "cellType", "gene1"])["gene2"].mean().reset_index()
    data_cell["group_cell"] = data_cell["group"] + "_" + data_cell["cellType"]
    new_df = data_cell.pivot_table(index="group_cell", columns="gene1", values="gene2", aggfunc="first")

    return new_df


def df_arrange_by_tissueCell(df):
    df["group_tissue_cell"] = df["group"] + "_" + df["tissue"] + df["cellType"]
    new_df = df.pivot_table(index="group_tissue_cell", columns="gene1", values="gene2", aggfunc="first")

    return new_df


# def draw_heatmap(df, title):
#     # 设置图形大小
#     height_per_row = 1
#     width_per_col = 1
#     fig_height = df.shape[0] * height_per_row
#     fig_width = df.shape[1] * width_per_col if df.shape[1] > 3 else 4
#     plt.figure(figsize=(fig_width + 1, fig_height + 1))  #, layout="constrained"
#
#     sns.heatmap(df, annot=True, cmap="YlGnBu", fmt=".0f", square=True)
#     # sns.clustermap(df, cmap="YlGnBu", standard_scale=1)  # row cluster
#
#     plt.rcParams["font.family"] = fig_format["font_family"]
#     plt.title(" Heatmap of marker-related gene number",
#               fontsize=fig_format["title_size"],
#               fontweight="bold")  # 添加标题
#     plt.ylabel("Marker",
#                fontsize=fig_format["font_size"],
#                )  # 设置y轴标签
#     plt.xlabel("Group Tissue",
#                fontsize=fig_format["font_size"],
#                )  # 设置x轴标签
#     plt.tick_params(axis='both',
#                     labelsize=fig_format["label_size"])
#     plt.yticks(rotation=0)
#     plt.xticks(rotation=90)
#     plt.tight_layout()
#     plt.savefig(f"{title}", dpi=300)

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
    ax.set_ylabel("Marker", fontsize=fig_format["font_size"])
    ax.set_xlabel("Group Tissue", fontsize=fig_format["font_size"])
    ax.tick_params(axis='both', labelsize=fig_format["label_size"])

    # 设置轴标签角度
    plt.setp(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.get_yticklabels(), rotation=0)

    # 保存图像时使用bbox_inches='tight'来自动调整边距
    plt.savefig(f"{title}", dpi=300, bbox_inches='tight')


def _excute_tissue(path, file, markers, fmt):
    title = path + f"agingMarkers-related_tissue_heatmap.{fmt}"
    data = pd.read_csv(file, sep=",", header=0)
    data_slice = data[data["is_marker"] == "yes"]

    tissue_df = df_arrange_by_tissue(data_slice)
    tissue_agg_df = add_blank_marker(tissue_df, markers).astype(float).T
    draw_heatmap(tissue_agg_df, title)


def _excute_cell(path, file, markers, fmt):
    title1 = path + f"agingMarkers-related_tissue_heatmap.{fmt}"
    title2 = path + f"agingMarkers-related_cell_heatmap.{fmt}"
    title3 = path + f"agingMarkers-related_tissueCell_heatmap.{fmt}"

    data = pd.read_csv(file, sep=",", header=0)
    data_slice = data[data["is_marker"] == "yes"]

    tissue_df = df_arrange_by_tissue(data_slice)
    tissue_agg_df = add_blank_marker(tissue_df, markers).astype(float).T

    cell_df = df_arrange_by_cellType(data_slice)
    cell_agg_df = add_blank_marker(cell_df, markers).astype(float).T

    tissue_cell_df = df_arrange_by_tissueCell(data_slice)
    tissue_cell_agg_df = add_blank_marker(tissue_cell_df, markers).astype(float).T

    draw_heatmap(tissue_agg_df, title1)
    draw_heatmap(cell_agg_df, title2)
    draw_heatmap(tissue_cell_agg_df, title3)


def excute(path, level, markers, fmt="png"):
    check_folder(path)
    file = path + "marker_degree.csv"

    match level:
        case "tissue":
            _excute_tissue(path, file, markers, fmt)
        case "cell":
            _excute_cell(path, file, markers, fmt)
        case _:
            raise TypeError("missing argument `level`=`tissue` or `cell`")


with open("./config.json") as j:
    config = json.load(j)

fig_save_path = config["fig_save_path"]
folder_dict = config["folder_dict"]
markers = config["markers"]
corr_shrefold = config["corr_shrefold"]
p_shrefold = config["p_shrefold"]
fig_format = config["fig_config"]
fmt = fig_format["fmt"]

"""
Aging markers-related gene number heatmap <- median age
"""
path1 = fig_save_path + f"{folder_dict['tm']}corr{corr_shrefold}_p{p_shrefold}/"
path2 = fig_save_path + f"{folder_dict['cm']}corr{corr_shrefold}_p{p_shrefold}/"
"""
Aging markers-related gene number heatmap <- 40-60 age
"""
path3 = fig_save_path + f"{folder_dict['t46']}corr{corr_shrefold}_p{p_shrefold}/"
path4 = fig_save_path + f"{folder_dict['c46']}corr{corr_shrefold}_p{p_shrefold}/"

excute(path1, level="tissue", markers=markers, fmt=fmt)
excute(path2, level="cell", markers=markers, fmt=fmt)
excute(path3, level="tissue", markers=markers, fmt=fmt)
excute(path4, level="cell", markers=markers, fmt=fmt)
