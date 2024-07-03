from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette

import networkx as nx
import os
import pandas as pd
import numpy as np
from utils import log, common


def draw_heatmap(df,
                 save_title: str = "heatmap.png",
                 fig_format: dict | None = None,
                 color_marker: bool = False,
                 marker_list: list | None = None,
                 fig_title: str = "Heatmap of gene pair connection number"):
    """ Heatmap
    if color_marker is True, please pass marker_list argument, otherwise Userwarning will raise up.
    :param: df:dataframe
    :param: title: title with full save path and format
    :param: color_marker: color x stick label or not
    :param: marker_list: the list of color x stick label
    :param: fig_title: figure title in heatmap.
    """
    # 设置图形大小
    height_per_row = 0.5
    width_per_col = 0.5
    fig_height = df.shape[0] * height_per_row
    fig_width = df.shape[1] * width_per_col if df.shape[1] > 3 else 4
    cbar_kws = {"shrink": 0.3}  # color bar缩小到原来的一半
    # 使用plt.subplots来创建图像和轴
    fig, ax = plt.subplots(figsize=(fig_width + 1, fig_height + 1), constrained_layout=True)

    # 绘制heatmap
    sns.heatmap(df, annot=True, cmap="YlGnBu", fmt=fig_format["fmt"], square=True, ax=ax, cbar_kws=cbar_kws)

    # 设置字体
    plt.rcParams["font.family"] = fig_format["font_family"]

    # 设置标题、轴标签和其他参数
    ax.set_title(fig_title, fontsize=fig_format["title_size"], fontweight="bold")
    ax.set_ylabel("Group", fontsize=fig_format["font_size"])
    ax.set_xlabel("Candidate Marker", fontsize=fig_format["font_size"])
    ax.tick_params(axis='both', labelsize=fig_format["label_size"])

    # 设置轴标签角度
    plt.setp(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.get_yticklabels(), rotation=0)

    # 设置标签颜色
    if color_marker:
        if not marker_list:
            raise UserWarning("no marker list was defined")
        for label in ax.get_xticklabels():
            if label in marker_list:
                label.set_color("tomato")
        marker_patch = mpatches.Patch(color='tomato', label='Known markers')
        plt.legend(handles=[marker_patch], loc='upper right', bbox_to_anchor=(1.1, 1.05), frameon=False)

        # if len(ax.get_xticklabels()) > 22:  # marker has 22
        #     for label in ax.get_xticklabels()[:22]:
        #         label.set_color("tomato")
        #     for label in ax.get_xticklabels()[22:]:
        #         label.set_color("royalblue")
        #     marker_patch = mpatches.Patch(color='tomato', label='Markers')
        #     candidate_patch = mpatches.Patch(color='royalblue', label='Candidates')
        #     plt.legend(handles=[marker_patch, candidate_patch], loc='upper right', bbox_to_anchor=(1, 1), frameon=False)

    # 保存图像时使用bbox_inches='tight'来自动调整边距
    plt.savefig(f"{save_title}", dpi=300, bbox_inches='tight')
    plt.close()


def _dendro(Z, y: np.ndarray | None = None, no_plot: bool = True, show_leaf_counts: bool = True, leaf_font_size=3,
            draw_hline: int | None = None, ax: Any = None):
    dendro = dendrogram(Z,
                        orientation='top',
                        labels=y,
                        ax=ax,
                        leaf_rotation=90,
                        leaf_font_size=leaf_font_size,
                        color_threshold=draw_hline,
                        distance_sort='descending',
                        show_leaf_counts=show_leaf_counts,
                        no_plot=no_plot, )
    return dendro


def draw_dendrogram(Z,
                    y: np.ndarray,
                    draw_hline: int,
                    no_plot: bool = True,
                    save_figure: bool = False,
                    save_path: str = None,
                    save_figure_hier: str = None):
    colors = cm.rainbow(np.linspace(0, 1, draw_hline))
    set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in colors])

    # scipy绘制树状图
    plt.figure(figsize=(120, 70), constrained_layout=True)
    dendro = _dendro(Z, y, draw_hline=draw_hline, no_plot=no_plot)
    plt.axhline(draw_hline, color='k')

    if save_figure:
        plt.savefig(f"{save_path}{save_figure_hier}.pdf")

    plt.show()
    set_link_color_palette(None)
    plt.close()
    return dendro


def draw_dendro_with_heatmap(Z,
                             y,
                             pivot_df,
                             draw_hline: int,
                             font_size: int | None = 3,
                             dendro_show_leaf: bool = True,
                             hmp_xticklabels: bool = True,
                             save_figure: bool = False,
                             save_path: str = None,
                             save_figure_hier: str = None, ):
    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(100, 8),
    #                                gridspec_kw={'height_ratios': [1, 3]}) #加了sharex=True的话dendro会显示不全
    fig = plt.figure(figsize=(120, 70), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=(1, 1), left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    colors = cm.rainbow(np.linspace(0, 1, draw_hline))
    set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in colors])

    dn = dendrogram(Z, ax=ax1, labels=y.to_numpy(),
                    leaf_rotation=90,
                    show_leaf_counts=dendro_show_leaf,
                    leaf_font_size=font_size,
                    color_threshold=draw_hline, )
    ax1.axhline(draw_hline, color='k')
    sorted_idx = dn['leaves']
    sorted_df = pivot_df.iloc[sorted_idx]
    hmp_data = sorted_df.T
    sns.heatmap(hmp_data, ax=ax2, cmap='viridis', cbar=False, xticklabels=hmp_xticklabels)

    if save_figure:
        plt.savefig(f"{save_path}{save_figure_hier}.pdf", bbox_inches='tight')
    plt.tight_layout()
    plt.show()
    set_link_color_palette(None)
    plt.close()


class network:
    def __init__(self, save_path):
        self.path = save_path
        self.log = log.logger()

    def save_node_degree(self, node_degree: dict) -> pd.DataFrame:
        """generate node degree file"""
        node_degree_df = (pd.DataFrame.from_dict(node_degree, orient="index")
        .reset_index()
        .set_axis(
            ["node", "node_degree"], axis=1))
        common.save_csv(node_degree_df,
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
        # nx.draw_networkx_nodes(G, pos=pos, nodelist=target_nodes, alpha=0.7,
        #                        node_color=[node_degree[n] for n in target_nodes],
        #                        node_shape="o", node_size=[node_degree[n] * 100 for n in target_nodes], cmap="YlGn")
        nx.draw_networkx_nodes(G, pos=pos, nodelist=target_nodes, alpha=0.7,
                               node_color=[node_degree[n] for n in target_nodes],
                               node_shape="o", node_size=300, cmap="YlGn")
        # nx.draw_networkx_nodes(G, pos=pos, nodelist=blank_nodes, alpha=0.7, node_color="red", node_shape="v", node_size=500)

        # draw edges
        nx.draw_networkx_edges(G, pos, edge_color="grey", width=edges["normConn"])

        # draw labels
        # node_labels = {node: node for node in G.nodes}
        # nx.draw_networkx_labels(G, pos, labels=node_labels, )

        # draw tissue node labels
        tissue_node_labels = {n: n for n in tissue_nodes}
        nx.draw_networkx_labels(G, pos, labels=tissue_node_labels, )

        # draw top 3 target node labels
        # value_list = node_degree.values()
        # threshold = sorted(value_list)[-5:]
        target_node_labels = {n: n for n in target_nodes if node_degree[n] >= 3}
        nx.draw_networkx_labels(G, pos, labels=target_node_labels, )

        self.save_net_as_graphml(G, os.path.join(self.path, net_name))
        self.save_net_as_figure(os.path.join(self.path, net_name))
        plt.close()

    def save_net_as_figure(self, name, format="pdf"):
        plt.savefig(f"{name}.{format}")
        self.log.debug(f"save network {name}.{format} 🍺🍺")

    def save_net_as_graphml(self, G, name):
        nx.write_graphml(G, f"{name}.graphml")
        self.log.debug(f"save network {name}.graphml 🍺🍺🍺")

# %%
