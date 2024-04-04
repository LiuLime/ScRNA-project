import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import networkx as nx
import os
import pandas as pd


def draw_heatmap(df, title, fig_format):
    """ Heatmap

    :param: df:dataframe
    :param: title: title with full save path
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

    # 设置标签颜色
    if len(ax.get_xticklabels()) > 22:  # marker has 22
        for label in ax.get_xticklabels()[:22]:
            label.set_color("tomato")
        for label in ax.get_xticklabels()[22:]:
            label.set_color("royalblue")
        # 设置标签注释
        # fig.text(1.02,0.16, 'X-label:', color='black', fontsize=fig_format["font_size"], verticalalignment='center')
        # fig.text(1.02,0.14, '-- Markers', color='tomato', fontsize=fig_format["font_size"], verticalalignment='center')
        # fig.text(1.02,0.12, '-- Candidates', color='royalblue', fontsize=fig_format["font_size"], verticalalignment='center')

        marker_patch = mpatches.Patch(color='tomato', label='Markers')
        candidate_patch = mpatches.Patch(color='royalblue', label='Candidates')
        plt.legend(handles=[marker_patch, candidate_patch], loc='upper right', bbox_to_anchor=(1, 1), frameon=False)

    # 保存图像时使用bbox_inches='tight'来自动调整边距
    plt.savefig(f"{title}", dpi=300, bbox_inches='tight')
    plt.close()


class drawNetwork:
    def __init__(self,
                 level: str,
                 save_folder: str,
                 highlight_nodes: list = None,
                 save_fig: bool = True,
                 save_object: bool = False):

        self.level = level
        self.save_folder = save_folder
        self.highlight_nodes = highlight_nodes
        self.groupby_levels = {"tissue": ["study", "group", "tissue"], "cell": ["study", "group", "tissue", "cellType"]}
        self.groupby = self.groupby_levels.get(self.level, None)

        self.save_fig = save_fig
        self.save_object = save_object

    def draw_network(self, nodes, save_title):
        save_path = os.path.join(self.save_folder, save_title)

        G = nx.Graph(name=save_title)
        G.add_edges_from(nodes)

        if self.highlight_nodes:
            node_colors = ["orange" if node in self.highlight_nodes else "blue" for node in G.nodes()]
            node_labels = {node: node for node in self.highlight_nodes if node in G.nodes()}
            pos = nx.spring_layout(G)
            nx.draw_networkx_nodes(G, pos=pos, node_color=node_colors, node_size=15, alpha=0.7)
            nx.draw_networkx_edges(G, pos, edge_color="grey")
            nx.draw_networkx_labels(G, pos, labels=node_labels)
        else:
            nx.draw(G, with_labels=True, node_size=15, width=0.5, alpha=0.7)

        if self.save_fig:
            self.save_graph_as_figure(save_path)

        if self.save_object:
            self.save_graph_object(G, save_path)

        plt.close()
        return G

    @staticmethod
    def save_graph_as_figure(name, format="pdf"):
        plt.savefig(f"{name}.{format}")

    @staticmethod
    def save_graph_object(G, name):
        nx.write_graphml(G, f"{name}.graphml")

    @staticmethod
    def generate_edge_attributes(x):
        return tuple([x["gene1"], x["gene2"], {"weight": x["cor_pearson.binMean"]}])

    def reshape_matrix_from_df(self, df: pd.DataFrame):
        """
        Reshape dataframe into the format for drawing network, save reshaped df as 'network_matrix.csv'
        :param:
            * df: dataframe
            * markers: aging marker list

        """

        df.loc[:, "edge_attr"] = df.apply(lambda x: drawNetwork.generate_edge_attributes(x), axis=1)
        if not self.groupby:
            raise TypeError("invalid 'level' argument")
        reshaped_matrix = df.groupby(by=self.groupby)["edge_attr"].apply(lambda x: x.tolist()).reset_index()
        reshaped_matrix.to_csv(f"{self.save_folder}/network_matrix.csv", index=False)
        return reshaped_matrix

    def draw_network_from_df(self, df, markers):
        self.highlight_nodes = markers
        reshaped_matrix = self.reshape_matrix_from_df(df)

        # draw network
        for _, data in reshaped_matrix.iterrows():
            save_title = "|".join([data[g] for g in self.groupby])
            self.draw_network(nodes=data["edge_attr"], save_title=save_title)
