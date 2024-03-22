import matplotlib.pyplot as plt
import seaborn as sns


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

    # 设置标签颜色
    if len(ax.get_xticklabels()) > 22:  # marker has 22
        for label in ax.get_xticklabels()[:22]:
            label.set_color("tomato")
        for label in ax.get_xticklabels()[22:]:
            label.set_color("royalblue")

    # 保存图像时使用bbox_inches='tight'来自动调整边距
    plt.savefig(f"{title}", dpi=300, bbox_inches='tight')
    plt.close()
