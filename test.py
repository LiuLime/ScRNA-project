# # %%
# from sqlalchemy import create_engine, Column, Integer, String, ForeignKey
# from sqlalchemy.orm import relationship, sessionmaker, declarative_base
#
# Base = declarative_base()
#
#
# # 定义子表模型
# class Child(Base):
#     __tablename__ = 'child'
#
#     id = Column(Integer, primary_key=True)
#     name = Column(String(30))
#     # 假设有其他字段
#
#
# # 定义总表模型，其中包含指向子表的外键
# class Parent(Base):
#     __tablename__ = 'parent'
#
#     id = Column(Integer, primary_key=True)
#     child_id = Column(Integer, ForeignKey('child.id'))
#
#     # 通过relationship建立与Child模型的联系，这将允许我们方便地访问子表的实例
#     child = relationship("Child")
#
#
# # 创建数据库引擎
# engine = create_engine('mysql+mysqlconnector://root@localhost:3306/test')
#
# # 创建所有表
# Base.metadata.create_all(engine)
#
# # 创建会话
# Session = sessionmaker(bind=engine)
# session = Session()
#
# # 示例：创建一个Child实例和一个与之关联的Parent实例
# # new_child = Child(name ="child 1")
# # new_parent = Parent(child=new_child)
# # session.add(new_parent)
#
# new_child = [Child("Child 1"), Child("Child2")]
# new_parent=Parent(child=new_child)
# session.add_all(new_parent)
#
#
# session.commit()
#
# # 查询示例：获取所有Parent对象及其关联的Child对象
# parents = session.query(Parent).all()
# for parent in parents:
#     print(f'Parent ID: {parent.id}, Child Name: {parent.child.name}')

# %%
from database import db
from sqlalchemy import text, create_engine

# def add_id(table):
#     with engine.connect() as conn:
#         conn.execute(text(f"ALTER TABLE {table} DROP PRIMARY KEY"))
#         conn.execute(text(f"ALTER TABLE {table} ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY FIRST"))
#         conn.commit()
#
#
# db_list = ["scrna_mt", "scrna_mc", "scrna_4060t", "scrna_4060c"]
# for db_path in db_list:
#     engine = create_engine(f'mysql+mysqlconnector://root@localhost:3306/{db_path}')
#     db = db.databaseHelper(db_path)
#     db.init()
#     _, tasks = db.get_available_tasks()
#     print("available_tasks", tasks)
#     # print(db.get_table_cutoff('s8g2t20',0.7,1))
#     for task in tasks:
#         table = task["Table"]
#         add_id(table)

# %%
# from utils.log import logger
# log=logger()
# from database import db
# db_list = ["scrna_mt", "scrna_mc", "scrna_4060t", "scrna_4060c"]
# db = db.databaseHelper(db_list[2])
# db.init()
#
# log.debug(db.get_correlation_records_above_threshold('s8g2t20',0.8))

# %%
# from utils.log import logger
#
# log = logger()
# log.error("test error")
# log.debug("test debug")
#
# #%%
# import pandas as pd
# test = [{'id': 191983, 'gene': 'A2ML1', 'log10p': '4.52676230992907', 'abbr_id': 's3g2t4'},
#         {'id': 192095, 'gene': 'AC024293.1', 'log10p': '25.3048788411024', 'abbr_id': 's3g2t4'}]
#
# data=pd.DataFrame(test)

# %%
# test=[{'id': 2, 'Table': 's12g2t41c25'}, {'id': 3, 'Table': 's3g2t4c13'}]
# print(list(zip(*test[0].items())))
# %%
# import pandas as pd
#
# path = "/Users/liuyuting/WorkSpace/ScRNA_project/mysql/20230328_median_10humanStudies-20230410_40-60_5humanStudies.1ageGrpGE50cells/ageGrp_byMedian/tissue-cellType_level/log10p.csv"
# data = pd.read_csv(path, sep=",",header=0)
#
# count=0
# for idx, row in data.iterrows():
#     if row["log10p"] =='#NAME?':
#         count+=1
# print(count)

# %%
# from MySQL import sql
# import pandas as pd
# test_p=pd.DataFrame([{"gene":"a","log10p":3,"abbr_id":"s1"},
#                      {"gene":"b","log10p":33,"abbr_id":"s2"}])
#
# s=sql('test')
# with s:
#     # tasks=s.get_avalible_studies()
#     # print(tasks)
#     # s.create_log10p_table(test_p)
#     table=s.screen_by_corr_log10p('s12g2t38',1,1)

# %%
import networkx as nx
import matplotlib.pyplot as plt

# 创建一个图
G = nx.Graph()

# 添加节点，可以附带节点属性
G.add_node('A', ratio=0.5)
G.add_node('B', ratio=0.3)
G.add_node('C', ratio=0.2)

# 添加边
G.add_edge('A', 'B')
G.add_edge('B', 'C')
G.add_edge('A', 'C')

# 绘图位置
pos = nx.spring_layout(G)  # 生成一个布局

# 绘制网络图
nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=3000)

# 添加节点比例标签
for node, (x, y) in pos.items():
    ratio = G.nodes[node]['ratio']
    plt.text(x, y + 0.1, s=f'ratio={ratio}', bbox=dict(facecolor='white', alpha=0.5), horizontalalignment='center')

# 显示图形
plt.show()
# %%
import pandas as pd
# pd.options.display.max_columns = 20
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# rng = np.random.RandomState(seed=5)
# ints = rng.randint(1, 11, size=(3, 2))
# a = ["A", "B", "C"]
# b = ["D", "A", "E"]
# df = pd.DataFrame(ints, columns=["weight", "cost"])
# df[0] = a
# df["b"] = b
# G = nx.from_pandas_edgelist(df, 0, "b", ["weight", "cost"])

edges = pd.DataFrame(
    {
        "source": [0, 1, 2],
        "target": [2, 2, 3],
        "weight": [3, 4, 15],
        "color": ["red", "blue", "yellow"],
        "width": [1, 4, 5]

    }
)
G = nx.from_pandas_edgelist(edges, edge_attr=True)
pos = nx.spring_layout(G)
tissue_nodes = edges["source"]
target_nodes = edges["target"]
nx.draw_networkx_nodes(G, pos=pos, nodelist=tissue_nodes, node_size=30, alpha=0.7, node_color="black", node_shape="v")
nx.draw_networkx_nodes(G, pos=pos, nodelist=target_nodes, node_size=10, alpha=0.7, node_color="yellow", node_shape="o")

nx.draw_networkx_edges(G, pos, edge_color=edges["color"], width=edges["width"])
node_labels = {node: node for node in G.nodes}
nx.draw_networkx_labels(G, pos, labels=node_labels)
edge_labels = nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

# %% test pyvis
from pyvis.network import Network
import networkx as nx

nx_graph = nx.cycle_graph(10)
nx_graph.nodes[1]['title'] = 'Number 1'
nx_graph.nodes[1]['group'] = 1
nx_graph.nodes[3]['title'] = 'I belong to a different group!'
nx_graph.nodes[3]['group'] = 10
nx_graph.add_node(20, size=20, title='couple', group=2)
nx_graph.add_node(21, size=15, title='couple', group=2)
nx_graph.add_edge(20, 21, weight=5)
nx_graph.add_node(25, size=25, label='lonely', title='lonely node', group=3)
nt = Network('500px', '500px')
# populates the nodes and edges data structures
nt.from_nx(nx_graph)
nt.show('nx.html')

# %%
import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()
#
# G.add_edge("a", "b", weight={"weight":0.6,"line_wide":1})
# G.add_edge("a", "c", weight={"weight":0.2,"line_wide":2})
# G.add_edge("c", "d", weight={"weight":0.1,"line_wide":1.5})
# G.add_edge("c", "e", weight={"weight":0.7,"line_wide":6})
# G.add_edge("c", "f", weight={"weight":0.9,"line_wide":5})
# G.add_edge("a", "d", weight={"weight":0.3,"line_wide":4})

G.add_edge("a", "b", weight_=0.6, length=1)
G.add_edge("a", "c", weight_=0.2, length=2)
G.add_edge("c", "d", weight_=0.1, length=1.5)
G.add_edge("c", "e", weight_=0.7, length=6)
G.add_edge("c", "f", weight_=0.9, length=5)
G.add_edge("a", "d", weight_=0.3, length=4)
from pyvis.network import Network

nt = Network('500px', '500px', filter_menu=True)
nt.from_nx(G, edge_weight_transf=G.weight_)
nt.show('nx2.html')
# elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight_"] > 0.5]
# esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight_"] <= 0.5]
#
# pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility
#
# # nodes
# nx.draw_networkx_nodes(G, pos, node_size=700)
#
# # edges
# nx.draw_networkx_edges(G, pos, edgelist=elarge, width=6)
# nx.draw_networkx_edges(
#     G, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
# )
#
# # node labels
# nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")
# # edge weight labels
# edge_labels = nx.get_edge_attributes(G, "weight_")
# nx.draw_networkx_edge_labels(G, pos, edge_labels)
#
# ax = plt.gca()
# ax.margins(0.08)
# plt.axis("off")
# plt.tight_layout()
# plt.show()

# %%

import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()

G.add_edge("a", "b", weight=0.6)
G.add_edge("a", "c", weight=0.2)
G.add_edge("c", "d", weight=0.1)
G.add_edge("c", "e", weight=0.7)
G.add_edge("c", "f", weight=0.9)
G.add_edge("a", "d", weight=0.3)

elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > 0.5]
esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]

pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility

# nodes
nx.draw_networkx_nodes(G, pos, node_size=700)

# edges
nx.draw_networkx_edges(G, pos, edgelist=elarge, width=6)
nx.draw_networkx_edges(
    G, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
)

# node labels
nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")
# edge weight labels
edge_labels = nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels)

ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.tight_layout()
plt.show()
# %%
import pandas as pd
import re

data = pd.DataFrame({"gene1": ["ACTB.1", "ACTG1.1", "LINC00936"],
                     "gene2": [3, 5, 6],
                     "is_marker": ["no", "no", "no"],
                     "study": ["s8g2t16", "s8g2t16", "s8g2t16"]})


def map_non_coding_gene(x):
    pattern1 = "\."
    pattern2 = "^LINC\d+"
    pattern3 = "^(RPS|RPL)[0-9]{3,}"
    pattern4 = "^RF[0-9]{3,}"
    flag = False
    if re.search(pattern1, x):
        print(x, 1)
        flag = True
    elif re.search(pattern2, x):
        print(x, 2)
        flag = True
    elif re.search(pattern3, x):
        print(x, 3)
        flag = True
    elif re.search(pattern4,x):
        flag = True
    return flag


non_coding_gene = data[data["gene1"].map(map_non_coding_gene)]
coding_gene = data.drop(non_coding_gene.index)
#%%
import pandas as pd
import numpy as np
a=pd.Series([2,3,5,-6])
np.abs(a).max()

#%%
import os
for root,folder,file in os.walk("./03figure/ageGrp_byMedian/"):
    print(root,folder,file)