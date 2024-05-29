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
# import pandas as pd
# import re
#
# data = pd.DataFrame({"gene1": ["ACTB.1", "ACTG1.1", "LINC00936"],
#                      "gene2": [3, 6, 6],
#                      "is_marker": ["no", "no", "no"],
#                      "study": ["s8g2t16", "s8g2t16", "s8g2t16"]})
# gene_list = pd.Series(["ACTG1.1", "LINC00936"])
# new_data = data[data["gene1"].isin(gene_list)]
#
#
# def map_non_coding_gene(x):
#     pattern1 = "\."
#     pattern2 = "^LINC\d+"
#     pattern3 = "^(RPS|RPL)[0-9]{3,}"
#     pattern4 = "^RF[0-9]{3,}"
#     flag = False
#     if re.search(pattern1, x):
#         print(x, 1)
#         flag = True
#     elif re.search(pattern2, x):
#         print(x, 2)
#         flag = True
#     elif re.search(pattern3, x):
#         print(x, 3)
#         flag = True
#     elif re.search(pattern4, x):
#         flag = True
#     return flag


# non_coding_gene = data[data["gene1"].map(map_non_coding_gene)]
# coding_gene = data.drop(non_coding_gene.index)


# %%
# import re
#
# map_dict = {"s1": "a", "g1": "h", "t1": "eye", "c1": "cell"}
#
#
# def parse_string(input_string):
#     pattern = r'^(s\d+)?(g\d+)(t\d+)(c\d+)?$'
#     match = re.match(pattern, input_string)
#
#     if match:
#         s, g, t, c = match.groups()  # 这将返回一个包含结果的元组，如 ('s1', 'g2', 't2', 'c3') 或 ('s1', 'g2', 't2', None)
#         print(s, g, t, c)
#         full_id = "|".join([map_dict.get(s), map_dict.get(g), map_dict.get(t), map_dict.get(c)])
#
#         return full_id
#     else:
#         return "No match found."
#
#
# test = pd.Series(["s1g1t1c1", "s1g1t1", "g1t1c1"])
# test.applymap(parse_string, args=(map_dict))
# %%
# import pandas as pd
# import drawNetwork
#
# p = drawNetwork.preprocess(filepath=None)
# df = pd.DataFrame({"node": ["ACTB.1", "ACTG1.1", "LINC00936", "AAA", "BBB"],
#                    "node_degree": [3, 6, 6, 6, 2]})
# # count=df.groupby(by=["node_degree"]).count().sort_values(by="node_degree",ascending=False)
# # max_rows=4
# # # df["cumulative_sum"]=count["node"].cumsum()
# # count["cumulative_sum"]=count["node"].cumsum()
#
# # cutoff=count[count['cumulative_sum'] <= max_rows]['cumulative_sum'].max()
# # selected_df = df[df.groupby("node_degree")['cumulative_sum'].transform('max') <= cutoff]
#
# p.cutoff_by_maxrows(df, sort_columns="node_degree", max_rows=4)
# # %%
# import pandas as pd
#
# # 创建示例数据
# data = {
#     'int_column': [1, 2, 2, 3, 3, 3, 4, 4, 4, 4]
# }
# df = pd.DataFrame(data)
#
# # 对整数列进行统计
# value_counts = df['int_column'].value_counts().sort_index(ascending=False)
#
# # 确定哪些组可以完全包含在输出中
# cumulative_count = 0
# included_groups = []
#
# for value, count in value_counts.items():
#     if cumulative_count + count <= 5:  # 例如这里的5可以替换为100或其他限制
#         cumulative_count += count
#         included_groups.append(value)
#     else:
#         break
#
# # 筛选出可以完全包含的组
# result = df[df['int_column'].isin(included_groups)]
#
# print(result)
#
# # %%
# import pandas as pd
#
# def jaccard_index(set1: set, set2: set) -> float:
#     """计算两个集合的Jaccard 相似度指数"""
#     if not set1 or not set2:
#         return 0.0
#     # 计算两个集合的交集
#     intersection = set1.intersection(set2)
#     # 计算两个集合的并集
#     union = set1.union(set2)
#     # 计算Jaccard指数
#     jaccard_index = len(intersection) / len(union)
#     return jaccard_index
#
# # 创建示例数据
# data1 = {
#     'int_column': [1, 4, 4, 4, 2, 2, 3, 3, 4, 4],
#     'test': ["a", "D", "d", "D", "B", "B", "C", "C", "C", "D"]
# }
# data2={'int_column': [2, 4, 5, 4, 2, 2, 3, 5, 3, 4],
# 'test': ["a", "D", "d", "D", "B", "B", "C", "C", "C", "D"]
# }
# df1 = pd.DataFrame(data1)
# df2 = pd.DataFrame(data2)
#
# # print(df[df['int_column'] > 3])
# test_df1 = df1.groupby(by="int_column")["test"].apply(set)
# test_df2 = df2.groupby(by="int_column")["test"].apply(set)
#
# # 计算Jaccard指数
# jaccard_scores = {}
# for int_val in test_df1.index:
#
#     set1 = test_df1.get(int_val, set())
#     set2 = test_df2.get(int_val, set())
#
#     jaccard_scores[int_val] = jaccard_index(set1, set2)
#
# # 将Jaccard指数添加到df1
# new_df=df1.groupby(by="int_column")["test"].count().reset_index()
# new_df['jaccard_index'] = new_df['int_column'].map(jaccard_scores)
#
# print(new_df.reset_index)
#
# # %%
# from utils import common
# import json
# import Filter
# import pandas as pd
#
# with open("config.json") as p:
#     config = json.load(p)
#
# stringdb = config["stringdb_filepath"]
#
# link = common.read_file(stringdb["link"])
# info = common.read_file((stringdb["info"]))
# info_map = info.set_index("#string_protein_id")["preferred_name"].to_dict()
#
# filter_link = link[link["combined_score"] >= 900]
# filter_link.loc[:, "protein1"] = filter_link["protein1"].map(info_map)
# filter_link.loc[:, "protein2"] = filter_link["protein2"].map(info_map)
#
#%%
import pandas as pd
import networkx as nx

edges = pd.DataFrame(
    {
        "source": [0, 1, 2],
        "target": [2, 2, 3],
        "weight": [3, 4, 15],
        "color": ["red", "blue", "yellow"],
        "width":[1, 4, 5]

    }
)
G = nx.from_pandas_edgelist(edges, edge_attr=True)
pos = nx.spring_layout(G)
tissue_nodes=edges["source"]
target_nodes=edges["target"]
node_degree = dict(G.degree())
draw_node_labels={n:n for n in target_nodes if node_degree[n] > 1}

nx.draw_networkx_nodes(G, pos=pos, nodelist=tissue_nodes,node_size=30, alpha=0.7,node_color="black", node_shape="v")
nx.draw_networkx_nodes(G, pos=pos, nodelist=target_nodes,node_size=10, alpha=0.7,node_color="yellow", node_shape="o")

nx.draw_networkx_edges(G, pos, edge_color=edges["color"], width=edges["width"])
node_labels={node:node for node in G.nodes}

nx.draw_networkx_labels(G, pos, labels=draw_node_labels)
edge_labels=nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
#
# #%%
# import pandas as pd
#
#
# # 创建示例数据
# data1 = {
#     'int_column': [1, 4, 4, 4, 2, 2, 3, 3, 4, 4, 4],
#     'test': ["a", "D", "d", "D", "B", "B", "C", "C", "C", "D", "d"]
# }
# data2={'int_column': [2, 4, 4, 4, 2, 2, 3, 5, 3, 4],
#        'test': ["a", "D", "d", "D", "B", "B", "C", "C", "C", "D"]
#        }
# df1 = pd.DataFrame(data1)
# df2 = pd.DataFrame(data2)
# df1 = df1.set_index(["int_column","test"])
# df2 = df2.set_index(["int_column","test"])
# print(pd.Series(df1.index.isin(df2.index)))
# new_df = df1.copy().reset_index()
# a=pd.Series(df1.index.isin(df2.index))
# new_df["is_same"]= a
#%%
from utils import common
import pandas as pd
import drawNetwork
data1 = {
        'int_column': [1, 4, 4],
        'test': ["g1t4", "g2t5c7", "s2g1t10c10"]
    }
df1 = pd.DataFrame(data1)
df1=df1["test"].apply(drawNetwork.map_id)
print(df1)