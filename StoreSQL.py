import json
import os
import re

import MySQL

with open("./config.json") as j:
    config = json.load(j)

raw_data_folder = config["rawDataFolder"]

for _, folder, file in os.walk(raw_data_folder):
    folder_list = [f for f in folder if not f.startswith('.')]
    break
tissue_list = [t.split(".")[0] for t in folder_list]
tissue_list = set(tissue_list)
sql = MySQL.sql()


# create tables
# with sql as sql:
#     for t in tissue_list:
#         t = re.sub(r"[.+-]", "_", t)
#         cell = t+"_cell"
#         sql.createTable(tableName=t)
#         sql.createTable(tableName=cell)

# store donars content
# donar_file = config["donar_file"]
# with sql as sql:
#     sql.loadTable(donar_file, table_name="donars")
# %%
def read_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
        lines = [l.strip("\n").replace("\t", " ") for l in lines]
    return lines


tissue_test = "/Users/liuyuting/WorkSpace/ScRNA_project/mysql/17_human_TabulaSapiens/Bladder"
file_list = os.listdir(tissue_test)

# %%
# store tissue content:
"""没写完，只存了bladder"""
with sql as sql:
    for file in file_list:
        if not file.startswith("."):
            file_path = os.path.join(tissue_test, file)
            gene_exp = read_file(file_path)
            # pass through blank file
            if len(gene_exp) == 1 and gene_exp[0] == "":
                continue

            for record in gene_exp:
                r = record.split(" ")
                tissue = "Bladder"
                gene_name = file.replace(".txt", "")
                cell_id = r[0]
                gene_expression = r[1]
                age_group = r[2]
                donar = re.search(r"TSP[0-9]+", cell_id).group()
                # print(gene_name, cell_id, gene_expression, age_group, donar)
                sql.store(tableName=tissue,
                          tissue=tissue,
                          cellid=cell_id,
                          donar=donar,
                          agetype=age_group,
                          genename=gene_name,
                          expression=gene_expression
                          )

print("finish")