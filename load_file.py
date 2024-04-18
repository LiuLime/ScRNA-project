"""transform pvalue file and correlation gene pair file as long table and load to mysql database
2024/4/12
"""

import math

from MySQL import sql
import pandas as pd
from utils import common
from utils.log import logger
import os
import json

with open("./config.json") as c:
    config = json.load(c)


def read_file(path):
    with open(path, "r") as p:
        delimiter = common.detect_delimiter(path)
        data = pd.read_csv(p, delimiter=delimiter, header=0)

    return data



class preprocessor:
    def __init__(self):
        """load source file in mysql database"""
        self.sql = None
        self.log = logger()

    def wide_to_long_p(self, folder, file):
        p_file = read_file(os.path.join(folder, file))
        new_p = pd.melt(p_file, id_vars=["gene"], value_vars=p_file.columns[1:])
        new_p = new_p.set_axis(["gene", "full_id", "log10p"], axis=1)

        unique_full_ids = new_p["full_id"].unique().tolist()
        id_map = preprocessor.fullID_to_abbrID(unique_full_ids)  # convert full id to abbrevation id

        new_p["abbr_id"] = new_p["full_id"].map(id_map)
        new_p = new_p.drop("full_id", axis=1)

        return new_p

    def CreateAbbrId(self, row: pd.Series):
        """return abbrevation id for each full study name"""
        # s = MySQL.sql("scrna_base")
        # with s:
        study_id = self.sql.search("study_id", "study_id", f"study = '{row['study']}'")[0][0]
        group_id = self.sql.search("group_id", "group_id", f"health_group = '{row['group']}'")[0][0]
        tissue_id = self.sql.search("tissue_id", "tissue_id", f"tissue = '{row['tissue']}'")[0][0]
        cell_id = ""
        try:
            cell_id = self.sql.search("cell_id", "cell_id", f"cell = '{row['cellType']}'")[0][0]
        except KeyError:
            pass
        abbr_id = "".join([study_id, group_id, tissue_id, cell_id])
        return abbr_id

    def fill_abbr_id(self, file_path):
        """ create abbr_id for each study full title according to mysql database `study_id`,`group_id`,`tissue_id`,
        `cell_id`"""

        data = pd.read_csv(file_path, delimiter=common.detect_delimiter(file_path), header=0)
        # print(data.columns)
        data["full_id"] = data.apply(lambda x: ":".join(x.astype(str)), axis=1)

        self.sql = sql("scrna_base")
        with self.sql:
            data["abbr_id"] = data.apply(self.CreateAbbrId, axis=1)
        return data

    @staticmethod
    def replace_file(file):
        return file.replace("allGenePairs_binMeanPearsonCor_mtx_pLT0.2_1ageGrpGE50genes.", "").replace(".txt",
                                                                                                       "").replace("|",
                                                                                                                   ":")

    @staticmethod
    def fullID_to_abbrID(unique_full_ids) -> dict:
        """return full_id to abbr_id map"""
        sql_base = sql("scrna_base")

        with sql_base:
            query = f"SELECT full_id, abbr_id FROM abbr_id WHERE full_id IN ({','.join(['%s'] * len(unique_full_ids))})"
            results = sql_base.execute_query(query, unique_full_ids)
        id_map = {full_id: abbr_id for full_id, abbr_id in results}
        return id_map

    def wide_to_long_corr_batch(self, folder, db):
        """convert correlation matrix into long matrix.
        Batch process files in folder.
        Save dataframe as table in mysql database
        :param: folder:path-like
        :param: db: mysql database name for saving tables
        """
        # create id map dict for title(full_id) to abbr_id
        file_list = os.listdir(folder)
        unique_full_ids = [preprocessor.replace_file(file) for file in file_list]
        # print("unique full ids-----",unique_full_ids)
        id_map = preprocessor.fullID_to_abbrID(unique_full_ids)
        print("studies in this project-----", len(id_map))
        # create table
        sql_corr = sql(db)

        # fill data into tables
        for file, full_id in zip(file_list, unique_full_ids):
            message = f"{id_map[full_id]} {full_id}"
            self.log.debug(f"process corr df--->{message}")
            path = os.path.join(folder, file)
            with open(path, "r") as f:
                data = pd.read_csv(f, delimiter=common.detect_delimiter(path), header=0).set_index("Unnamed: 0")
            df_long = data.stack().reset_index()  # Not include NA value pairs
            df_long.columns = ['gene1', 'gene2', 'cor_pearson.binMean']

            df_list = df_long.to_dict(orient="list")
            df_list_record = list(zip(*df_list.values()))  # mysql table record content, each row as tuple
            columns = f"{tuple(df_list.keys())}".replace('\'', "`")  # mysql table column name
            # print("Input column names----", columns)
            trunksize = 5000
            truncktime = math.ceil(len(df_long) / trunksize)
            with sql_corr:
                query1 = f"CREATE table if not exists {id_map[full_id]}(`gene1` varchar(30), `gene2` varchar(30), `cor_pearson.binMean` DOUBLE,PRIMARY KEY (`gene1`, `gene2`))"
                sql_corr.execute_query(query1, return_object=False)
                for times in range(truncktime):
                    batch = df_list_record[times * trunksize: (times + 1) * trunksize]
                    value = ','.join(str(i) for i in batch)
                    query2 = f"INSERT IGNORE INTO {id_map[full_id]}{columns} VALUES {value}"
                    sql_corr.execute_query(query2, return_object=False)


save_path = "/Users/liuyuting/WorkSpace/ScRNA_project/mysql/"
# p-value file path
ageGrp_byMedian_tissue = config["ageGrp_byMedian_tissue"]
ageGrp_by4060_tissue = config["ageGrp_by40-60_tissue"]
tissuePvalueFile = config["tissuePvalueFile"]

ageGrp_byMedian_cell = config["ageGrp_byMedian_cell"]
ageGrp_by4060_cell = config["ageGrp_by40-60_cell"]
cellPvalueFile = config["cellPvalueFile"]
# study list folder
tissueList = config["tissueList"]
cellTypeList = config["cellTypeList"]
# correlation file
CorrFolder = config["CorrFolder"]

# p folder
folder_list = [ageGrp_byMedian_tissue, ageGrp_byMedian_cell, ageGrp_by4060_tissue, ageGrp_by4060_cell]
pfile_list = [tissuePvalueFile, cellPvalueFile, tissuePvalueFile, cellPvalueFile]
study_list = [tissueList, cellTypeList, tissueList, cellTypeList]

id_path_list = [ageGrp_byMedian_tissue + tissueList, ageGrp_byMedian_cell + cellTypeList]  # median age study id file
db_list = ["scrna_mt", "scrna_mc", "scrna_4060t", "scrna_4060c"]
corr_path_list = [group + CorrFolder for group in folder_list]  # correlation folder in 4 groups
folder_db = {"scrna_mt": ageGrp_byMedian_tissue + CorrFolder,
             "scrna_mc": ageGrp_byMedian_cell + CorrFolder,
             "scrna_4060t": ageGrp_by4060_tissue + CorrFolder,
             "scrna_4060c": ageGrp_by4060_cell + CorrFolder}

""" 保存所有study对应id"""
# for folder, file in list(zip(folder_list[2:], study_list[2:])):
#     p = preprocessor()
#     id_df = p.fill_abbr_id(os.path.join(folder, file)).rename(columns={"group": "health_group"})
#     df_list = id_df.to_dict(orient="list")
#     df_list_record = list(zip(*df_list.values()))
#     value = ','.join(str(i) for i in df_list_record)
#     s = sql("scrna_base")
#     with s:
#         query = f"INSERT IGNORE INTO `abbr_id` (`study`,`health_group`,`tissue`,`cellType`,`full_id`,`abbr_id`) VALUES {value}"
#         s.execute_query(query, return_object=False)

# id_df.to_csv(folder + "study_abbr_id_list.csv", sep=",", index=False)

"""保存p-value的long table"""
# for folder, file, db in list(zip(folder_list, pfile_list, db_list)):
#     print(file, db)
#     p = preprocessor()
#     new_p = p.wide_to_long_p(folder, file)
#
#     s = sql(db)
#     with s:
#         s.create_log10p_table()
#         s.df_to_sql(new_p, 'log10p')

"""保存correlation的long table"""
# for db, corr_path in folder_db.items():
#     print(f"conduct {db} database from path {corr_path}")
#     p = preprocessor()
#     p.wide_to_long_corr_batch(corr_path, db)

"""Create Index"""
# for db_task in db_list:
#     s = sql(db_task)
#     with s:
#         tasks = s.get_avalible_studies()
#         s.create_index(tasks)

"""join p file and corr file"""
for db_task in db_list:
    s = sql(db_task)
    with s:
        tasks = s.get_avalible_studies()
        s.join_corr_log10p(tasks)
# %%
