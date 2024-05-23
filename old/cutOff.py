import math

import MySQL
import pandas as pd
from utils import common
from utils.log import logger
from database import db
import os
import json


def save_json(data, path):
    with open(path, "w") as j:
        json.dump(data, j)


class QueryDB:
    def __init__(self, database):
        self.db = db.databaseHelper(database)
        self.db.init()
        self.log = logger()

    def set_filter(self, p_threshold, corr_threshold):
        """init p threshold, correlation threshold and log10p records above p threshold in the database"""
        self.p_threshold = p_threshold
        self.corr_threshold = corr_threshold
        self.log10p_dict_list = self.get_log10p_dict_in_all_studies()

    def get_avaliable_studies(self) -> list:
        """extract all abbr_id tables in the database"""
        _, tasks = self.db.get_available_tasks()
        tasks = [i["Table"] for i in tasks]
        self.log.debug(f"available_studies--->{len(tasks)}:{tasks}")
        return tasks

    def get_log10p_dict_in_all_studies(self) -> list[dict]:
        """extract list of log10p records above p_shreshold in the database"""

        flag, self.log10p_dict_list = self.db.get_log10p_records_above_threshold(self.p_threshold)
        self.log.debug(f"log10p_dict init 🎉🎉🎉～～～～～{flag}")
        return self.log10p_dict_list

    def get_corr_dict_in_one_study(self, abbr_id) -> list[dict]:
        """get dict of gene pair info if over corr_threshold and p_threshold in single study"""
        new_corr_dict_list = []
        flag, corr_dict_list = self.db.get_correlation_records_above_threshold(abbr_id, self.corr_threshold)
        temp_gene_dict = {record["gene"]: record["log10p"] for record in self.log10p_dict_list if
                          record["abbr_id"] == abbr_id}  # 得到在当前study中的gene到log10p的映射mapping
        if flag:
            self.log.debug(f"--- get corr_dict_in_one_study {abbr_id} ---> flag {flag}")

            for record in corr_dict_list:
                # print("record waiting search", record)
                gene1_log10p = temp_gene_dict.get(record["gene1"], None)
                gene2_log10p = temp_gene_dict.get(record["gene2"], None)
                if gene1_log10p and gene2_log10p:
                    record["log10p1"] = gene1_log10p
                    record["log10p2"] = gene2_log10p
                    record["abbr_id"] = abbr_id
                    # print("record searched", record)
                    new_corr_dict_list.append(record)
                else:
                    pass
        return new_corr_dict_list

    def get_corr_dict_in_all_study(self, path, save_json_file=True) -> list[dict]:
        """get dict of gene pair info if over corr_threshold and p_threshold in all studies, save json file if True"""

        tasks = self.get_avaliable_studies()
        corr_dict_in_all_study = []
        count = 0
        for study in tasks:
            self.log.debug(f"get corr_dict {study} ---> START")
            corr_dict_in_all_study += self.get_corr_dict_in_one_study(study)
            if save_json_file:
                save_json(corr_dict_in_all_study, path + "screen.json")
                count += 1
                self.log.debug(f"save corr_dict in json ---> count:{count}/{len(tasks)} study:{study}")
        self.log.debug(f"finish corr_dict SHOW top2 ---> {corr_dict_in_all_study[0:2]}")
        return corr_dict_in_all_study


# %%
with open("./config.json") as c:
    config = json.load(c)

figSavePath = config["figSavePath"]
db_list = config["db_list"]
corr_cutoffs = config["corr_cutoffs"]
log10p_abs_cutoffs = config["log10p_abs_cutoffs"]
save_title_df = "screen.csv"

for database, save_path in zip(db_list, figSavePath.values()):
    logger().log.debug(f"USE {database}, save path--->{save_path}")
    query = QueryDB(database)
    for p in log10p_abs_cutoffs:
        for c in corr_cutoffs:
            folder = f"log10p_{p}_corr_{c}/"
            common.create_folder(save_path + folder)
            query.set_filter(p, c)
            corr_dict = query.get_corr_dict_in_all_study(save_path + folder, save_json_file=True)

            # df = pd.DataFrame(corr_dict)
            # df.to_csv(save_path + folder + save_title_df)
