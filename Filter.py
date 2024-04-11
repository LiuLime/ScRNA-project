"""
根据相关性和p-value筛出dataframe切片，并进行可视化

"""

import os
import json
import pandas as pd
import utils as us
import DrawTool
from MySQL import sql
from utils import common
from utils.log import logger


class dataFilter:
    def __init__(self, db, corr_threshold, p_threshold):
        self.db = db
        self.corr_threshold = corr_threshold
        self.p_threshold = p_threshold
        self.s = sql(self.db)

    def filter_data_by_criteria(self):
        """Filter dataframe by correalation value and pvalue, return dataframe slice"""
        with self.s:
            tasks = self.s.get_avalible_studies()
            # self.s.join_corr_log10p(tasks)




