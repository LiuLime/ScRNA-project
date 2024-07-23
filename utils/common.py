import os
import re
import pandas as pd
import csv
from utils import log

logger = log.logger()


def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)


def read_file(path) -> pd.DataFrame:
    with open(path, "r") as p:
        delimiter = detect_delimiter(path)
        data = pd.read_csv(p, header=0, sep=delimiter)
    return data


def save_csv(df: pd.DataFrame, path):
    df.to_csv(path, index=False, sep=",")
    logger.debug(f"save dataframe {path} 🍺")


def replace_name(path, oldchar=":", newchar="|") -> None:
    """ replace file names in folder
    param:
    path: folder path
    oldchar: old character
    newchar: new replace character
    """
    # 指定原始文件夹路径
    folder_path = path

    # 获取文件夹中的文件名列表
    file_list = os.listdir(folder_path)

    # 循环遍历文件名列表
    for filename in file_list:
        # 构建原始文件的完整路径
        old_filepath = os.path.join(folder_path, filename)

        # 构建新的文件名
        new_filename = filename.replace(oldchar, newchar)

        # 构建新的文件路径
        new_filepath = os.path.join(folder_path, new_filename)

        # 修改文件名
        os.rename(old_filepath, new_filepath)

        # 打印修改后的文件路径
        print(f"文件名已修改：{new_filepath}")


def detect_delimiter(file_path):
    """检测文件分隔符"""
    with open(file_path, 'r') as file:
        sample = file.readline()  # 读取文件的第一行来测试
        sniffer = csv.Sniffer()
        sniffer.preferred = [';', ',', '\t', ' ']
        dialect = sniffer.sniff(sample)
        return dialect.delimiter


class converter:
    def __init__(self):
        study_id = read_file("./01datasource/study_id.csv")
        group_id = read_file("./01datasource/group_id.csv")
        tissue_id = read_file("./01datasource/tissue_id.csv")
        cell_id = read_file("./01datasource/cell_id.csv")

        self.study_id_dict = study_id.set_index("study_id")["study"].to_dict()
        self.group_id_dict = group_id.set_index("group_id")["health_group"].to_dict()
        self.tissue_id_dict = tissue_id.set_index("tissue_id")["tissue"].to_dict()
        self.cell_id_dict = cell_id.set_index("cell_id")["cell"].to_dict()

    @staticmethod
    def parse_string(input_string):
        pattern = r'^(s\d+)?(g\d+)(t\d+)(c\d+)?$'
        match = re.match(pattern, input_string)

        if match:
            s, g, t, c = match.groups()  # 这将返回一个包含结果的元组，如 ('s1', 'g2', 't2', 'c3') 或 ('s1', 'g2', 't2', None)
            return s, g, t, c
        else:
            return "No match found."

    def convert_abbrid_to_fullname(self, abbr_id: str):
        """convert abbr id to full id"""

        study, group, tissue, cell = converter.parse_string(abbr_id)
        study_full_id = self.study_id_dict.get(study)
        group_full_id = self.group_id_dict.get(group)
        tissue_full_id = self.tissue_id_dict.get(tissue)
        cell_full_id = self.cell_id_dict.get(cell)
        full_id = "|".join([i for i in [study_full_id, group_full_id, tissue_full_id, cell_full_id] if i is not None])
        return full_id


class process_title:
    def __init__(self, level):
        self.level = level
        self.title = None
        self.study = None
        self.group = None
        self.tissue = None
        self.cellType = None

    def extract_info(self, title, titleSegment="|"):
        """conduct single file title, return study, group, tissue (and cellType) info from title
        :param:
        * title
        * titleSegment, default = '|'

        :return: self.study, self.group, self.tissue, (self.cellType)
        """
        title = title.replace("allGenePairs_binMeanPearsonCor_mtx_pLT0.2_1ageGrpGE50genes.", "")
        self.title = title.replace(".txt", "")
        title_segment = title.split(titleSegment)
        levels = {"tissue": -2, "cell": -3}
        cut_off = levels.get(self.level, None)
        if cut_off is None:
            raise TypeError(f"missing argument 'level':'tissue' or 'cell'")
        self.study = ":".join(title_segment[:cut_off])
        details = title_segment[cut_off:]

        if self.level == "tissue":
            self.group, self.tissue = details
            return self.study, self.group, self.tissue
        elif self.level == "cell":
            self.group, self.tissue, self.cellType = details
            return self.study, self.group, self.tissue, self.cellType
