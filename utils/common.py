import os

import pandas as pd
import csv


def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)


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
