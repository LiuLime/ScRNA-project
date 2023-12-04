

import os

# 指定原始文件夹路径
folder_path = '/Users/liuyuting/ScRNA_project/test/'

# 获取文件夹中的文件名列表
file_list = os.listdir(folder_path)

# 循环遍历文件名列表
for filename in file_list:
    # 构建原始文件的完整路径
    old_filepath = os.path.join(folder_path, filename)

    # 构建新的文件名
    new_filename = filename.replace(':', '|')

    # 构建新的文件路径
    new_filepath = os.path.join(folder_path, new_filename)

    # 修改文件名
    os.rename(old_filepath, new_filepath)

    # 打印修改后的文件路径
    print(f"文件名已修改：{new_filepath}")

