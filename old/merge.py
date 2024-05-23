import os, sys
import pandas as pd


def merge_data_in_folder(filepath):
    # 读取文件夹内所有txt文件
    filenames = [f for f in os.listdir(filepath) if f.endswith(".txt")]

    dfs = []
    label_dict = {}

    for filename in filenames:
        file_path = os.path.join(filepath, filename)
        if os.path.getsize(file_path) >= 2:
            gene_name = filename.replace(".txt", "")

            temp_df = pd.read_csv(file_path, sep="\t", header=None, names=['cell_id', gene_name, 'label'])  # names自定义列
            label_dict.update(temp_df.set_index('cell_id')['label'].to_dict())  # 更新label字典

            temp_df = temp_df.set_index('cell_id').drop(columns='label')  # 将cell_id设为索引并删除label列
            dfs.append(temp_df)
        else:
            pass

    merged_df = pd.concat(dfs, axis=1)  # 根据cell_id索引进行横向拼接

    # 创建label数据帧并与主数据帧合并
    label_df = pd.DataFrame(list(label_dict.items()), columns=['cell_id', 'label']).set_index('cell_id')
    merged_df = pd.concat([merged_df, label_df], axis=1).reset_index().rename(columns={'index': 'cell_id'})
    
    print('merge success')

    return merged_df

def main(argv):
    FOLDER_PATH = argv[0]
    SAVE_PATH = argv[1]
    cell_types = [foldername for foldername in os.listdir(FOLDER_PATH) if not foldername.startswith('.')]
    for cell_type in cell_types:
        print(f"{cell_type} start process")
        filepath = os.path.join(FOLDER_PATH, cell_type)
        merged_df = merge_data_in_folder(filepath)
        merged_df.to_csv(f'{SAVE_PATH}{cell_type}_merge.csv', index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python ./main.py {DATA_FOLDER_PATH} {SAVE_PATH}')
        print('Example: DATA_FOLDER_PATH = /Users/human_TabulaSapiens/exp/ SAVE_PATH = /Users/human_TabulaSapiens/02results/merge_df/')
        sys.exit(1)
    main(sys.argv[1:])