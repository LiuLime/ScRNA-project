#!/bin/bash

# 对当前目录下的每个文件进行操作
for file in *_merge.csv; do
    # 使用'.'作为分隔符提取文件名的前缀
    prefix=$(echo $file | cut -d'.' -f1)

    # 如果子文件夹不存在，则创建它
    [ ! -d "$prefix" ] && mkdir "$prefix"

    # 将文件移动到相应的子文件夹中
    mv "$file" "$prefix/"
done