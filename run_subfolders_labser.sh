#!/bin/bash

PYTHON_SCRIPT_PATH="/home/lyt/ScRNA_project/model/lr_batch.py"
counter=0

for dir in */; do
    abs_dir_path=$(realpath "$dir")
    save_path="/home/lyt/ScRNA_project/reports/"
    file_end="merge.csv"

    # 启动Python进程并将输出重定向到日志文件
    nohup python3 "$PYTHON_SCRIPT_PATH" "$abs_dir_path" "$save_path" "$file_end" > "output_$!.log" 2>&1 &

    counter=$((counter+1))
    if [[ $counter -eq 10 ]]; then
        wait
        counter=0
    fi
done


