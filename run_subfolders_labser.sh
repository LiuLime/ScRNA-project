#!/bin/bash

PYTHON_SCRIPT_PATH="/home/lyt/ScRNA_project/model/lr_batch.py"
counter=0

for dir in */; do
    abs_dir_path=$(realpath "$dir")
    save_path="/home/lyt/ScRNA_project/reports/"
    file_end="merge.csv"

    # 启动Python进程
    nohup python3 "$PYTHON_SCRIPT_PATH" "$abs_dir_path" "$save_path" "$file_end" > /dev/null 2>&1 &

    # 获取最近的后台进程ID
    pid=$!
    
    # 使用进程ID命名日志文件
    log_file="output_${pid}.log"
    nohup python3 "$PYTHON_SCRIPT_PATH" "$abs_dir_path" "$save_path" "$file_end" > "$log_file" 2>&1 &

    counter=$((counter+1))
    if [[ $counter -eq 10 ]]; then
        wait
        counter=0
    fi
done

