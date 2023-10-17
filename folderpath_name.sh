for dir in */; do
    # 获取当前子文件夹的绝对路径
    abs_dir_path=$(realpath "$dir")
    echo "${abs_dir_path}"
done
