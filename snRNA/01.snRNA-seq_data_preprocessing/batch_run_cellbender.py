import os

path_list = [
    "/share/org/BGI/bgi_baiyq/lai_yiwei/xiangyuchen/rawmatrix_data1/W202504160000371",
    "/share/org/BGI/bgi_baiyq/lai_yiwei/xiangyuchen/rawmatrix_data1/W202504160000361",
    "/share/org/BGI/bgi_baiyq/lai_yiwei/xiangyuchen/rawmatrix_data1/W202504170002525",
    "/share/org/BGI/bgi_baiyq/lai_yiwei/xiangyuchen/rawmatrix_data1/W202504160000313",
    "/share/org/BGI/bgi_baiyq/lai_yiwei/xiangyuchen/rawmatrix_data1/W202504160000303"]

# 遍历 path_list 中的每个路径
for target_path in path_list:
    # 查找当前路径下的 h5ad 文件
    h5ad_file = None
    for file in os.listdir(target_path):
        if file.endswith(".h5ad"):
            h5ad_file = os.path.join(target_path, file)
            break  # 假设每个子路径中只有一个 h5ad 文件，找到后退出循环

    if h5ad_file:
        output_file = os.path.join(target_path, "cellbender_output.h5")
        os.chdir(target_path)
        cellbender_command = (
            f"nohup cellbender remove-background "
            f"--input {h5ad_file} "
            f"--output {output_file} "
            f"--projected-ambient-count-threshold 5 "
            f"--checkpoint-mins 120 "
            f"--learning-rate 0.0001 "
            f"--training-fraction 0.9 "
            f"--low-count-threshold 20 "
            f"--epochs 150 "
            f"--cuda &"
        )
        print(cellbender_command)
        os.system(cellbender_command)
        os.chdir(target_path)
