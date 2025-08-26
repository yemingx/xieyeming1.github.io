#!/bin/bash

# 检查是否为 CentOS 7 系统
if [ -f "/etc/os-release" ]; then
    . /etc/os-release
    if [ "$ID" != "centos" ] || [ "$VERSION_ID" != "7" ]; then
        echo "当前系统不是 CentOS 7 x86，不支持此脚本"
        exit 1
    fi
else
    echo "未能检测到系统信息文件 /etc/os-release"
    exit 1
fi

# 检测系统版本和内核
detect_system() {
    echo
    echo "感谢使用咕噜云 | 阮绘科技一键换源脚本"
    echo
    distributor_id=$(cat /etc/centos-release)
    echo "系统: $distributor_id"
    echo
    jiagou=$(uname -m)
    echo "系统架构: $jiagou"
    echo
    kernel_version=$(uname -r)
    echo "内核版本: $kernel_version"
    echo
}

detect_system

# 软件源数组
declare -a repo_urls=(
    "http://107.149.212.83/CentOS/CentOS-ali.repo|阿里源"
    "http://107.149.212.83/CentOS/CentOS-Base.repo|腾讯源【推荐使用腾讯源比较快】"
    "http://107.149.212.83/CentOS/CentOS-hw.repo|华为源"
    "http://107.149.212.83/CentOS/CentOS-kernel.repo|kernel源"
    "http://107.149.212.83/CentOS/CentOS-gulu.repo|咕噜源[保有备用速度较慢]"
)

# 备份文件夹路径
backup_dir="/etc/yum.repos.d/backup"

# 创建备份文件夹（如果不存在）
if [ ! -d "$backup_dir" ]; then
    sudo mkdir -p "$backup_dir"
    echo "已创建备份文件夹：$backup_dir"
fi

# 显示可用的软件源列表
echo "可用的软件源列表"
echo
for ((i=0; i<${#repo_urls[@]}; i++)); do
    name=$(echo ${repo_urls[$i]} | cut -d '|' -f 2)
    echo "$i. $name"
done

# 提示用户选择一个源
echo
read -p "请选择要使用的软件源的编号: " repo_index

# 验证用户输入的编号是否有效
if [[ $repo_index =~ ^[0-9]+$ ]] && [ $repo_index -ge 0 ] && [ $repo_index -lt ${#repo_urls[@]} ]; then
    # 获取选择的源地址
    repo_url=$(echo ${repo_urls[$repo_index]} | cut -d '|' -f 1)

    # 检测源是否正常
    if ping -c 1 $(echo $repo_url | awk -F/ '{print $3}') &> /dev/null; then
        echo "源 $repo_url 可正常连接。"
    else
        echo "源 $repo_url 无法连接，请检查网络连接。"
        exit 1
    fi

    # 备份当前的源文件
    for file in /etc/yum.repos.d/*.repo; do
        if [ -f "$file" ]; then
            sudo mv "$file" "$backup_dir/"
            echo "已移动备份文件：$file 到 $backup_dir"
        fi
    done

    # 备份所有可能的文件（非.repo格式），排除备份文件夹本身
    for file in /etc/yum.repos.d/*; do
        if [ -f "$file" ] && [[ "$file" != *.repo ]] && [[ "$file" != "$backup_dir"/* ]]; then
            sudo mv "$file" "$backup_dir/"
            echo "已移动非.repo文件：$file 到 $backup_dir"
        fi
    done

    # 判断系统上是否存在 wget 或 curl
    if command -v wget >/dev/null 2>&1; then
        download_cmd="wget -O /etc/yum.repos.d/CentOS-Base.repo $repo_url"
    elif command -v curl >/dev/null 2>&1; then
        download_cmd="curl -o /etc/yum.repos.d/CentOS-Base.repo $repo_url"
    else
        echo
        echo -e "\e[1;31m未安装 wget 和 curl，请安装其中一个后再继续执行\e[0m"
        echo
        exit 1
    fi

    # 执行下载命令
    sudo $download_cmd

    # 清除缓存
    sudo yum clean all

    # 生成缓存
    sudo yum makecache
    echo
   
    echo -e "\e[1;31mYUM源已更换重建缓存完成，不需要更新系统到此结束脚本即可\e[0m"
    # 提示用户是否更新系统
    echo 
    distributor_id=$(cat /etc/centos-release)
   
    
    echo -e "\e[1;31m当前系统版本: $distributor_id  ，选择更新会更新到最新版7.9，不需要更新终止脚本即可\e[0m"
    echo
    read -p "是否要更新到最新系统？不更新输入n或者按ctrl+c即可(y/n): " update_choice
    if [[ "$update_choice" == "y" || "$update_choice" == "Y" ]]; then
    echo
        sudo yum update -y
        echo
        #echo -e "\e[1;32m系统更新完成\e[0m"
        echo
    else
        echo "系统未更新。"
    fi

    # 清理残余文件
    sudo package-cleanup --oldkernels --count=1 -y
    sudo yum autoremove -y
    
    echo
    echo -e "\e[1;31m恭喜已完成配置\e[0m"
    echo
    cat /etc/centos-release
    echo
else
    echo -e "\e[1;31m无效的选择，请重新运行脚本并选择正确的编号\e[0m"
fi
