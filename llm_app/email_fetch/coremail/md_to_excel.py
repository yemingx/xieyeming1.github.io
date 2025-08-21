import re
import pandas as pd

def md_to_excel(md_path, excel_path):
    with open(md_path, 'r') as f:
        lines = f.readlines()
    
    # 提取表格行（跳过分隔行）
    table_lines = [line.strip() for line in lines if '|' in line and not re.match(r'^\|?-+\|', line)]
    
    # 构建二维列表
    data = []
    for line in table_lines:
        cols = [col.strip() for col in line.split('|') if col.strip()]
        data.append(cols)
    
    # 转为 DataFrame 并导出
    df = pd.DataFrame(data[1:], columns=data[0])
    df.to_excel(excel_path, index=False)

# 调用示例
md_to_excel("tmp/wangjn_bgi.com_250710.md", "output.xlsx")


def md_to_csv(md_path, tsv_path):
    with open(md_path, 'r') as f:
        lines = f.readlines()
    
    # 提取表格行（跳过分隔行）
    table_lines = [line.strip() for line in lines if '|' in line and not re.match(r'^\|?-+\|', line)]
    
    # 构建二维列表
    data = []
    for line in table_lines:
        cols = [col.strip() for col in line.split('|') if col.strip()]
        data.append(cols)
    
    # 转为 DataFrame 并导出
    df = pd.DataFrame(data[1:], columns=data[0])
    df.to_csv(tsv_path, index=False)

# 调用示例
md_to_csv("tmp/wangjn_bgi.com_250710.md", "output.csv")

import csv

def md_to_tsv(md_path, tsv_path):
    with open(md_path, 'r') as f:
        lines = [line.strip() for line in f if '|' in line][2:]  # 跳过表头和分隔行
    
    with open(tsv_path, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        for line in lines:
            # 清洗数据：去除首尾管道符，分割列
            cols = [col.strip() for col in line.split('|') if col.strip()]
            writer.writerow(cols)

# 调用示例
md_to_tsv("tmp/wangjn_bgi.com_250710.md", "output.tsv")