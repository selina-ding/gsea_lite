
"""
1. 解析GMT文件
2. 读取GEO的TSV文件中的差异表达数据
"""

import os
import pandas as pd

def read_gmt(gmt_file_path):
    """
    解析GMT文件

    返回一个字典：
    key = 通路名称
    value = 该通路含有的基因列表
    """
    gene_sets = {} 
    
    # 对于 file 中的每一行 line 进行处理
    with open(gmt_file_path, 'r', encoding='utf-8') as file:
        for line in file: 
            line = line.strip() 

            if not line: 
                continue 
            
            # 按 \t分割，得到一个列表items 
            items = line.split('\t') 
             
            if len(items) < 2: 
                continue 
            
            # 提取列表中的通路名称、基因
            pathway_name = items[0]
            genes = items[2:] # items[1] 通常是描述信息，忽略 
            
            # 清洗基因列表并存入字典
            valid_genes = []            

            for g in genes:                   
                if g.strip():                
                    valid_genes.append(g.strip()) 

            gene_sets[pathway_name] = valid_genes 
            
    return gene_sets 


# 修改 io.py 中的 read_geo 函数
def read_geo(geo_file_path):
    """
    读取GEO的TSV文件中的关键数据
    返回:经过过滤的差异表达数据dataframe
    """
    # 读取数据
    df = pd.read_csv(geo_file_path, sep='\t', encoding='utf-8') 
    # 验证是否包含核心列
    required_columns = ['Gene.symbol', 'logFC', 'P.Value'] 
    for col in required_columns: 
        if col not in df.columns: 
            raise ValueError(f"error: 缺失列 {col}") 
    # 去除空行,并加上 .copy() 避免后续报错            
    df = df.dropna(subset=required_columns).copy() 
    
    # 清洗逻辑：将包含 '///' 的字符串只保留第一个基因符号
    df['Gene.symbol'] = df['Gene.symbol'].astype(str).apply(lambda x: x.split('///')[0].strip())
    
    return df