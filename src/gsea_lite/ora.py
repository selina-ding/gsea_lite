from pathlib import Path
import pandas as pd
from scipy.stats import hypergeom

def extract_significant_genes(df, p_thresh=0.05, lfc_thresh=1.0):
    """
    根据给定的p-value和|log2FC|阈值,从差异表达矩阵中筛选出显著差异基因DEGs
    """

    condition = (df['P.Value'] < p_thresh) & (df['logFC'].abs() > lfc_thresh)
    
    genes_filtered = df.loc[condition, 'Gene.symbol'] 
    
    sig_genes = set(genes_filtered.dropna().unique()) # 去除缺失值和重复基因，转换为集合,便于后续的交集运算
    
    return sig_genes

def do_ora(df, gene_sets, p_thresh=0.05, lfc_thresh=1.0):
    """
    ORA分析,计算每个通路的 p-value
    """
    
    # 1. 定义全局背景基因库
    universe_genes = set(df['Gene.symbol'].dropna().unique())
    M = len(universe_genes) 
    
    # 2. 提取显著差异基因集
    deg_genes = extract_significant_genes(df, p_thresh, lfc_thresh)
    deg_genes = deg_genes.intersection(universe_genes) # 确保差异基因一定在背景库中
    N = len(deg_genes) 

    print(f"--- 背景基因库说明 ---")
    print(f"提取到全局背景基因 (Universe) 数量: {M}")
    print(f"提取到显著差异基因 (DEGs) 数量: {N}")
    print(f"原理解释: 仅使用在当前测序平台/微阵列中实际检测到的 {M} 个基因作为超几何检验的背景，而不是整个人类基因组，能有效避免因技术偏倚导致的假阳性。")
    print(f"----------------------")

    if N == 0:
        print("Warning: no significant genes found. ")
        return pd.DataFrame(columns=[
            'pathway', 'P.Value', 'Count', 'pathway_size', 'GeneRatio', 'hits_genes'
        ])

    # 3. 遍历 GMT 字典中的每一个通路，进行超几何检验
    results = []
    
    for pathway_name, pathway_genes in gene_sets.items():
        pathway_set = set(pathway_genes)
        
        pathway_in_universe = pathway_set.intersection(universe_genes)
        n = len(pathway_in_universe)
        
        if n == 0:
            continue
            
        # 计算交集 k：既属于该通路，又是显著差异的基因
        hit_genes = deg_genes.intersection(pathway_in_universe)
        k = len(hit_genes)
        raw_p = hypergeom.sf(k - 1, M, n, N)
        
        results.append({
            'pathway': pathway_name,
            'P.Value': raw_p,
            'Count': k,                                  
            'pathway_size': n,
            'GeneRatio': k / N if N > 0 else 0,          
            'hits_genes': ",".join(list(hit_genes))
        })
        
    ora_result_df = pd.DataFrame(results)
    
    # 按照原始 p-value 从小到大排列（最显著的排在最前）
    if not ora_result_df.empty:
        ora_result_df = ora_result_df.sort_values(by='P.Value').reset_index(drop=True)
        
    return ora_result_df