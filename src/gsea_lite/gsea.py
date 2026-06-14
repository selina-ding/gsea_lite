# gsea.py

import pandas as pd
import numpy as np
from pathlib import Path

def generate_ranked_list(deg_data):
    """
    Generate a ranked gene list for GSEA.

    Parameters
    ----------
    deg_data : pandas.DataFrame
        Differential expression result table.
        Required columns:
            - Gene.symbol
            - logFC

    Returns
    -------
    pandas.DataFrame
        Gene list sorted by logFC (descending) without NaNs and duplicates.
    """

    # 检查输入表中是否包含必要列
    required_columns = ["Gene.symbol", "logFC"]
    for col in required_columns:
        if col not in deg_data.columns:
            raise ValueError(f"Missing required column: {col}")

    # 剔除没有基因名字（NaN）的非特异性行
    cleaned_data = deg_data.dropna(subset=["Gene.symbol"]).copy()

    # 3. 处理重复基因：按 logFC 的绝对值（变化幅度）从大到小排序，保留每个基因表达变化最显著的那一个探针
    cleaned_data["abs_logFC"] = cleaned_data["logFC"].abs()
    cleaned_data = cleaned_data.sort_values(by="abs_logFC", ascending=False)
    cleaned_data = cleaned_data.drop_duplicates(subset=["Gene.symbol"], keep="first")

    # 按 log2FoldChange 从大到小排，上调最显著的基因位于顶部，下调最显著的基因位于底部
    ranked_gene_list = (
        cleaned_data.sort_values(by="logFC", ascending=False)# 按照 log2FoldChange（Log2倍数变化）这一列进行降序排列（ascending=False）
        [["Gene.symbol", "logFC"]]# 只保留Gene_Symbol和log2FoldChange两列
        .reset_index(drop=True)# 按照log2FoldChange排序后原表格行的索引被打乱，drop=True丢弃原来的旧索引，重新从 0开始顺序编号。
    )

    return ranked_gene_list


def calculate_es(ranked_gene_list, pathway_genes):
    """
    Calculate simplified GSEA Enrichment Score (ES).

    Parameters
    ----------
    ranked_gene_list : pandas.DataFrame
        Output of generate_ranked_list()

    pathway_genes : set
        Genes belonging to a biological pathway

    Returns
    -------
    float
        Enrichment Score (ES)
    """

    # 保证 pathway_genes 为集合， 集合在底层基于哈希表实现。
    # Python查找集合时直接用哈希算法，不需要遍历，时间复杂度是O(1)
    pathway_genes = set(pathway_genes)

    # Ranked List总基因数
    total_genes = len(ranked_gene_list)

    # 实际出现在排序列表中的通路基因数
    pathway_hits = len(
        pathway_genes.intersection(
            ranked_gene_list["Gene.symbol"]
        )
    )

    # 如果没有任何通路基因出现在列表中，则无法进行富集分析
    if pathway_hits == 0:
        return 0.0

    # 极端情况保护
    if pathway_hits == total_genes:
        return 0.0

    # 命中通路基因时增加的步长
    hit_score = 1.0 / pathway_hits

    # 未命中通路基因时减少的步长
    miss_score = 1.0 / (total_genes - pathway_hits)

    # 当前运行富集分数
    current_es = 0.0

    # 记录最大偏离值
    max_es = 0.0

    # 从排序列表顶部开始遍历 (修复 Bug 1: 修改为 "Gene.symbol")
    for gene in ranked_gene_list["Gene.symbol"]:

        # 命中通路基因
        if gene in pathway_genes:
            current_es += hit_score
        # 非通路基因
        else:
            current_es -= miss_score

        # 记录绝对值最大的偏离位置（可能是正富集也可能是负富集）
        if abs(current_es) > abs(max_es):
            max_es = current_es

    return max_es


def run_gsea(deg_data, gene_sets):
    """
    Run simplified GSEA for all pathways.

    Parameters
    ----------
    deg_data : pandas.DataFrame
        Differential expression table

    gene_sets : dict
        Example:
        {
            "Apoptosis": {"TP53", "BAX", "CASP3"},
            "Cell_Cycle": {"CDK1", "CCNB1"}
        }

    Returns
    -------
    pandas.DataFrame
        Pathway enrichment results
    """
    # Step 1: 生成标准的基因排序列表
    ranked_gene_list = generate_ranked_list(deg_data)

    results = []

    # Step 2: 对每个通路计算ES
    for pathway_name, pathway_genes in gene_sets.items():
        es = calculate_es(ranked_gene_list, pathway_genes)

        results.append({
            "Pathway": pathway_name,
            "ES": es,
            "Gene_Count": len(pathway_genes)
        })

    # Step 3: 按ES绝对值降序排序，最显著富集的通路排在最前
    results_df = pd.DataFrame(results)
    results_df = (
        results_df
        .sort_values(by="ES", key=abs, ascending=False)
        .reset_index(drop=True)
    )

    return results_df