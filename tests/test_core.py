import pytest
import pandas as pd
import numpy as np
from gsea_lite.stats import bh_fdr_correction
from gsea_lite.ora import extract_significant_genes, do_ora
from gsea_lite.gsea import generate_ranked_list, calculate_es

# 1. 检验 stats 模块中的 BH-FDR 校正

def test_bh_fdr_correction():
    """测试1: 验证 Benjamini-Hochberg 算法的正确性与单调性"""
    p_values = [0.01, 0.04, 0.03, 0.20]
    fdrs = bh_fdr_correction(p_values)
    
    assert len(fdrs) == 4
    assert np.all(fdrs >= 0) and np.all(fdrs <= 1)
    # 原始 p-value 小的，校正后的 FDR 理论上也应该更小或相等
    assert fdrs[0] <= fdrs[1]

# 2. 检验 ora 模块中的富集检验逻辑
def test_extract_significant_genes():
    """测试2: 能否正确根据 p-value 和 logFC 阈值筛选差异基因"""
    mock_df = pd.DataFrame({
        "Gene.symbol": ["G1", "G2", "G3", "G4"],
        "P.Value": [0.01, 0.20, 0.005, 0.04],
        "logFC": [2.0, 1.5, -0.5, -3.0]
    })
    # G1: 显著上调; G2: P值太大; G3: logFC太小; G4: 显著下调
    sig_genes = extract_significant_genes(mock_df, p_thresh=0.05, lfc_thresh=1.0)
    assert sig_genes == {"G1", "G4"}

def test_do_ora_empty_protection():
    """测试3: 当没有一个基因符合显著差异时，程序应安全触发警告保护并返回空表"""
    mock_df = pd.DataFrame({
        "Gene.symbol": ["G1", "G2"],
        "P.Value": [0.5, 0.8],
        "logFC": [0.1, 0.2]
    })
    gene_sets = {"PATHWAY_1": ["G1", "G2"]}
    
    res = do_ora(mock_df, gene_sets, p_thresh=0.05, lfc_thresh=1.0)
    assert res.empty
    assert list(res.columns) == ['pathway', 'P.Value', 'Count', 'pathway_size', 'GeneRatio', 'hits_genes']
# 3. 检验 gsea 模块中的简化 ES 计算逻辑
def test_calculate_es_logic():
    """测试4: 验证简化 GSEA 在极端的正/负富集情况下的得分方向"""
    # 构造一个已经降序排列好的假基因排序列表 (顶部全为正，底部全为负)
    ranked_list = pd.DataFrame({
        "Gene.symbol": ["MA", "MB", "MC", "MD"],
        "logFC": [3.0, 2.0, -1.0, -4.0]
    })
    
    # 情况 A: 如果通路基因全集中在顶部（正富集），ES 应该大于 0
    es_positive = calculate_es(ranked_list, {"MA", "MB"})
    assert es_positive > 0
    
    # 情况 B: 如果通路基因全集中在底部（负富集），ES 应该小于 0
    es_negative = calculate_es(ranked_list, {"MC", "MD"})
    assert es_negative < 0
    
    # 情况 C: 如果没有任何通路基因命中列表，ES 应该为 0.0
    es_zero = calculate_es(ranked_list, {"UNKNOWN_GENE"})
    assert es_zero == 0.0