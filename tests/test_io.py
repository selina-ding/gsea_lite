import pytest
import pandas as pd
from gsea_lite.io import read_gmt, read_geo

# 1. 测试 read_gmt 函数
def test_read_gmt_normal(tmp_path):
    """测试1: 能够正确解析标准的 GMT 格式文件"""
    # 动态创建一个临时的假 gmt 文件
    gmt_file = tmp_path / "test.gmt"
    gmt_content = (
        "PATHWAY_A\tdesc_a\tGENE1\tGENE2\tGENE3\n"
        "PATHWAY_B\tdesc_b\tGENE2\tGENE4\n"
    )
    gmt_file.write_text(gmt_content, encoding="utf-8")
    
    # 运行函数并断言
    gene_sets = read_gmt(gmt_file)
    assert len(gene_sets) == 2
    assert gene_sets["PATHWAY_A"] == ["GENE1", "GENE2", "GENE3"]
    assert gene_sets["PATHWAY_B"] == ["GENE2", "GENE4"]

def test_read_gmt_with_empty_and_spaces(tmp_path):
    """测试2: 测试 GMT 文件中包含空行和多余空格时的清洗能力"""
    gmt_file = tmp_path / "spaces.gmt"
    gmt_content = (
        "PATHWAY_C\tdesc_c\tGENE1  \t  GENE2\n"  # 基因前后带空格
        "\n"                                     # 存在空行
        "PATHWAY_D\tdesc_d\t\n"                  # 通路不含任何基因
    )
    gmt_file.write_text(gmt_content, encoding="utf-8")
    
    gene_sets = read_gmt(gmt_file)
    assert "PATHWAY_C" in gene_sets
    assert gene_sets["PATHWAY_C"] == ["GENE1", "GENE2"]  # 应该被 strip 干净
    assert gene_sets["PATHWAY_D"] == []                  # 空列表

# 2. 测试 read_geo 函数
def test_read_geo_normal(tmp_path):
    """测试3: 能够正确读取标准的 GEO 差异表达 TSV 文件"""
    geo_file = tmp_path / "geo.tsv"
    geo_content = "Gene.symbol\tlogFC\tP.Value\nGAPDH\t2.5\t0.001\nTP53\t-1.2\t0.04\n"
    geo_file.write_text(geo_content, encoding="utf-8")
    
    df = read_geo(geo_file)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert "Gene.symbol" in df.columns

def test_read_geo_missing_column(tmp_path):
    """测试4: 如果缺少核心列(如 logFC),必须抛出 ValueError 异常"""
    geo_file = tmp_path / "bad_geo.tsv"
    geo_content = "Gene.symbol\tP.Value\nGAPDH\t0.001\n"  # 缺少 logFC 列
    geo_file.write_text(geo_content, encoding="utf-8")
    
    with pytest.raises(ValueError, match="error: 缺失列 logFC"):
        read_geo(geo_file)