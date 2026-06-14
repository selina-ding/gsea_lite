"""
gsea-lite: 轻量级基因集富集分析工具

本项目包含了两种经典的富集分析方法：过表达分析 (ORA) 和 简化版基因集富集分析 (GSEA)。
提供从数据读取、超几何检验、FDR多重假设检验校正到结果可视化的完整工作流。
"""

__version__ = "0.1.0"
__author__ = "丁奕新，夏奕琳婷，于润湉"

# 1. 数据读取模块
from .io import read_gmt, read_geo

# 2. 数据预处理模块 (多重检验)
from .stats import bh_fdr_correction

# 3. ORA 分析模块
from .ora import extract_significant_genes, do_ora

# 4. GSEA 分析模块
from .gsea import generate_ranked_list, calculate_es, run_gsea

# 5. 可视化模块
from .visualize import plot_barplot, plot_dotplot, plot_gsea_barplot

__all__ = [
    # IO
    "read_gmt",
    "read_geo",
    
    # Stats
    "bh_fdr_correction",
    
    # ORA
    "extract_significant_genes",
    "do_ora",
    
    # GSEA
    "generate_ranked_list",
    "calculate_es",
    "run_gsea",
    
    # Plotting
    "plot_barplot",
    "plot_dotplot",
    "plot_gsea_barplot"
]
