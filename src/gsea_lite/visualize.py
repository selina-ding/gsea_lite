import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def plot_barplot(enrichment_results: pd.DataFrame, top_n: int = 10, output_path: str = "barplot.png"):
    """绘制 Top N 通路的富集条形图"""
    if enrichment_results.empty:
        print("警告：没有富集结果可用于绘制柱状图。")
        return
    
    # 1. 排序与筛选：按 FDR 从小到大排序，截取 Top N
    top_df = enrichment_results.sort_values(by="FDR").head(top_n).copy()
    
    # 2. 数据转换：计算 -log10(FDR)，让越显著的通路柱子越长
    top_df['-log10(FDR)'] = -np.log10(top_df['FDR'] + 1e-10) # 加上一个极小值防止 log(0)
    
    # 3. 绘图核心：初始化画布，绘制横向条形图
    # 3. 绘图核心：初始化画布，绘制横向条形图
    plt.figure(figsize=(10, 6))
    # 根据警告提示，增加 hue='pathway' 并关闭冗余的图例 legend=False
    sns.barplot(x='-log10(FDR)', y='pathway', data=top_df, palette='Reds_r', hue='pathway', legend=False)
    
    # 4. 图表修饰与输出
    plt.title("Top Enriched Pathways")
    plt.xlabel(r"$-\log_{10}(\text{FDR})$")
    plt.ylabel("Pathway")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_dotplot(enrichment_results: pd.DataFrame, top_n: int = 10, output_path: str = "dotplot.png"):
    """绘制多维映射的气泡图"""
    if enrichment_results.empty:
        print("警告：没有富集结果可用于绘制气泡图。")
        return

    # 1. 排序与筛选
    top_df = enrichment_results.sort_values(by="FDR").head(top_n).copy()
    
    # 2. 绘图核心
    plt.figure(figsize=(10, 7))
    scatter = sns.scatterplot(
        x='GeneRatio',        # X轴：命中比例
        y='pathway',          # Y轴：通路名称
        size='Count',         # 气泡大小：命中基因数
        hue='FDR',            # 气泡颜色：FDR值
        data=top_df,
        palette='viridis',
        sizes=(50, 400)       # 控制气泡大小范围
    )
    
    # 3. 图表修饰与输出
    plt.title("Pathway Enrichment Dotplot")
    plt.xlabel("Gene Ratio")
    plt.ylabel("Pathway")
    
    # 移动图例到图表外部防止遮挡
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_gsea_barplot(gsea_results: pd.DataFrame, top_n: int = 10, output_path: str = "gsea_es_barplot.png"):
    """绘制简化 GSEA 结果中绝对 ES 最大的 Top N 通路条形图"""
    if gsea_results.empty:
        print("警告：没有 GSEA 结果可用于绘制条形图。")
        return

    required_columns = {"Pathway", "ES"}
    missing_columns = required_columns - set(gsea_results.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns for GSEA plot: {sorted(missing_columns)}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    top_gsea = (
        gsea_results.assign(abs_ES=lambda df: df["ES"].abs())
        .sort_values("abs_ES", ascending=False)
        .head(top_n)
        .sort_values("ES")
    )

    colors = np.where(top_gsea["ES"] >= 0, "#ba3b46", "#2f6f9f")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(top_gsea["Pathway"], top_gsea["ES"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Enrichment score (ES)")
    ax.set_ylabel("Pathway")
    ax.set_title("Top simplified GSEA results by absolute ES")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved GSEA plot to:", output_path)
