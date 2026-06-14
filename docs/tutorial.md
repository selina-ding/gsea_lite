# 🚀 gsea-lite 使用教程 (Tutorial)

欢迎使用 `gsea-lite`！本工具专为 BIO2502 课程设计，旨在提供一个轻量级、易于理解的基因集富集分析工作流。本教程将引导你从安装到完成第一次分析。

## 📦 1. 安装指南

我们推荐在专用的 Python 虚拟环境（如 Conda）中安装和使用本工具。

确保你已经下载了本项目的源代码，并在终端中进入了项目根目录（即 `pyproject.toml` 所在的目录）：

```bash
# 进入项目目录
cd path/to/gsea-lite

# 以开发者模式安装包及其依赖项 (包括绘图和测试工具)
pip install -e .[dev]
```

安装完成后，你可以在终端直接输入 `gsea-lite -h` 来验证是否安装成功。

## 🛠️ 2. 准备数据

在运行分析之前，你需要准备两个标准格式的文件：

1. **表达量差异分析结果 (DEG 表格)**：一个以制表符分隔的 `.tsv` 或 `.csv` 文件。**必须**包含以下列：
   - `Gene.symbol`: 基因名称（如 TP53, GAPDH）
   - `logFC`: Log2 转换后的差异倍数
   - `P.Value`: 差异表达的显著性 P 值
2. **基因集文件 (GMT 格式)**：可以从 MSigDB 或其他数据库下载的通路文件（例如 `h.all.v2026.1.Hs.symbols.gmt`）。

*提示：在项目的 `examples/example_data/` 目录下，我们为你提供了一组测试数据用于练习。*

## 💻 3. 命令行 (CLI) 快速入门

`gsea-lite` 提供了极其优雅的命令行接口，让你无需编写代码即可快速出图。

### 模式 A: 过表达分析 (ORA)

ORA 分析会根据你设定的阈值筛选出显著差异基因，并计算它们在各个通路中的富集情况。

**基础命令：**

Bash

```
gsea-lite ora --geo examples/example_data/GSE42568.top.table.tsv \
              --gmt examples/example_data/h.all.v2026.1.Hs.symbols.gmt \
              --outdir ./my_results
```

**进阶用法（自定义阈值）：** 你可以通过参数调整筛选差异基因的严格程度：

Bash

```
# 设定 p-value < 0.01 且 |logFC| > 1.5
gsea-lite ora --geo input.tsv --gmt input.gmt --p-thresh 0.01 --lfc-thresh 1.5 --outdir ./results
```

**输出结果**：程序将在 `outdir` 目录下生成定量的 CSV 结果表 `ora_results_fdr.csv`，并自动绘制 `ora_barplot.png` (条形图) 和 `ora_dotplot.png` (气泡图)。

### 模式 B: 简化版基因集富集分析 (GSEA)

GSEA 不需要硬性阈值，它通过评估基因在排序列表中的整体分布来计算富集分数 (Enrichment Score)。

**运行命令：**

Bash

```
gsea-lite gsea --geo examples/example_data/GSE42568.top.table.tsv \
               --gmt examples/example_data/h.all.v2026.1.Hs.symbols.gmt \
               --outdir ./my_results
```

**输出结果**：程序将在指定目录下生成 `gsea_results.csv`，其中包含了按 ES 绝对值降序排列的通路富集结果。

## 🐍 4. 在 Python / Jupyter Notebook 中使用 (API)

如果你希望在 Python 脚本中进行更灵活的探索（例如你们的 `demo.ipynb`），可以直接导入我们的核心模块。

Python

```
import pandas as pd
from pathlib import Path
from gsea_lite import (
    read_gmt, read_geo, 
    do_ora, bh_fdr_correction, 
    run_gsea, 
    plot_barplot, plot_dotplot
)

# 1. 加载数据
geo_data = read_geo("examples/example_data/GSE42568.top.table.tsv")
gene_sets = read_gmt("examples/example_data/h.all.v2026.1.Hs.symbols.gmt")

# 2. 运行 ORA 分析
print("Running ORA...")
ora_res = do_ora(geo_data, gene_sets, p_thresh=0.05, lfc_thresh=1.0)

# 3. 添加 FDR 校正
ora_res['FDR'] = bh_fdr_correction(ora_res['P.Value'])

# 4. 可视化
plot_dotplot(ora_res, top_n=10, output_path="custom_dotplot.png")
plot_barplot(ora_res, top_n=10, output_path="custom_barplot.png")

# 5. 运行 GSEA 分析
print("Running GSEA...")
gsea_res = run_gsea(geo_data, gene_sets)
display(gsea_res.head())
```

## ❓ 5. 常见问题 (FAQ)

**Q: 为什么运行 ORA 提示“未找到显著差异基因”？**

> A: 可能是因为你的数据差异不明显，导致没有基因满足默认的 `P.Value < 0.05` 且 `|logFC| > 1.0` 阈值。你可以尝试在命令中加上 `--p-thresh 0.1 --lfc-thresh 0.5` 放宽条件试试。

**Q: 我的输入文件没有 logFC 列怎么办？**

> A: 本工具要求输入数据必须经过标准的差异表达分析（如 DESeq2, edgeR 或 limma）。请确保在预处理阶段生成了这些统计指标。
