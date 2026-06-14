import argparse
from pathlib import Path

# 导入对应的函数
from gsea_lite import (
    read_gmt, read_geo, 
    do_ora, bh_fdr_correction, 
    run_gsea, 
    plot_barplot, plot_dotplot
)

def main():
    # 1. 创建主解析器
    parser = argparse.ArgumentParser(
        prog="gsea-lite",
        description="简化的基因集富集分析 (GSEA) 与过表达分析 (ORA) 工具"
    )

    # 2. 创建子命令解析器
    subparsers = parser.add_subparsers(
        title="分析模式",
        dest="command",
        help="请选择要运行的分析类型: ora 或 gsea"
    )
    subparsers.required = True 

    # 子命令 1: ORA (过表达分析)
    parser_ora = subparsers.add_parser("ora", help="运行 ORA 分析与可视化")
    parser_ora.add_argument("--geo", required=True, help="输入的 GEO 差异表达分析结果文件路径 (例如 GSE42568.top.table.tsv)")
    parser_ora.add_argument("--gmt", required=True, help="输入的 GMT 格式基因集文件路径")
    parser_ora.add_argument("--outdir", default="./results", help="输出结果的目录 (默认: ./results)")
    parser_ora.add_argument("--p-thresh", type=float, default=0.05, help="ORA P-value 阈值 (默认: 0.05)")
    parser_ora.add_argument("--lfc-thresh", type=float, default=1.0, help="ORA logFC 阈值 (默认: 1.0)")

    # 子命令 2: GSEA (基因集富集分析)
    parser_gsea = subparsers.add_parser("gsea", help="运行简化版 GSEA 分析")
    parser_gsea.add_argument("--geo", required=True, help="输入的 GEO 差异表达分析结果文件路径")
    parser_gsea.add_argument("--gmt", required=True, help="输入的 GMT 格式基因集文件路径")
    parser_gsea.add_argument("--outdir", default="./results", help="输出结果的目录 (默认: ./results)")

    # 3. 解析终端参数
    args = parser.parse_args()

    # 4. 环境准备：确保输出目录存在
    output_dir = Path(args.outdir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    geo_file = Path(args.geo).resolve()
    gmt_file = Path(args.gmt).resolve()

    # 提前校验文件是否存在
    if not geo_file.exists():
        print(f"❌ 错误: 找不到 GEO 文件 {geo_file}")
        return
    if not gmt_file.exists():
        print(f"❌ 错误: 找不到 GMT 文件 {gmt_file}")
        return

    # 5. 执行具体的命令逻辑
    if args.command == "ora":
        run_ora_workflow(geo_file, gmt_file, output_dir, args.p_thresh, args.lfc_thresh)
    elif args.command == "gsea":
        run_gsea_workflow(geo_file, gmt_file, output_dir)


def run_ora_workflow(geo_file, gmt_file, output_dir, p_thresh, lfc_thresh):
    """提取出来的 ORA 业务流"""
    print(f"🚀 开始执行 ORA 分析流程...")
    print(f"📂 读取 GEO 数据: {geo_file.name}")
    print(f"📂 读取 GMT 数据: {gmt_file.name}")
    
    geo_data = read_geo(geo_file)
    gene_sets = read_gmt(gmt_file)
    
    print(f"\n🔬 使用阈值 P.Value < {p_thresh}, |logFC| > {lfc_thresh} 进行 ORA 分析...")
    ora_results = do_ora(geo_data, gene_sets, p_thresh=p_thresh, lfc_thresh=lfc_thresh)
    
    if not ora_results.empty:
        # 添加 FDR 校正
        ora_results['FDR'] = bh_fdr_correction(ora_results['P.Value'])
        cols = ['pathway', 'P.Value', 'FDR', 'Count', 'pathway_size', 'GeneRatio', 'hits_genes']
        ora_results = ora_results[cols]
        
        # 导出结果
        ora_csv_path = output_dir / "ora_results_fdr.csv"
        ora_results.to_csv(ora_csv_path, index=False)
        print(f"✅ ORA 分析完成，表格已保存至: {ora_csv_path}")
        
        # 绘图
        print("📊 正在绘制 ORA 可视化图表...")
        plot_barplot(ora_results, top_n=10, output_path=output_dir / "ora_barplot.png")
        plot_dotplot(ora_results, top_n=10, output_path=output_dir / "ora_dotplot.png")
        print("✅ 绘图完成！")
    else:
        print("⚠️ 未找到显著差异基因，跳过导出和绘图。")


def run_gsea_workflow(geo_file, gmt_file, output_dir):
    """提取出来的 GSEA 业务流"""
    print(f"🚀 开始执行简化版 GSEA 分析流程...")
    
    geo_data = read_geo(geo_file)
    gene_sets = read_gmt(gmt_file)
    
    gsea_results = run_gsea(geo_data, gene_sets)
    
    if not gsea_results.empty:
        gsea_csv_path = output_dir / "gsea_results.csv"
        gsea_results.to_csv(gsea_csv_path, index=False)
        print(f"✅ GSEA 分析完成，表格已保存至: {gsea_csv_path}")
    else:
        print("⚠️ GSEA 结果为空。")


if __name__ == "__main__":
    main()