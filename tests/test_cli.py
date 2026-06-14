import pytest
import subprocess
from pathlib import Path

def test_cli_help_message():
    """测试1: 在终端输入 gsea-lite -h 能够正确返回帮助文档"""
    result = subprocess.run(["gsea-lite", "-h"], capture_output=True, text=True)
    
    # 断言执行状态为 0 (正常退出)
    assert result.returncode == 0
    # 断言打印出来的帮助信息里包含关键参数说明
    assert "gsea-lite" in result.stdout
    assert "ora" in result.stdout
    assert "gsea" in result.stdout

def test_cli_missing_args():
    """测试2: 运行子命令时如果缺失必需参数，必须报错拦截"""
    # 运行 ora 但不给任何文件路径
    result = subprocess.run(["gsea-lite", "ora"], capture_output=True, text=True)
    
    # 状态码不为 0 说明被成功拦截报错
    assert result.returncode != 0
    assert "the following arguments are required" in result.stderr

def test_cli_file_not_found():
    """测试3: 当输入不存在的假文件路径时，系统能够提示文件找不到而安全退出"""
    result = subprocess.run(
        ["gsea-lite", "ora", "--geo", "non_existent_geo.tsv", "--gmt", "non_existent.gmt"], 
        capture_output=True, 
        text=True
    )
    # cli.py 里的 if not geo_file.exists(): return 机制测试
    assert "❌ 错误: 找不到 GEO 文件" in result.stdout
