# 1. 基础镜像：使用你指定的 Python 3.10 版本，-slim 表示精简版（体积小，无冗余系统库）
FROM python:3.10-slim

# 2. 设置工作目录为容器内的 /app
WORKDIR /app

# 3. 设置环境变量：防止 Python 打印日志缓冲，并在某些环境下强制 UTF-8 编码
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# 4. 升级基础 pip 包
RUN pip install --no-cache-dir --upgrade pip setuptools build

# 5. 为了利用 Docker 的层缓存（Layer Cache）加速后续构建，先仅复制项目元数据文件
COPY pyproject.toml README.md ./

# 6. 复制项目源代码
COPY src/ ./src/

# 7. 安装项目及其所有依赖项（numpy, pandas, scipy, matplotlib, seaborn）
# 这里 pip 会自动读取 pyproject.toml 并配置好终端命令 gsea-lite
RUN pip install --no-cache-dir .

# 8. 定义容器的主启动程序，对应 pyproject.toml 中的 [project.scripts]
ENTRYPOINT ["gsea-lite"]

# 9. 默认参数：如果用户只运行容器但不加参数，默认显示帮助文档
CMD ["--help"]