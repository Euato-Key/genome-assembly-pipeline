# 基因组组装 Pipeline (Genome Assembly Pipeline)

这是一个基于 Bash 脚本的自动化基因组组装流程，专为 Ubuntu 24 环境设计。它集成了目前最主流的组装工具，支持从原始数据到染色体级别组装的全过程。

## 目录
- [项目特点](#项目特点)
- [快速开始](#快速开始)
- [目录结构](#目录结构)
- [流程详解](#流程详解)
- [配置说明](#配置说明)
- [文档索引](#文档索引)
- [环境验证脚本](#环境验证脚本)
- [常见问题](#常见问题)

## 项目特点
- **模块化设计**：各步骤解耦，便于维护和调试。
- **配置分离**：所有参数通过 `config/config.sh` 统一管理。
- **环境隔离**：基于 Conda 管理依赖，避免环境冲突。
- **日志记录**：全程记录运行日志，方便排查。
- **易于扩展**：代码结构清晰，方便后续封装为 Docker。

## 快速开始

### 1. 初始化环境
```bash
bash scripts/setup_tools.sh   # 安装 Conda、创建环境并下载必要工具
# 注意：jellyfish 有时需要手动安装一次（`conda install -c bioconda jellyfish`），
# 之后可通过 `bash scripts/validate_environment.sh` 确认所有工具已就绪。
```
如果需要下载测试数据进行演练：
```bash
bash scripts/setup_tools.sh --download-test-data
```

### 2. 环境验证（可选）
在正式运行之前，建议先执行环境验证脚本，检查所有必备工具、数据文件和数据库是否齐全：
```bash
bash scripts/validate_environment.sh
```
若脚本输出 `环境验证完成`，说明可以安全运行后续流程。

### 3. 修改配置
编辑 `config/config.sh`，设置您的数据路径和运行参数。
```bash
vim config/config.sh
```
> 重点修改 `HIFI_READS`、`WGS_R1/R2`、`HIC_R1/R2`、`KRAKEN_DB`、`BUSCO_DB` 等路径。

### 4. 运行组装
```bash
bash scripts/run_assembly.sh
```
程序会自动按顺序执行：质控 → 基因组调查 → 组装 → 评估 → Hi‑C 比对 → Scaffolding → 可视化预处理。

### 5. 查看结果
- **日志**：`logs/` 目录下的时间戳日志记录每一步的详细信息。
- **最终组装**：`results/04_assembly/<PROJECT_NAME>_flye/assembly.fasta`（Contig）和 `results/07_yahs/<PROJECT_NAME>_scaffolded.fasta`（Scaffold）。
- **BUSCO 报告**（若开启）：`results/05_busco/busco_result/short_summary.*`
- **Hi‑C 可视化**（若开启）：`results/08_juicer/` 生成的 `.hic` 文件可在 Juicebox 中查看。

## 目录结构
```
genome-assembly-pipeline/
├── config/                # 配置文件
│   └── config.sh
├── scripts/               # 脚本集合
│   ├── setup_tools.sh     # 环境安装与验证脚本
│   ├── run_assembly.sh    # 主流程脚本
│   └── validate_environment.sh  # 环境验证脚本（新增）
├── docs/                  # 文档说明
│   ├── prerequisites.md   # 前置知识说明
│   ├── data_download.md   # 数据下载教程
│   └── parameter_tuning.md# 参数调节指南
├── logs/                  # 运行日志（自动生成）
├── results/               # 分析结果（自动生成）
├── tools/                 # Conda 环境及工具安装目录
├── sequence_files/        # 原始测序数据存放目录（建议）
└── README.md              # 项目说明文档（此文件）
```

## 流程详解
本 Pipeline 包含以下步骤：
1. **Fastp** – 对 WGS 与 Hi‑C 数据进行质控。
2. **Jellyfish + GenomeScope** – 基于 WGS 数据进行 K‑mer 计数，估算基因组大小与杂合度。
3. **Kraken2** –（可选）检测并过滤可能的污染序列。
4. **Flye** – 使用 HiFi 数据进行高质量 Contig 组装。
5. **BUSCO** –（可选）评估组装完整性。
6. **Chromap** – 将 Hi‑C 数据比对到组装好的 Contigs。
7. **YAHS** – 利用 Hi‑C 比对信息进行 Scaffolding，得到染色体级别的组装。
8. **Juicer** – 生成 `.hic` 文件，用于 Juicebox 可视化和人工纠错。

## 配置说明
请参考 `config/config.sh` 中的注释。主要配置项包括：
- `PROJECT_NAME`、`THREADS`、`HIFI_READS`、`HIC_R1/R2`、`KRAKEN_DB`、`BUSCO_DB` 等。
- 详细的参数调节建议请查看 `docs/parameter_tuning.md`。

## 文档索引
- **[前置知识说明](docs/prerequisites.md)** – 了解组装原理、测序技术及常用术语。
- **[数据下载教程](docs/data_download.md)** – 学习如何获取公共数据库和测试数据。
- **[参数调节指南](docs/parameter_tuning.md)** – 根据不同基因组特性调节各工具参数的建议。

---

## 更新日志

### 2025-12-02
- 将组装器从Hifiasm更改为Flye，解决在某些环境中出现的"Illegal instruction"错误
- 添加了HiFi数据格式转换的解决方案
- 更新了组装结果路径说明
- 添加了常见问题解答部分

## 常见问题

### Q: 为什么使用Flye而不是Hifiasm进行组装？
A: 在某些环境中，Hifiasm可能会出现"Illegal instruction"错误，这通常是由于CPU架构不兼容或编译问题导致的。Flye是一个优秀的替代方案，特别适合处理HiFi数据，并且在各种环境中都有良好的兼容性。

### Q: 如何评估组装质量？
A: 可以通过以下方式评估组装质量：
1. 使用BUSCO评估基因完整性
2. 检查N50和L50统计值
3. 查看总组装长度与预估基因组大小的比较
4. 分析contig数量和长度分布

### Q: 如果HiFi数据格式不正确怎么办？
A: 如果HiFi数据没有标准的FASTA头部，可以使用以下Python脚本进行转换：
```python
# 读取原始文件
with open('原始文件.fasta', 'r') as f:
    content = f.read()

# 创建新的FASTA文件
with open('转换后文件.fasta', 'w') as out:
    lines = content.split('\n')
    seq_count = 0
    for line in lines:
        if not line.strip():
            continue
        if line.startswith('m84208'):  # 根据实际情况修改
            parts = line.split('\t')
            if len(parts) >= 10:
                seq = parts[9]
                out.write(f'>read_{seq_count}\n')
                out.write(seq + '\n')
                seq_count += 1
```

### Q: 如何解决Kraken2数据库未找到的问题？
A: 需要下载并配置Kraken2数据库。可以执行以下步骤：
1. 创建数据库目录：`mkdir -p /workspace/database/kraken_db`
2. 下载标准数据库：`kraken2-build --download-library bacteria --download-library viral --download-library human --db /workspace/database/kraken_db`
3. 构建数据库：`kraken2-build --build --db /workspace/database/kraken_db`

---
*Created by Antigravity Agent*
