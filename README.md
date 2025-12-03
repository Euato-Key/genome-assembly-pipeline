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

## 项目特点
- **模块化设计**：各步骤解耦，便于维护和调试。
- **配置分离**：所有参数通过 `config/config.sh` 统一管理。
- **环境隔离**：基于 Conda 管理依赖，避免环境冲突。
- **日志记录**：全程记录运行日志，方便排查。
- **易于扩展**：代码结构清晰，方便后续封装为 Docker。

## 快速开始

### 1. 初始化环境

#### 一键安装（推荐）
```bash
bash scripts/setup_tools.sh
```
**脚本功能：**
- ✅ 自动检测并安装 Conda/Mamba（使用 Mamba 加速安装）
- ✅ 创建独立的运行环境 `genome_assembly_env`
- ✅ **自动编译 hifiasm**：从源码编译以解决云环境常见的 CPU 指令集兼容性问题（Illegal instruction）
- ✅ 自动配置 Shell 环境（conda init）

#### 快速安装（仅限云 IDE）
如果您在腾讯云 IDE 等环境中，也可以使用：
```bash
bash scripts/setup_tools_fast.sh
```

### 2. 激活环境
安装完成后，请重新加载 Shell 配置并激活环境：
```bash
source ~/.zshrc   # 如果使用 bash，请运行 source ~/.bashrc
conda activate genome_assembly_env
```

### 3. 环境验证
在正式运行之前，**强烈建议**运行验证脚本，检查工具是否安装正确且能正常运行：
```bash
bash scripts/validate_environment.sh
```
该脚本会检查：
- Conda 环境状态
- 核心工具（fastp, jellyfish, samtools 等）是否存在
- **hifiasm 是否能正常运行**（关键！）
- 可选工具状态

若脚本输出 `✅ 环境验证通过！`，说明一切就绪。

### 3. 修改配置
编辑 `config/config.sh`，设置您的数据路径和运行参数。
```bash
vim config/config.sh
```
> 重点修改 `HIFI_DATA`、`WGS_R1/R2`、`HIC_R1/R2`、`KRAKEN_DB`、`BUSCO_DB` 等路径。

### 4. 运行组装
```bash
bash scripts/run_assembly.sh
```
程序会自动按顺序执行：质控 → 基因组调查 → 组装 → 评估 → Hi‑C 比对 → Scaffolding → 可视化预处理。

### 5. 查看结果
- **日志**：`logs/` 目录下的时间戳日志记录每一步的详细信息。
- **最终组装**：`results/04_assembly/<PROJECT_NAME>.p_ctg.fa`（Contig）和 `results/07_yahs/<PROJECT_NAME>_scaffolded.fasta`（Scaffold）。
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
4. **Hifiasm** – 使用 HiFi 数据进行高质量 Contig 组装。
5. **BUSCO** –（可选）评估组装完整性。
6. **Chromap** – 将 Hi‑C 数据比对到组装好的 Contigs。
7. **YAHS** – 利用 Hi‑C 比对信息进行 Scaffolding，得到染色体级别的组装。
8. **Juicer** – 生成 `.hic` 文件，用于 Juicebox 可视化和人工纠错。

## 配置说明
请参考 `config/config.sh` 中的注释。主要配置项包括：
- `PROJECT_NAME`、`THREADS`、`HIFI_DATA`、`HIC_R1/R2`、`KRAKEN_DB`、`BUSCO_DB` 等。
- 详细的参数调节建议请查看 `docs/parameter_tuning.md`。

## 文档索引
- **[前置知识说明](docs/prerequisites.md)** – 了解组装原理、测序技术及常用术语。
- **[数据下载教程](docs/data_download.md)** – 学习如何获取公共数据库和测试数据。
- **[参数调节指南](docs/parameter_tuning.md)** – 根据不同基因组特性调节各工具参数的建议。
- **[hifiasm CPU 兼容性](docs/hifiasm_cpu_compatibility.md)** – 解决 CPU 指令集不兼容问题的详细指南。
- **[云 IDE 快速部署](docs/cloud_ide_setup.md)** – 腾讯云 IDE 等云环境的快速部署指南。

---
*Created by Antigravity Agent*
