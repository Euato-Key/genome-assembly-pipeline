# Genome Assembly Pipeline

一个模块化、可重复的染色体级别基因组组装流程，使用 PacBio HiFi 和 Hi-C 测序数据。

A modular and reproducible pipeline for chromosome-level genome assembly using PacBio HiFi and Hi-C data.

## 📋 项目概述

本项目实现了一个完整的基因组组装流程，包括：

- 📊 **数据质控**：使用 fastp 进行高质量数据过滤
- 🔬 **基因组调查**：使用 Jellyfish 和 GenomeScope2 估算基因组大小
- 🧹 **污染检测**：使用 Kraken2 检测并去除污染序列
- 🧬 **基因组组装**：使用 Hifiasm 进行高质量组装
- ✅ **质量评估**：使用 BUSCO 和 QUAST 评估组装质量
- 🧲 **染色体挂载**：使用 Chromap 和 YaHS 进行 Hi-C 挂载
- 👁️ **可视化**：生成 Juicer/JBAT 兼容文件用于手动校正

## 🎯 核心特性

- ✨ **模块化设计**：每个步骤独立运行，支持断点续跑
- 🔄 **可重复性**：完整的日志记录和状态追踪
- ⚙️ **灵活配置**：通过 YAML 配置文件轻松定制参数
- 📈 **进度监控**：实时日志和状态记录
- 🛠️ **环境管理**：基于 Conda 的环境管理，易于部署

## 📁 项目结构

```
genome-assembly-pipeline/
├── assembly_pipeline.sh      # 主流程脚本
├── config.yaml               # 配置文件
├── setup_environment.sh      # 环境安装脚本
├── scripts/
│   ├── utils.sh             # 工具函数库
│   └── validate.sh          # 环境验证脚本
├── public/
│   └── assembly-参考文档.sh  # 参考脚本
├── RequiredTechnology.md     # 技术栈说明
└── README.md                 # 本文档
```

## 🚀 快速开始

### 1. 环境准备

#### 系统要求

- **操作系统**: Linux (推荐 Ubuntu 20.04/22.04 或 CentOS 7/Rocky Linux 8+)
- **CPU**: 64核或更多
- **内存**: 512GB 起步，推荐 1TB+
- **存储**: 10TB+ 高速存储 (SSD/NVMe)

#### 安装 Conda 环境

```bash
# 1. 安装 Miniconda (如果尚未安装)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 2. 运行环境安装脚本
bash setup_environment.sh

# 3. 激活环境
conda activate genome_assembly

# 重要说明：运行主流程前必须确保genome_assembly环境已激活

# 4. 验证环境
bash scripts/validate.sh
```

### 2. 配置文件设置

编辑 `config.yaml` 文件，设置你的项目参数：

```yaml
# 项目基本信息
project:
  name: "MyGenome"              # 修改为你的项目名
  species: "species_name"       # 物种名称

# 输入数据路径
input:
  hifi_data: "/path/to/hifi.fasta"     # HiFi 数据
  hic_data: "/path/to/hic_raw"         # Hi-C 数据目录
  wgs_data: "/path/to/wgs_raw"         # WGS 数据（可选）

# 数据库路径
databases:
  kraken_db: "/path/to/kraken_db"      # Kraken2 数据库
  busco_db: "/path/to/busco_lineage"   # BUSCO 数据库

# 计算资源
resources:
  threads: 64                           # CPU 线程数
  memory_gb: 512                        # 内存大小
```

### 3. 运行流程

**重要：运行流程前必须确保已激活正确的环境！**

```bash
# 先确保环境已激活
conda activate genome_assembly

# 然后运行流程
# 运行完整流程
bash assembly_pipeline.sh run

# 从断点续跑
bash assembly_pipeline.sh resume

# 清理中间文件
bash assembly_pipeline.sh clean
```

## 📊 流程详解

### 流程图

```
HiFi 数据 ──┐                    ┌──> BUSCO 评估
            ├──> Kraken2 ──> Hifiasm ──┤
WGS 数据 ───┘    污染检测    组装      └──> QUAST 评估
                                              ↓
Hi-C 数据 ──> fastp ──> Chromap ──> YaHS ──> Juicer ──> JBAT 手动校正
             质控      比对       挂载     可视化
```

### 各步骤说明

#### 步骤 1: 数据质控 (01_qc)
- 使用 `fastp` 对 Hi-C、WGS 数据进行质控
- 过滤低质量 reads，去除接头序列
- 生成质控报告

#### 步骤 2: 基因组调查 (02_genome_survey)
- 使用 `jellyfish` 进行 k-mer 计数
- 使用 `GenomeScope2` 估算基因组大小、杂合度等参数

#### 步骤 3: 污染检测 (03_kraken)
- 使用 `Kraken2` 检测潜在污染
- 根据分类结果过滤序列
- 保留目标物种数据

#### 步骤 4: 基因组组装 (04_hifiasm)
- 使用 `Hifiasm` 进行 contig 组装
- 整合 Hi-C 数据（如果可用）
- 输出高质量 contigs

#### 步骤 5-6: 质量评估 (05_busco, 06_quast)
- `BUSCO`: 评估基因组完整性
- `QUAST`: 评估组装连续性 (N50, N90 等)

#### 步骤 7: Hi-C 比对 (07_chromap)
- 使用 `Chromap` 快速比对 Hi-C reads
- 生成 BAM 文件

#### 步骤 8: 染色体挂载 (08_yahs)
- 使用 `YaHS` 进行染色体级别挂载
- 生成 AGP 和 FASTA 文件

#### 步骤 9: 可视化准备 (09_juicer)
- 生成 Juicer 格式文件
- 创建 .hic 文件用于 JBAT 可视化和手动校正

## 🔧 进阶使用

### 自定义参数

在 `config.yaml` 中修改各步骤参数：

```yaml
steps:
  run_kraken: true                   # 是否运行 Kraken2
  kraken_confidence: 0.5             # Kraken2 置信度
  hifiasm_options: "--primary -l 1"  # Hifiasm 参数

qc_params:
  fastp:
    min_length: 145                  # 最小读长
  jellyfish:
    kmer_size: 19                    # k-mer 大小
```

### 跳过特定步骤

修改步骤控制参数：

```yaml
steps:
  run_qc: false           # 跳过质控（如果数据已质控）
  run_genome_survey: false # 跳过基因组调查
```

### 查看日志

所有日志保存在 `logs/` 目录：

```bash
# 查看主日志
tail -f logs/assembly_*.log

# 查看状态日志
cat logs/status.log
```

## 📚 所需数据库

### Kraken2 数据库

```bash
# 下载预构建数据库
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
tar -xzvf k2_standard_20230605.tar.gz -C /path/to/kraken_db/
```

### BUSCO 数据库

```bash
# 下载物种特异性数据库
wget https://busco-data.ezlab.org/v5/data/lineages/laurasiatheria_odb10.2024-01-08.tar.gz
tar -xzvf laurasiatheria_odb10.2024-01-08.tar.gz -C /path/to/busco_db/
```

## 🐛 故障排查

### 常见问题

**Q: 内存不足错误**
```bash
# 减少线程数
export THREADS=16

# 或在 config.yaml 中修改
resources:
  threads: 16
```

**Q: Conda 环境激活失败**
```bash
# 初始化 conda
conda init bash
source ~/.bashrc

# 重新激活
conda activate genome_assembly
```

**Q: 某个步骤失败需要重跑**
```bash
# 删除该步骤的状态记录
sed -i '/step_name/d' logs/status.log

# 重新运行
bash assembly_pipeline.sh resume
```

## 📖 参考文献

- **Hifiasm**: Cheng et al. (2021) Nature Methods
- **YaHS**: Zhou et al. (2023) Bioinformatics  
- **BUSCO**: Manni et al. (2021) Molecular Biology and Evolution
- **Chromap**: Zhang et al. (2021) Nature Communications

## 📄 许可证

MIT License

## 🤝 贡献

欢迎提交 Issue 和 Pull Request！

## 📧 联系方式

如有问题，请提交 Issue 或联系项目维护者。

---

**注意**: 本流程需要大量计算资源，建议在 HPC 集群或高性能服务器上运行。
