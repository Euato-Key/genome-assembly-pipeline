# 快速开始指南

本指南将帮助你在 10 分钟内开始运行基因组组装流程。

## 前置条件检查

确保你的系统满足以下要求：

```bash
# 检查系统资源
free -g                  # 至少 512GB 内存
nproc                    # 至少 64 核心
df -h                    # 至少 10TB 可用空间
```

## 步骤 1: 克隆或下载项目

```bash
cd /path/to/your/workspace
# 如果使用Git
git clone <repository-url> genome-assembly-pipeline
cd genome-assembly-pipeline

# 或者直接进入项目目录
cd genome-assembly-pipeline
```

## 步骤 2: 安装软件环境

```bash
# 运行环境安装脚本
bash setup_environment.sh

# 等待安装完成（约 10-30 分钟）
# 安装完成后激活环境
conda activate genome_assembly
```

## 步骤 3: 准备数据

### 3.1 组织你的数据文件

```bash
# 创建数据目录
mkdir -p data/raw/{hifi,hic,wgs}

# 将数据文件复制或链接到相应目录
# HiFi数据（FASTA格式）
cp /source/hifi.fasta data/raw/hifi/

# Hi-C数据（配对的FASTQ格式）
cp /source/hic/*R1*.fastq.gz data/raw/hic/
cp /source/hic/*R2*.fastq.gz data/raw/hic/

# WGS数据（可选）
cp /source/wgs/*.fastq.gz data/raw/wgs/
```

### 3.2 下载必需的数据库

```bash
# 创建数据库目录
mkdir -p databases/{kraken2,busco}

# 下载 Kraken2 标准数据库（约 60GB）
cd databases/kraken2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
tar -xzvf k2_standard_20230605.tar.gz
cd ../..

# 下载 BUSCO 数据库（选择合适的物种谱系）
cd databases/busco
# 对于哺乳动物
wget https://busco-data.ezlab.org/v5/data/lineages/laurasiatheria_odb10.2024-01-08.tar.gz
tar -xzvf laurasiatheria_odb10.2024-01-08.tar.gz
cd ../..
```

## 步骤 4: 配置流程参数

```bash
# 复制示例配置文件
cp config.example.yaml config.yaml

# 编辑配置文件
nano config.yaml  # 或使用 vim、vi 等编辑器
```

**最小必需配置**：

```yaml
project:
  name: "MyGenome"              # 修改为你的项目名

input:
  hifi_data: "data/raw/hifi/hifi.fasta"
  hic_data: "data/raw/hic"
  wgs_data: "data/raw/wgs"      # 可选

databases:
  kraken_db: "databases/kraken2/standard"
  busco_db: "databases/busco/laurasiatheria_odb10"

resources:
  threads: 64                    # 根据你的CPU核心数调整
```

## 步骤 5: 验证环境

```bash
# 运行验证脚本
bash scripts/validate.sh

# 应该看到绿色的 ✓ 标记表示检查通过
```

## 步骤 6: 运行流程

```bash
# 激活 conda 环境（如果尚未激活）
conda activate genome_assembly

# 开始运行完整流程
bash assembly_pipeline.sh run

# 流程将自动运行所有步骤，可能需要数小时到数天
```

## 步骤 7: 监控进度

在另一个终端窗口中：

```bash
# 实时查看日志
tail -f logs/assembly_*.log

# 查看步骤状态
cat logs/status.log

# 检查某个步骤的输出
ls -lh 04_hifiasm/
```

## 如果中断了怎么办？

流程支持断点续跑：

```bash
# 从上次中断的地方继续
bash assembly_pipeline.sh resume
```

## 结果文件位置

```
工作目录/
├── 04_hifiasm/
│   └── MyGenome.contigs.fa      # 组装的contigs
├── 05_busco/
│   └── MyGenome/                # BUSCO评估结果
├── 08_yahs/
│   └── MyGenome_scaffolds_final.agp  # 染色体级别scaffold
└── 09_juicer/
    └── MyGenome.hic             # 用于JBAT可视化
```

## 下一步

1. 查看 BUSCO 和 QUAST 评估结果
2. 使用 Juicebox Assembly Tools (JBAT) 进行手动校正
3. 查看详细文档：`README.md`

## 快速测试（使用小数据集）

如果你想先测试流程，可以使用小数据集：

```bash
# 创建测试配置
cp config.yaml config.test.yaml

# 修改为测试数据
nano config.test.yaml

# 运行测试
CONFIG_FILE=config.test.yaml bash assembly_pipeline.sh run
```

## 常见问题

**Q: 内存不足**
```bash
# 减少线程数和内存使用
# 在 config.yaml 中修改
resources:
  threads: 32
```

**Q: 缺少某个命令**
```bash
# 重新安装环境
conda env remove -n genome_assembly
bash setup_environment.sh
```

**Q: 想跳过某些步骤**
```bash
# 在 config.yaml 中设置
steps:
  run_genome_survey: false  # 跳过基因组调查
  run_kraken: false         # 跳过污染检测
```

需要帮助？查看完整文档或提交 Issue。
