# 项目结构说明

本文档详细说明了基因组组装流程项目的文件组织结构。

## 📂 目录结构

```
genome-assembly-pipeline/
│
├── 📄 主要脚本
│   ├── assembly_pipeline.sh          # 主流程脚本（核心文件）
│   ├── setup_environment.sh          # 环境安装脚本
│   └── config.yaml                   # 配置文件（需根据实际情况修改）
│
├── 📋 配置文件
│   ├── config.yaml                   # 实际使用的配置文件
│   └── config.example.yaml           # 配置文件示例
│
├── 📚 文档
│   ├── README.md                     # 项目主文档
│   ├── QUICKSTART.md                 # 快速开始指南
│   ├── PROJECT_STRUCTURE.md          # 本文件
│   └── RequiredTechnology.md         # 技术栈要求说明
│
├── 🛠️ 辅助脚本 (scripts/)
│   ├── utils.sh                      # 工具函数库
│   └── validate.sh                   # 环境验证脚本
│
├── 📖 参考资料 (public/)
│   ├── assembly-参考文档.sh           # 参考脚本
│   └── AssemblyFlowchart.png         # 流程图
│
├── 🧬 测试数据 (sequence_files/)
│   ├── HiFi.fasta                    # HiFi测试数据
│   ├── Hi-C_R1.fastq                 # Hi-C测试数据（R1）
│   ├── Hi-C_R2.fastq                 # Hi-C测试数据（R2）
│   ├── WGS_R1.fastq                  # WGS测试数据（R1）
│   └── WGS_R2.fastq                  # WGS测试数据（R2）
│
└── 🗂️ 运行时生成的目录（自动创建）
    ├── logs/                         # 日志文件
    ├── 01_qc/                        # 质控结果
    ├── 02_genome_survey/             # 基因组调查结果
    ├── 03_kraken/                    # Kraken2污染检测结果
    ├── 04_hifiasm/                   # Hifiasm组装结果
    ├── 05_busco/                     # BUSCO评估结果
    ├── 06_quast/                     # QUAST评估结果
    ├── 07_chromap/                   # Chromap比对结果
    ├── 08_yahs/                      # YaHS挂载结果
    └── 09_juicer/                    # Juicer可视化文件
```

## 📝 核心文件说明

### 1. 主流程脚本

#### `assembly_pipeline.sh`
- **功能**: 主流程控制脚本
- **用途**: 协调所有组装步骤的执行
- **使用方法**:
  ```bash
  bash assembly_pipeline.sh run      # 运行完整流程
  bash assembly_pipeline.sh resume   # 从断点续跑
  bash assembly_pipeline.sh clean    # 清理中间文件
  ```

#### `setup_environment.sh`
- **功能**: 自动安装所需软件环境
- **用途**: 创建Conda环境并安装所有依赖软件
- **使用方法**:
  ```bash
  bash setup_environment.sh
  ```

### 2. 配置文件

#### `config.yaml`
- **功能**: 流程配置文件
- **内容**: 
  - 项目基本信息
  - 输入数据路径
  - 数据库路径
  - 计算资源配置
  - 各步骤参数设置
- **重要性**: ⭐⭐⭐⭐⭐ (必须根据实际情况修改)

#### `config.example.yaml`
- **功能**: 配置文件模板
- **用途**: 提供配置示例，复制后修改使用

### 3. 辅助脚本

#### `scripts/utils.sh`
- **功能**: 工具函数库
- **包含**:
  - 日志记录函数
  - 文件验证函数
  - 状态跟踪函数
  - 错误处理函数

#### `scripts/validate.sh`
- **功能**: 环境验证
- **检查内容**:
  - 必需命令是否存在
  - 系统资源是否充足
  - 配置文件是否正确

### 4. 文档

#### `README.md`
- 项目概述
- 安装指南
- 使用说明
- 故障排查

#### `QUICKSTART.md`
- 快速开始指南
- 最小配置示例
- 常见问题解答

#### `RequiredTechnology.md`
- 所需技术栈详细说明
- 硬件要求
- 软件依赖
- 学习路径建议

## 🔄 数据流程

```
输入数据
  ↓
config.yaml (配置)
  ↓
assembly_pipeline.sh (主流程)
  ↓
scripts/utils.sh (辅助函数)
  ↓
各步骤执行
  ↓
输出结果到编号目录
  ↓
logs/ (日志记录)
```

## 📊 输出目录说明

### `logs/`
- `assembly_*.log`: 主日志文件
- `status.log`: 步骤状态记录

### `01_qc/`
- 质控后的fastq文件
- fastp报告 (HTML/JSON)

### `02_genome_survey/`
- k-mer分析结果
- 基因组大小估算
- GenomeScope图表

### `03_kraken/`
- 物种分类报告
- 过滤后的HiFi数据

### `04_hifiasm/`
- **重要输出**: `<项目名>.contigs.fa` (组装的contigs)
- GFA图文件
- Hifiasm日志

### `05_busco/`
- BUSCO完整性评估报告
- 基因集完整度统计

### `06_quast/`
- 组装质量统计 (N50, N90等)
- 连续性评估报告

### `07_chromap/`
- Hi-C比对BAM文件
- Chromap索引

### `08_yahs/`
- **重要输出**: `<项目名>_scaffolds_final.agp` (染色体scaffold)
- YaHS bin文件

### `09_juicer/`
- **重要输出**: `<项目名>.hic` (JBAT可视化文件)
- Juicer中间文件

## 🎯 重要文件标识

### 必须修改的文件
- ✅ `config.yaml` - 必须根据实际数据路径修改

### 直接使用的文件
- ✅ `assembly_pipeline.sh` - 主流程脚本
- ✅ `setup_environment.sh` - 环境安装脚本
- ✅ `scripts/*.sh` - 辅助脚本

### 参考文档
- 📖 `README.md` - 主要参考
- 📖 `QUICKSTART.md` - 快速开始
- 📖 `RequiredTechnology.md` - 技术栈说明

## 💡 使用建议

1. **首次使用**:
   - 阅读 `QUICKSTART.md`
   - 运行 `setup_environment.sh`
   - 复制 `config.example.yaml` 为 `config.yaml`
   - 修改配置文件
   - 运行 `scripts/validate.sh` 验证环境

2. **正式运行**:
   - 准备好所有输入数据
   - 下载必需的数据库
   - 运行 `assembly_pipeline.sh run`

3. **监控进度**:
   - 查看 `logs/assembly_*.log`
   - 检查 `logs/status.log`

4. **结果分析**:
   - 查看 `05_busco/` 中的BUSCO报告
   - 查看 `06_quast/` 中的质量统计
   - 使用 `09_juicer/<项目名>.hic` 文件在JBAT中可视化

## 🔒 文件权限

所有 `.sh` 脚本已设置可执行权限：
```bash
chmod +x *.sh scripts/*.sh
```

## 📦 版本控制

项目使用Git进行版本控制：
- `.git/`: Git仓库
- `.gitignore`: 忽略文件配置

## ⚠️ 注意事项

1. **不要修改**:
   - `scripts/utils.sh` (除非你知道在做什么)
   - `assembly_pipeline.sh` 的核心逻辑

2. **可以修改**:
   - `config.yaml` (必须修改)
   - 各步骤的参数配置

3. **可以删除**:
   - `sequence_files/` (如果不需要测试数据)
   - `public/assembly-参考文档.sh` (参考脚本)

4. **不应该删除**:
   - 所有核心脚本
   - 配置文件
   - 文档文件

## 🆘 获取帮助

- 查看 `README.md` 中的故障排查部分
- 查看 `QUICKSTART.md` 中的常见问题
- 检查 `logs/` 目录中的日志文件
