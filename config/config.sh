#!/bin/bash
# ==============================================================================
# 基因组组装 Pipeline 配置文件
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 项目基础设置
# ------------------------------------------------------------------------------
# 项目名称 (用于输出文件命名)
export PROJECT_NAME="MyGenome"

# 工作目录 (默认为当前项目根目录)
export WORKDIR="$(pwd)"

# 线程数设置
export THREADS=2

# ------------------------------------------------------------------------------
# 2. 输入数据路径 (请根据实际情况修改)
# ------------------------------------------------------------------------------
# 原始数据目录
export DATA_DIR="${WORKDIR}/sequence_files"

# WGS 数据 (用于基因组调查)
export WGS_R1="${DATA_DIR}/WGS_R1.fastq"
export WGS_R2="${DATA_DIR}/WGS_R2.fastq"

# HiFi 数据 (用于组装)
export HIFI_READS="${DATA_DIR}/HiFi.fasta"

# Hi-C 数据 (用于挂载/Scaffolding)
export HIC_R1="${DATA_DIR}/Hi-C_R1.fastq"
export HIC_R2="${DATA_DIR}/Hi-C_R2.fastq"

# RNA-seq 数据 (可选)
export RNA_R1="${DATA_DIR}/dpil_RNA_seq_R1.fastq"
export RNA_R2="${DATA_DIR}/dpil_RNA_seq_R2.fastq"

# ------------------------------------------------------------------------------
# 3. 数据库路径
# ------------------------------------------------------------------------------
# Kraken2 数据库路径 (污染去除)
export KRAKEN_DB="${WORKDIR}/database/kraken_db"

# BUSCO 数据库路径 (质量评估)
# 例如: laurasiatheria_odb10
export BUSCO_DB="${WORKDIR}/database/busco_db"

# ------------------------------------------------------------------------------
# 4. 工具参数设置
# ------------------------------------------------------------------------------
# --- Fastp (质控) ---
export FASTP_MIN_LEN=145

# --- Jellyfish & GenomeScope (基因组调查) ---
export KMER_SIZE=19
export GENOME_SIZE_EST="1G" # 预估基因组大小，用于jellyfish hash size

# --- Kraken2 (污染去除) ---
export KRAKEN_CONFIDENCE=0.5
export KRAKEN_MIN_HIT_GROUPS=300

# --- Hifiasm (组装) ---
# 无特定参数，使用默认即可

# --- Chromap (Hi-C 比对) ---
# 预设模式
export CHROMAP_PRESET="hic"

# --- Juicer (可视化) ---
# Juicer tools jar 包路径 (安装脚本会自动下载)
export JUICER_JAR="${WORKDIR}/tools/juicer_tools.jar"
export JUICER_JVM_MEM="50G" # Java 内存设置

# ------------------------------------------------------------------------------
# 5. 系统路径设置
# ------------------------------------------------------------------------------
# Conda 环境名称
export CONDA_ENV_NAME="genome_assembly_env"
# 脚本日志目录
export LOG_DIR="${WORKDIR}/logs"

# ------------------------------------------------------------------------------
# 6. Pipeline 流程控制 (yes/no)
# ------------------------------------------------------------------------------
export RUN_FASTP="yes"      # 是否运行 Fastp 质控
export RUN_KRAKEN="yes"     # 是否运行 Kraken2 去污染
export RUN_BUSCO="yes"      # 是否运行 BUSCO 评估
export RUN_CHROMAP="yes"    # 是否运行 Chromap 
export RUN_YAHS="yes"       # 是否运行 Yahs 挂载
export RUN_JUICER="yes"     # 是否运行 Juicer 可视化

# ------------------------------------------------------------------------------
# 7. 输出目录设置
# ------------------------------------------------------------------------------
export FASTP_DIR="${WORKDIR}/01_fastp"
export SURVEY_DIR="${WORKDIR}/02_survey"
export KRAKEN_DIR="${WORKDIR}/03_kraken"
export HIFIASM_DIR="${WORKDIR}/04_hifiasm"
export BUSCO_DIR="${WORKDIR}/05_busco"
export CHROMAP_DIR="${WORKDIR}/06_chromap"
export YAHS_DIR="${WORKDIR}/07_yahs"
export JUICER_DIR="${WORKDIR}/08_juicer"

