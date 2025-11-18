#!/bin/bash
##############################################################################
# Environment Setup Script for Genome Assembly Pipeline
# 基因组组装流程环境安装脚本
##############################################################################

set -e

echo "=========================================="
echo "Genome Assembly Pipeline - Environment Setup"
echo "=========================================="

# 检查Conda是否安装
if ! command -v conda &> /dev/null; then
    echo "Error: Conda not found. Please install Miniconda or Anaconda first."
    echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# 环境名称和路径
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${PROJECT_DIR}/env"
ENV_NAME="genome_assembly"

# 检查环境是否已存在
if [[ -d "${ENV_DIR}" ]]; then
    echo "Environment already exists at: ${ENV_DIR}"
    read -p "Do you want to remove and recreate it? (y/n): " choice
    if [[ "$choice" == "y" ]]; then
        echo "Removing existing environment..."
        rm -rf "${ENV_DIR}"
    else
        echo "Exiting..."
        exit 0
    fi
fi

echo "Creating conda environment at: ${ENV_DIR}..."

# 创建环境并安装软件（使用-p指定路径）
conda create -p "${ENV_DIR}" -y -c conda-forge -c bioconda \
    python=3.9 \
    fastp \
    jellyfish \
    kraken2 \
    seqkit \
    hifiasm \
    chromap \
    yahs \
    busco=5.4.7 \
    quast \
    samtools \
    bcftools \
    parallel \
    pigz \
    yq

echo ""
echo "Installing additional Python packages..."
# 使用conda run在指定环境中运行pip
conda run -p "${ENV_DIR}" pip install biopython pandas matplotlib seaborn pyyaml

echo ""
echo "=========================================="
echo "Environment setup completed!"
echo "=========================================="
echo ""
echo "Environment location: ${ENV_DIR}"
echo ""
echo "To activate the environment, run:"
echo "  conda activate ${ENV_DIR}"
echo "  # or"
echo "  source activate ${ENV_DIR}"
echo ""
echo "To deactivate, run:"
echo "  conda deactivate"
echo ""
echo "Additional tools that need manual installation:"
echo "  1. GenomeScope2: https://github.com/tbenavi1/genomescope2.0"
echo "  2. Juicer Tools: https://github.com/aidenlab/juicer/wiki/Download"
echo "  3. Juicebox Assembly Tools (JBAT): https://github.com/aidenlab/Juicebox"
echo ""
echo "Database preparation:"
echo "  1. Kraken2 database: https://benlangmead.github.io/aws-indexes/k2"
echo "  2. BUSCO database: https://busco-data.ezlab.org/v5/data/lineages/"
echo ""
