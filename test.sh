#!/bin/bash
##############################################################################
# Quick Test Script
# 快速测试脚本 - 使用小数据集测试流程
##############################################################################

set -e

echo "=========================================="
echo "Genome Assembly Pipeline - Quick Test"
echo "=========================================="
echo ""

# 颜色定义
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# 检查是否在项目根目录
if [[ ! -f "assembly_pipeline.sh" ]]; then
    echo -e "${RED}Error: Please run this script from the project root directory${NC}"
    exit 1
fi

echo "Step 1: Checking test data..."
echo "------------------------------"

# 检查测试数据是否存在
if [[ -d "sequence_files" ]]; then
    echo -e "${GREEN}✓${NC} Test data found in sequence_files/"
    TEST_DATA_DIR="sequence_files"
else
    echo -e "${YELLOW}!${NC} No test data found, this is a dry run test"
    TEST_DATA_DIR=""
fi

echo ""
echo "Step 2: Checking environment..."
echo "------------------------------"

# 检查conda
if command -v conda &> /dev/null; then
    echo -e "${GREEN}✓${NC} Conda found"
else
    echo -e "${RED}✗${NC} Conda not found"
    exit 1
fi

# 检查环境是否存在
if conda env list | grep -q "genome_assembly"; then
    echo -e "${GREEN}✓${NC} genome_assembly environment found"
else
    echo -e "${YELLOW}!${NC} genome_assembly environment not found"
    echo "Please run: bash setup_environment.sh"
    exit 1
fi

echo ""
echo "Step 3: Creating test configuration..."
echo "------------------------------"

# 创建测试配置
cat > config.test.yaml << 'EOF'
project:
  name: "TestGenome"
  species: "test_species"
  description: "Quick test run"

resources:
  threads: 4
  memory_gb: 16

input:
  hifi_data: "sequence_files/HiFi.fasta"
  hic_data: "sequence_files"
  wgs_data: ""
  rnaseq_data: ""

databases:
  kraken_db: ""
  mammal_taxids: ""
  busco_db: ""
  busco_lineage: ""

software:
  use_conda: true
  conda_env: "genome_assembly"
  use_docker: false
  docker_image: ""

steps:
  run_qc: false
  run_fastp: false
  run_genome_survey: false
  run_kraken: false
  run_hifiasm: false
  run_busco: false
  run_quast: false
  run_chromap: false
  run_yahs: false
  run_juicer: false

qc_params:
  fastp:
    min_length: 100
    threads: 4
  jellyfish:
    kmer_size: 19
    hash_size: "100000000"
  genomescope:
    ploidy: 2

output:
  work_dir: "test_run"
  keep_intermediate: true
  log_dir: "test_run/logs"
EOF

echo -e "${GREEN}✓${NC} Test configuration created: config.test.yaml"

echo ""
echo "Step 4: Running validation..."
echo "------------------------------"

# 运行验证
bash scripts/validate.sh config.test.yaml || {
    echo -e "${YELLOW}!${NC} Validation warnings detected, but continuing..."
}

echo ""
echo "=========================================="
echo -e "${GREEN}Quick test completed!${NC}"
echo "=========================================="
echo ""
echo "Test configuration saved to: config.test.yaml"
echo ""
echo "To run a minimal test (if you have test data):"
echo "  1. Edit config.test.yaml and enable desired steps"
echo "  2. Run: CONFIG_FILE=config.test.yaml bash assembly_pipeline.sh run"
echo ""
echo "To run the full pipeline:"
echo "  1. Edit config.yaml with your real data paths"
echo "  2. Run: bash assembly_pipeline.sh run"
echo ""

# 清理
echo "Cleaning up test configuration..."
# rm -f config.test.yaml  # 保留配置文件供参考

echo "Test configuration preserved for reference."
