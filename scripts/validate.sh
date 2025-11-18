#!/bin/bash
##############################################################################
# Pipeline Validation Script
# 流程验证脚本 - 检查环境和配置
##############################################################################

set -euo pipefail

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "=========================================="
echo "Genome Assembly Pipeline - Validation"
echo "=========================================="

# 检查函数
check_command() {
    local cmd=$1
    if command -v "$cmd" &> /dev/null; then
        echo -e "${GREEN}[✓]${NC} $cmd found: $(command -v $cmd)"
        return 0
    else
        echo -e "${RED}[✗]${NC} $cmd not found"
        
        # 添加特定工具的修复建议
        case "$cmd" in
            jellyfish)
                echo -e "  ${YELLOW}Suggestion:${NC} Run 'conda install -n genome_assembly -y -c bioconda jellyfish'"
                ;;
            fastp)
                echo -e "  ${YELLOW}Suggestion:${NC} Run 'conda install -n genome_assembly -y -c bioconda fastp'"
                ;;
            hifiasm)
                echo -e "  ${YELLOW}Suggestion:${NC} Run 'conda install -n genome_assembly -y -c bioconda hifiasm'"
                ;;
            samtools)
                echo -e "  ${YELLOW}Suggestion:${NC} Run 'conda install -n genome_assembly -y -c bioconda samtools'"
                ;;
        esac
        return 1
    fi
}

check_optional_command() {
    local cmd=$1
    if command -v "$cmd" &> /dev/null; then
        echo -e "${GREEN}[✓]${NC} $cmd found (optional)"
    else
        echo -e "${YELLOW}[!]${NC} $cmd not found (optional)"
    fi
}

check_file() {
    local file=$1
    if [[ -f "$file" ]]; then
        echo -e "${GREEN}[✓]${NC} File exists: $file"
        return 0
    else
        echo -e "${RED}[✗]${NC} File not found: $file"
        return 1
    fi
}

check_dir() {
    local dir=$1
    if [[ -d "$dir" ]]; then
        echo -e "${GREEN}[✓]${NC} Directory exists: $dir"
        return 0
    else
        echo -e "${YELLOW}[!]${NC} Directory not found: $dir"
        return 1
    fi
}

echo ""
echo "1. Checking required commands..."
echo "----------------------------"

MISSING_CMDS=0

# 必需的命令
for cmd in bash conda fastp jellyfish hifiasm samtools; do
    check_command "$cmd" || ((MISSING_CMDS++))
done

# 可选的命令
echo ""
echo "2. Checking optional commands..."
echo "----------------------------"

for cmd in kraken2 seqkit chromap yahs busco quast juicer genomescope2 yq; do
    check_optional_command "$cmd"
done

echo ""
echo "3. Checking configuration..."
echo "----------------------------"

CONFIG_FILE="${1:-config.yaml}"
if check_file "$CONFIG_FILE"; then
    if command -v yq &> /dev/null; then
        echo "  Project name: $(yq eval '.project.name' $CONFIG_FILE)"
        echo "  Threads: $(yq eval '.resources.threads' $CONFIG_FILE)"
    else
        echo -e "${YELLOW}[!]${NC} yq not found, cannot parse YAML"
    fi
fi

echo ""
echo "4. Checking system resources..."
echo "----------------------------"

# CPU核心数
CPU_CORES=$(nproc)
echo "  CPU cores: $CPU_CORES"

# 内存
TOTAL_MEM_GB=$(free -g | awk 'NR==2 {print $2}')
echo "  Total memory: ${TOTAL_MEM_GB}GB"

if [[ $TOTAL_MEM_GB -lt 64 ]]; then
    echo -e "${YELLOW}[!]${NC} Warning: Genome assembly typically requires 512GB+ memory"
fi

# 磁盘空间
AVAIL_SPACE_GB=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
echo "  Available disk space: ${AVAIL_SPACE_GB}GB"

if [[ $AVAIL_SPACE_GB -lt 1000 ]]; then
    echo -e "${YELLOW}[!]${NC} Warning: Genome assembly typically requires 10TB+ storage"
fi

echo ""  
echo "=========================================="
if [[ $MISSING_CMDS -eq 0 ]]; then
    echo -e "${GREEN}Validation passed!${NC}"
    echo "You can now run the pipeline."
else
    echo -e "${RED}Validation failed!${NC}"
    echo "Missing $MISSING_CMDS required command(s)."
    echo ""
    echo -e "${YELLOW}Fixing instructions:${NC}"
    echo "1. Make sure you have activated the environment: 'conda activate genome_assembly'"
    echo "2. Run the following command to reinstall missing tools:"
    echo "   conda install -n genome_assembly -y -c bioconda jellyfish fastp hifiasm samtools"
    echo "3. If still failing, try recreating the environment:"
    echo "   bash setup_environment.sh"
    echo ""
    echo "Please install missing software before running the pipeline."
    exit 1
fi
echo "=========================================="

echo ""
echo "Additional diagnostics:"
echo "----------------------"
# 检查conda环境路径是否正确
echo -e "${YELLOW}Current conda environment:${NC} $(conda info --envs | grep '*' | awk '{print $1}')"

# 检查PATH环境变量是否包含conda环境
if [[ "$PATH" == *"envs/genome_assembly"* ]]; then
    echo -e "${GREEN}[✓]${NC} Conda environment path found in PATH"
else
    echo -e "${YELLOW}[!]${NC} Warning: Conda environment path not in PATH"
    echo -e "  ${YELLOW}Suggestion:${NC} Make sure to activate the environment properly"
fi
