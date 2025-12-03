#!/bin/bash
# ==============================================================================
# 脚本名称: validate_environment.sh
# 功能: 验证基因组组装流程的运行环境
# ==============================================================================

set -e

# 获取脚本所在目录的绝对路径
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="${PROJECT_ROOT}/config/config.sh"

# 加载配置
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "错误: 找不到配置文件 $CONFIG_FILE"
    exit 1
fi

# 日志函数
log_info() {
    echo -e "\033[32m[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1\033[0m"
}

log_error() {
    echo -e "\033[31m[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1\033[0m"
}

log_warn() {
    echo -e "\033[33m[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1\033[0m"
}

# 1. 检查 Conda 环境
log_info "=== 1. 检查 Conda 环境 ==="
if ! command -v conda &> /dev/null; then
    log_error "Conda 未安装！"
    exit 1
fi

# 检查环境是否存在
if conda info --envs | grep -q "${CONDA_ENV_NAME}"; then
    log_info "Conda 环境 '${CONDA_ENV_NAME}' 存在。"
else
    log_error "Conda 环境 '${CONDA_ENV_NAME}' 不存在！请先运行 scripts/setup_tools.sh"
    exit 1
fi

# 检查当前是否激活了正确的环境
CURRENT_ENV=$(conda info --envs | grep "*" | awk '{print $1}')
if [ "$CURRENT_ENV" != "$CONDA_ENV_NAME" ]; then
    log_warn "当前激活的环境是 '$CURRENT_ENV'，而不是 '$CONDA_ENV_NAME'。"
    log_warn "请运行: conda activate $CONDA_ENV_NAME"
    # 尝试在脚本中临时激活以进行后续检查
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV_NAME}"
fi

# 2. 检查关键工具
log_info "=== 2. 检查关键工具 ==="
TOOLS=(
    "fastp"
    "jellyfish"
    "samtools"
    "chromap"
    "yahs"
    "make"
    "pigz"
)

MISSING_TOOLS=0
for tool in "${TOOLS[@]}"; do
    if command -v "$tool" &> /dev/null; then
        echo -e "  [OK] $tool: $( $tool --version 2>&1 | head -n 1 | cut -c 1-50 )..."
    else
        echo -e "  \033[31m[FAIL] $tool 未找到\033[0m"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
    fi
done

# 3. 特别检查 hifiasm (编译版)
log_info "=== 3. 检查 hifiasm (编译版) ==="
if command -v hifiasm &> /dev/null; then
    HIFIASM_VER=$(hifiasm --version 2>&1)
    if [ $? -eq 0 ]; then
        echo -e "  [OK] hifiasm: $HIFIASM_VER"
        log_info "hifiasm 可正常运行 (CPU 指令集兼容)"
    else
        echo -e "  \033[31m[FAIL] hifiasm 无法运行 (可能是 CPU 指令集不兼容)\033[0m"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
    fi
else
    echo -e "  \033[31m[FAIL] hifiasm 未找到\033[0m"
    MISSING_TOOLS=$((MISSING_TOOLS + 1))
fi

# 4. 检查可选工具
log_info "=== 4. 检查可选工具 ==="
OPTIONAL_TOOLS=("genomescope2" "kraken2" "busco")
for tool in "${OPTIONAL_TOOLS[@]}"; do
    if command -v "$tool" &> /dev/null; then
        echo -e "  [OK] $tool"
    else
        echo -e "  \033[33m[WARN] $tool 未找到 (可选)\033[0m"
    fi
done

# 5. 总结
log_info "=== 验证结果 ==="
if [ $MISSING_TOOLS -eq 0 ]; then
    log_info "✅ 环境验证通过！所有核心工具均已安装且可运行。"
    log_info "您可以开始运行组装流程了: bash scripts/run_assembly.sh"
    exit 0
else
    log_error "❌ 环境验证失败！有 $MISSING_TOOLS 个核心工具缺失或无法运行。"
    log_error "请检查上述错误信息，并尝试重新运行 scripts/setup_tools.sh"
    exit 1
fi
