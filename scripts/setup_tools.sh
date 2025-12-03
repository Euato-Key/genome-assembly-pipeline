#!/bin/bash
# ==============================================================================
# 脚本名称: setup_tools.sh
# 功能: 自动化环境安装与配置脚本
# 特点: 
#   1. 使用 mamba 加速安装
#   2. 自动编译 hifiasm 以解决 CPU 指令集兼容性问题
#   3. 包含完整的环境验证步骤
# ==============================================================================

set -e  # 遇到错误立即退出

# ------------------------------------------------------------------------------
# 1. 初始化配置
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 2. 检查 Conda/Mamba 环境
# ------------------------------------------------------------------------------
check_conda() {
    if ! command -v conda &> /dev/null; then
        log_error "未检测到 Conda，请先安装 Miniconda 或 Anaconda"
        log_info "安装命令: wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh"
        exit 1
    fi
    log_info "Conda 已安装"
    
    # 强制更新 conda
    log_info "正在更新 conda 到最新版本..."
    conda update -n base -c defaults conda -y || log_warn "Conda 更新失败，将继续使用当前版本"
}

install_mamba() {
    log_info "检查 mamba..."
    if command -v mamba &> /dev/null; then
        log_info "mamba 已安装"
        PKG_MANAGER="mamba"
    else
        log_info "尝试安装 mamba 以加速包管理..."
        # 尝试安装 mamba，如果失败则降级使用 conda
        conda install -y -n base -c conda-forge mamba || {
            log_warn "mamba 安装失败，将使用 conda（速度较慢）"
            PKG_MANAGER="conda"
            return 0
        }
        PKG_MANAGER="mamba"
    fi
}

# ------------------------------------------------------------------------------
# 3. 创建基础环境
# ------------------------------------------------------------------------------
setup_env() {
    log_info "=== 开始配置环境 ==="
    
    # 配置镜像源（仅使用官方源，避免镜像同步问题）
    log_info "配置 conda 通道..."
    conda config --set channel_priority flexible
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
    # 删除旧环境
    if conda info --envs | grep -q "${CONDA_ENV_NAME}"; then
        log_info "删除旧环境 ${CONDA_ENV_NAME}..."
        conda env remove -n "${CONDA_ENV_NAME}" -y
    fi
    
    # 创建环境（仅安装核心工具，不含 hifiasm）
    log_info "创建环境并安装核心工具..."
    
    # 工具列表
    CORE_TOOLS=(
        "python=3.9"
        "fastp"
        "jellyfish"
        "samtools"
        "chromap"
        "yahs"
        "wget"
        "pigz"
        "make"           # 编译依赖
        "gcc_linux-64"   # 编译依赖
        "gxx_linux-64"   # 编译依赖
        "zlib"           # 编译依赖
        "git"            # 编译依赖
    )
    
    log_info "安装工具: ${CORE_TOOLS[*]}"
    $PKG_MANAGER create -n "${CONDA_ENV_NAME}" -y -c bioconda -c conda-forge "${CORE_TOOLS[@]}" || {
        log_error "环境创建失败"
        exit 1
    }
    
    # 激活环境
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV_NAME}"
    
    # 安装可选工具（失败不影响主流程）
    log_info "安装可选工具..."
    $PKG_MANAGER install -y -c bioconda genomescope2 || log_warn "genomescope2 安装失败（可选）"
    $PKG_MANAGER install -y -c bioconda kraken2 || log_warn "kraken2 安装失败（可选）"
    $PKG_MANAGER install -y -c bioconda busco || log_warn "busco 安装失败（可选）"
}

# ------------------------------------------------------------------------------
# 4. 编译安装 Hifiasm
# ------------------------------------------------------------------------------
install_hifiasm() {
    log_info "=== 开始编译安装 hifiasm ==="
    log_info "说明: 为了解决 CPU 指令集兼容性问题，我们将从源码编译 hifiasm"
    
    BUILD_DIR="${PROJECT_ROOT}/tools/hifiasm_build"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    # 1. 获取源码
    if [ ! -d "hifiasm" ]; then
        log_info "克隆 hifiasm 源码..."
        git clone https://github.com/chhylp123/hifiasm.git || {
            log_warn "Git 克隆失败，尝试下载 zip..."
            wget https://github.com/chhylp123/hifiasm/archive/refs/heads/master.zip -O hifiasm.zip
            unzip -q hifiasm.zip
            mv hifiasm-master hifiasm
        }
    fi
    
    cd hifiasm
    
    # 2. 编译
    log_info "执行编译 (make)..."
    make clean 2>/dev/null || true
    make -j$(nproc) || {
        log_error "编译失败，请检查编译环境"
        exit 1
    }
    
    # 3. 验证编译结果 (使用用户提供的测试用例)
    log_info "验证编译结果..."
    
    # 下载测试数据
    if [ ! -f "chr11-2M.fa.gz" ]; then
        log_info "下载测试数据 chr11-2M.fa.gz..."
        wget https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz || log_warn "测试数据下载失败，跳过验证"
    fi
    
    if [ -f "chr11-2M.fa.gz" ]; then
        log_info "运行测试组装..."
        ./hifiasm -o test -t4 -f0 chr11-2M.fa.gz 2> test.log || {
            log_error "hifiasm 测试运行失败！可能是 CPU 指令集仍然不兼容"
            exit 1
        }
        
        log_info "生成结果文件..."
        if [ -f "test.bp.p_ctg.gfa" ]; then
            awk '/^S/{print ">"$2;print $3}' test.bp.p_ctg.gfa > test.p_ctg.fa
            if [ -s "test.p_ctg.fa" ]; then
                log_info "✅ 验证成功！hifiasm 运行正常且生成了结果。"
            else
                log_warn "验证警告: 生成的 FASTA 文件为空"
            fi
        else
            log_error "验证失败: 未生成 GFA 文件"
            exit 1
        fi
    fi
    
    # 4. 安装到 Conda 环境
    CONDA_PREFIX=$(conda info --base)/envs/${CONDA_ENV_NAME}
    log_info "安装 hifiasm 到: ${CONDA_PREFIX}/bin/"
    cp hifiasm "${CONDA_PREFIX}/bin/"
    chmod +x "${CONDA_PREFIX}/bin/hifiasm"
    
    # 清理临时文件
    rm -f test* chr11-2M.fa.gz
    cd "${PROJECT_ROOT}"
}

# ------------------------------------------------------------------------------
# 5. 配置 Shell 集成
# ------------------------------------------------------------------------------
setup_shell() {
    log_info "=== 配置 Shell 集成 ==="
    
    # 检测当前 shell
    CURRENT_SHELL=$(basename "$SHELL")
    log_info "检测到 shell: $CURRENT_SHELL"
    
    # 初始化 conda
    if [ "$CURRENT_SHELL" = "bash" ] || [ "$CURRENT_SHELL" = "zsh" ]; then
        conda init "$CURRENT_SHELL" > /dev/null 2>&1 || log_info "conda init 可能已经执行过"
        log_info "已配置 $CURRENT_SHELL 的 conda 集成"
    else
        log_warn "未识别的 shell: $CURRENT_SHELL，请手动运行: conda init $CURRENT_SHELL"
    fi
}

# ------------------------------------------------------------------------------
# 主流程
# ------------------------------------------------------------------------------
main() {
    log_info "开始安装流程..."
    
    check_conda
    install_mamba
    setup_env
    install_hifiasm
    setup_shell
    
    log_info ""
    log_info "✅ 所有步骤完成！"
    log_info ""
    log_info "📋 下一步操作："
    log_info "1. 重新加载 shell 配置:"
    if [ "$CURRENT_SHELL" = "bash" ]; then
        log_info "   source ~/.bashrc"
    elif [ "$CURRENT_SHELL" = "zsh" ]; then
        log_info "   source ~/.zshrc"
    fi
    log_info "2. 激活环境:"
    log_info "   conda activate ${CONDA_ENV_NAME}"
    log_info "3. 验证安装:"
    log_info "   hifiasm --version"
}

main "$@"
