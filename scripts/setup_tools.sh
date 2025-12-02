#!/bin/bash
# ==============================================================================
# 脚本名称: setup_tools.sh
# 功能: 安装环境、工具及下载测试数据
# ==============================================================================

set -e # 遇到错误立即退出

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

# ------------------------------------------------------------------------------
# 1. 环境检查与 Conda/Mamba 安装
# ------------------------------------------------------------------------------
install_conda_and_mamba() {
    # 检查是否已安装 conda
    if command -v conda &> /dev/null; then
        log_info "Conda 已安装。"
    else
        log_info "未检测到 Conda，开始安装 Miniconda..."
        mkdir -p "${PROJECT_ROOT}/tools"
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "${PROJECT_ROOT}/tools/miniconda.sh"
        bash "${PROJECT_ROOT}/tools/miniconda.sh" -b -p "${PROJECT_ROOT}/tools/miniconda3"
        rm "${PROJECT_ROOT}/tools/miniconda.sh"
        
        # 初始化 conda
        source "${PROJECT_ROOT}/tools/miniconda3/bin/activate"
        conda init bash
        log_info "Miniconda 安装完成。"
    fi

    # 检查是否已安装 mamba
    if command -v mamba &> /dev/null; then
        log_info "Mamba 已安装。"
    else
        log_info "安装 Mamba 以加速包管理..."
        # 配置镜像源
        configure_conda_mirrors
        # 安装 mamba
        conda install mamba -c conda-forge -y
        log_info "Mamba 安装完成。"
    fi
}

# ------------------------------------------------------------------------------
# 配置镜像源
# ------------------------------------------------------------------------------
configure_conda_mirrors() {
    log_info "正在配置 Conda/Mamba 镜像源 (清华源)..."
    
    # 清理现有配置
    conda config --remove-key channels 2>/dev/null || true
    
    # 添加镜像源
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
    conda config --set show_channel_urls yes
    
    # 设置镜像源优先级
    conda config --set channel_priority flexible
}

# ------------------------------------------------------------------------------
# 2. 创建环境并安装工具 (使用 Mamba)
# ------------------------------------------------------------------------------
setup_env() {
    log_info "正在配置 Conda 环境: ${CONDA_ENV_NAME}..."
    
    # 配置镜像源
    configure_conda_mirrors

    # 检查环境是否存在
    if conda info --envs | grep -q "${CONDA_ENV_NAME}"; then
        log_info "检测到环境 ${CONDA_ENV_NAME} 已存在。"
        read -p "是否重新创建环境？(y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_info "正在删除旧环境..."
            conda env remove -n "${CONDA_ENV_NAME}" -y
        else
            log_info "使用现有环境。"
            return 0
        fi
    fi

    # 使用 Mamba 创建环境并安装软件 (速度更快)
    log_info "使用 Mamba 创建环境并安装工具 (这比 conda 快得多)..."
    
    mamba create -n "${CONDA_ENV_NAME}" -y \
        fastp \
        jellyfish \
        genomescope2 \
        kraken2 \
        hifiasm \
        flye \
        busco \
        chromap \
        yahs \
        samtools \
        openjdk=11 \
        parallel \
        r-base \
        wget \
        pigz
    
    # 如果 mamba 创建失败，回退到 conda
    if [ $? -ne 0 ]; then
        log_error "Mamba 创建环境失败，尝试使用 Conda..."
        conda create -n "${CONDA_ENV_NAME}" -y \
            fastp \
            jellyfish \
            genomescope2 \
            kraken2 \
            hifiasm \
            flye \
            busco \
            chromap \
            yahs \
            samtools \
            openjdk=11 \
            parallel \
            r-base \
            wget \
            pigz
    fi
    
    log_info "Conda 环境创建完成。"

    # 下载 Juicer Tools
    if [ ! -f "${JUICER_JAR}" ]; then
        log_info "正在下载 Juicer Tools..."
        mkdir -p "$(dirname "${JUICER_JAR}")"
        # 尝试多个镜像源
        juicer_urls=(
            "https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar"
            "https://github.com/aidenlab/juicer/releases/download/1.22.01/juicer_tools_1.22.01.jar"
        )
        
        for url in "${juicer_urls[@]}"; do
            log_info "尝试从 $url 下载..."
            if wget --timeout=30 -t 3 -O "${JUICER_JAR}.tmp" "$url"; then
                mv "${JUICER_JAR}.tmp" "${JUICER_JAR}"
                log_info "Juicer Tools 下载完成。"
                break
            else
                log_error "下载失败: $url"
            fi
        done
        
        if [ ! -f "${JUICER_JAR}" ]; then
            log_error "所有 Juicer Tools 下载源都失败，请手动下载。"
        fi
    fi
}

# ------------------------------------------------------------------------------
# 3. 验证安装
# ------------------------------------------------------------------------------
verify_install() {
    log_info "正在验证工具安装..."
    
    # 激活环境
    if [ -f "${PROJECT_ROOT}/tools/miniconda3/etc/profile.d/conda.sh" ]; then
        source "${PROJECT_ROOT}/tools/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
    fi
    
    # 尝试激活环境
    if conda activate "${CONDA_ENV_NAME}" 2>/dev/null; then
        log_info "成功激活环境 ${CONDA_ENV_NAME}"
    else
        log_error "无法激活环境 ${CONDA_ENV_NAME}"
        exit 1
    fi

    local tools=(fastp jellyfish kraken2 hifiasm flye busco chromap yahs samtools java parallel)
    local fail=0

    for tool in "${tools[@]}"; do
        if command -v "$tool" &> /dev/null; then
            echo -e "  [OK] $tool"
        else
            echo -e "  \033[31m[FAIL] $tool 未找到\033[0m"
            fail=1
        fi
    done

    if [ -f "${JUICER_JAR}" ]; then
        echo -e "  [OK] Juicer Tools Jar"
    else
        echo -e "  \033[31m[FAIL] Juicer Tools Jar 未找到\033[0m"
        fail=1
    fi

    if [ $fail -eq 0 ]; then
        log_info "所有工具验证通过！"
    else
        log_error "部分工具未安装成功，请检查日志。"
        exit 1
    fi
}

# ------------------------------------------------------------------------------
# 4. 下载测试数据
# ------------------------------------------------------------------------------
download_test_data() {
    log_info "开始下载测试数据..."
    mkdir -p "${PROJECT_ROOT}/test_data"
    cd "${PROJECT_ROOT}/test_data"

    # 这里我们下载一个非常小的测试集，例如 酵母 (S. cerevisiae) 或者 噬菌体
    # 为了演示，我们下载 E. coli 的公开数据 (体积较小)
    
    # 1. HiFi 数据
    if [ ! -f "HiFi.fasta" ]; then
        log_info "下载 HiFi 测试数据..."
        # 这是一个示例链接，实际使用时请替换为真实可用的测试数据链接
        # 这里使用 wget 模拟下载，实际上如果链接失效需要用户手动提供
        echo "正在模拟下载 HiFi 数据..."
        # wget -O HiFi.fasta https://example.com/ecoli_hifi.fasta
        touch HiFi.fasta # 占位
        echo ">seq1" > HiFi.fasta
        echo "ATGC" >> HiFi.fasta
        log_info "HiFi 测试数据下载完成 (模拟)。"
    fi

    # 2. Hi-C 数据
    if [ ! -f "Hi-C_R1.fastq" ]; then
        log_info "下载 Hi-C 测试数据..."
        touch Hi-C_R1.fastq Hi-C_R2.fastq # 占位
        log_info "Hi-C 测试数据下载完成 (模拟)。"
    fi

    # 3. 数据库 (Busco / Kraken)
    # 数据库通常很大，这里仅创建目录结构并提示
    mkdir -p "${PROJECT_ROOT}/database/kraken_db"
    mkdir -p "${PROJECT_ROOT}/database/busco_db"
    
    log_info "测试数据目录结构已建立: ${PROJECT_ROOT}/test_data"
    log_info "注意：真实的基因组数据库体积巨大 (几十GB)，脚本中未自动下载。"
    log_info "请参考 docs/data_download.md 手动下载完整数据库。"
}

# ------------------------------------------------------------------------------
# 主流程
# ------------------------------------------------------------------------------
main() {
    if [[ "$1" == "--download-test-data" ]]; then
        download_test_data
        exit 0
    fi

    install_conda_and_mamba
    setup_env
    verify_install
    
    log_info "环境设置完成！使用以下命令激活环境："
    log_info "conda activate ${CONDA_ENV_NAME}"
}

main "$@"